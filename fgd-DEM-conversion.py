

"""
***************************************************************************
* fgd-DEM-conversion.py: convert DEMs in FGD-GML format to GeoTiff        *
* AUTHOR(S): Daisuke Takahashi <dtakahshi42@gmail.com>                    *
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License version 3.       *
*                                                                         *
***************************************************************************
"""

from io import StringIO
from warnings import warn
from pathlib import Path
from itertools import chain
from zipfile import ZipFile
from argparse import ArgumentParser

from xml.etree.ElementTree import parse as parseXML
from pandas import read_csv, Series
from pandas.api.types import CategoricalDtype
import numpy as np

try:
    from joblib import Parallel, delayed
except ImportError:
    # use a fake version instead. no parallel prcessings
    def iterate(gen):
        for x in gen:
            pass
    def Parallel(n_jobs = 1):
        return iterate
    def delayed(f):
        return f
                
from osgeo import gdal
from osgeo import osr

# reference:
# 基盤地図情報ダウンロードデータ ファイル仕様書 第4.1版
# FGD_GMLSchema.xsd v4.1

def read_as_num(txt):
    return np.array([float(x) for x in txt.strip().split(' ') if x != ''])

class BBox:
    def __repr__(self):
        return ('BBox(%s %s)' % (self.low, self.high))
    
    def __init__(self, low, high):
        self.high = np.array(high)
        self.low = np.array(low)
        return
    
    def merge(self, bb):
        return BBox(np.minimum(self.low, bb.low),
                    np.maximum(self.high, bb.high))
    
    def size(self):
        return self.high - self.low

# DEM構成点種別（列挙型）
# DEM構成点の種別。(p.18)
# 種別が「データなし」の場合，この属性値には ”-9999.” が設定される。(p.18)
DEM_point_category = CategoricalDtype([
    'データなし',  # 0
    '地表面',     # 1
    '表層面',     # 2
    '海水面',     # 3
    '内水面',     # 4
    'その他'      # 5
])
NODATA_code = Series(['データなし'], dtype=DEM_point_category).cat.codes[0]

ns = {
    'fgd': 'http://fgd.gsi.go.jp/spec/2008/FGD_GMLSchema',
    'gml': 'http://www.opengis.net/gml/3.2'
}

cell_type = {
    'Cell_type': {
        'np': np.byte,
        'gdal': gdal.GDT_Byte
    },
    'Altitude': {
        'np': np.float32,
        'gdal': gdal.GDT_Float32
    }
}

class demPatch:
    def __init__(self, et):
        # Is this "coverage" unique for each file?
        # Maybe the definition is,
        # coverage [1..*] : CV_DiscreteGridPointCoverage
        # このDEM区画を構成する被覆情報。
        # CV_DiscreteGridPointCoverageは， JIS X7123で定義された離散型グリッド点被覆を
        # 表現するためのクラスである。(p.16)
        # allowing multipart files...
        # I assume its uniqueness, and give a warning if it is not unique.
        coverages = et.findall('./fgd:DEM/fgd:coverage', ns)
        if len(coverages) != 1:
            warn('multiple coverages found, but only the first occurence is used')
        coverage = coverages[0]
        georange = coverage.find('./gml:boundedBy/gml:Envelope', ns)
        self.extent = BBox(read_as_num(georange.find('./gml:lowerCorner', ns).text),
                           read_as_num(georange.find('./gml:upperCorner', ns).text))
        d_index = coverage.find('./gml:gridDomain/gml:Grid/gml:limits/gml:GridEnvelope',
                                ns)
        # axisNames : Sequence<CharacterString>
        # グリッドセルの座標軸の名称。
        # 基盤地図情報では，”x y” と定義する。x軸は経度の正方向，y軸は緯度の正方向を意
        # 味している。
        # セル番号 (m, n) で，mはx軸の値，nはy軸の値である。(p.17)
        #
        # 表2-1 DEM種別ごとのextent設定値 (p.17)
        # DEM種別                 low属性値 high属性値
        # 5m メッシュ(数値地形) :  [0, 0]   [224, 149]
        # 5m メッシュ(標高) :      [0, 0]   [224, 149]
        # 10m メッシュ(標高) :     [0, 0]   [1124, 749]
        # 10m メッシュ(火山標高) : [0, 0]   [1124, 749]
        indexRange = np.array([read_as_num(d_index.find('./gml:low', ns).text),
                               read_as_num(d_index.find('./gml:high', ns).text)]).astype(int)
        # +1 to include endpoints
        shape = tuple(indexRange[1,:] - indexRange[0,:] + 1)
        cell_count = np.prod(shape)
        #
        # The default value is 'nodata'
        self.state_array = np.full(cell_count, NODATA_code,
                                   dtype=cell_type['Cell_type']['np'])
        self.cell_array = np.full(cell_count, -9999.0,
                                  dtype=cell_type['Altitude']['np'])
        #
        # なお，先頭部分で連続した構成点の値が存在しない場合は，valuesにおける値の指定
        # を省略することができる。その場合，実際に構成点の値が開始するグリッドセルを
        # startSequenceで指定する。(p.18) <-- startPoint?
        startPoint = read_as_num(
            coverage.find('./gml:coverageFunction/gml:GridFunction/gml:startPoint',
                    ns).text
        ).astype(int)
        npaddings = startPoint[0] + startPoint[1]*shape[0]
        #
        # Read a data table:
        dem_data_path = './gml:rangeSet/gml:DataBlock/gml:tupleList'
        tbl = read_csv(StringIO(coverage.find(dem_data_path, ns).text),
                       header=None,
                       dtype={
                           'Cell_type': DEM_point_category,
                           'Altitude': cell_type['Altitude']['np']
                       },
                       names=['Cell_type', 'Altitude'])
        #
        catcodes = tbl['Cell_type'].cat.codes.astype(cell_type['Cell_type']['np'])
        self.state_array[npaddings:(npaddings + tbl.shape[0])] = catcodes
        self.cell_array[npaddings:(npaddings + tbl.shape[0])]  = tbl['Altitude']
        #
        # また，末尾部分で連続した構成点の値が存在しない場合は，valuesにおける値の指定
        # を省略することができる。valueの配列で設定された値の数がグリッドセルの末尾に
        # 達しない場合，その後ろは省略されている。(p.18)
        #
        # 基盤地図情報では，type属性値=”Linear”，scanDirection属性値=”+x –y” と設定する。
        # この設定値は，先頭セルは北西端にあって，配列順序がx軸の正方向（西→東の順）
        # へ順に並んでおり，東端に達すると次に，y軸の負方向（北→南の順）に進む方式で
        # 南東端に至る配列であることを示している。(p.18)
        #
        # numpy array is row major
        self.cell_array.shape = (shape[1], shape[0])
        self.state_array.shape = (shape[1], shape[0])
        return

class DEMRasterizer:
    def __init__(self, zip_file_name):
        self.zipfile_path = Path(zip_file_name)
        #
        zp = ZipFile(str(self.zipfile_path))
        zp_filelist = [x.filename for x in zp.infolist()]
        patches = [demPatch(parseXML(zp.open(x))) for x in zp_filelist]
        ext = patches[0].extent
        for x in patches:
            ext = ext.merge(x.extent)
        self.extent = ext
        #
        self.lat = [ext.low[0], ext.high[0]]
        self.lon = [ext.low[1], ext.high[1]]
        #
        # Because text-serialized float numbers are almost always different from
        # original numbers, we need to accept some numerical inaccuracies.
        # I hope that sizes of these patches do not vary so much.
        average_bbox_size = sum([x.extent.size() for x in patches])/len(patches)
        #
        grid_size = np.round(self.extent.size()/average_bbox_size).astype(int)
        patch_locations = [(grid_size[0] - int(x[0]) - 1, int(x[1])) # (N-S, E-W)
                           for x in [np.round((x.extent.low - self.extent.low)/average_bbox_size)
                                     for x in patches]]
        #
        pshape = patches[0].cell_array.shape
        image_dimension = (pshape*grid_size).astype(int)
        #
        raster = np.full(tuple(image_dimension), -9999.0,
                         dtype=cell_type['Altitude']['np'])
        gtype = np.full(tuple(image_dimension), NODATA_code,
                        dtype=cell_type['Cell_type']['np'])
        for ((loc_y, loc_x), p) in zip(patch_locations, patches):
            start_y = loc_y*pshape[0]   # N
            start_x = loc_x*pshape[1]   # W
            end_y = start_y + pshape[0] # S (not included)
            end_x = start_x + pshape[1] # E (not included)
            # Numpy stores data in row major order. The first index is a row index.
            raster[start_y:end_y,:][:,start_x:end_x] = p.cell_array
            gtype[start_y:end_y,:][:,start_x:end_x] = p.state_array
        #
        self.raster = {
            'Altitude': raster,
            'Cell_type': gtype
        }
        return
    
    def write_raster(self, raster_type='Altitude', dest_dir='.'):
        fpath = (Path(dest_dir)
                 / (Path(self.zipfile_path.name).stem + '_' + raster_type)).with_suffix('.tif')
        #
        raster = self.raster[raster_type]
        #
        # DEM区画は，標高を測量し，又は算定した地点の集合体であり，経緯度によって四角く分割された１メッシュに
        # おける標高値の分布を数値標高モデルとして表現するためのクラスである。その１メッシュ内をグリッドによっ
        # て格子状に分割したセルに対して，標高値が割り当てられる。ここで表現される数値標高モデルは，JIS X7123
        # に準拠した被覆としての形式・構成となっている。
        # DEM区画は，経緯度により四角く分割された１メッシュを表し，そのメッシュ内がさらにグリッドにより格子状に
        # 分割されたセルに標高値が割り当てられる。(p.15)
        # --> an "extent" in the dataset seems to be an extent of rectangular grid cells (not points)
        (xmin, ymin, xmax, ymax) = [min(self.lon), min(self.lat),
                                    max(self.lon), max(self.lat)]
        dx = (xmax - xmin)/raster.shape[1]
        dy = (ymax - ymin)/raster.shape[0]
        geotransform = (xmin, dx, 0, ymax, 0, -dy)
        destination = gdal.GetDriverByName('GTiff').Create(fpath.name,
                                                           raster.shape[1],
                                                           raster.shape[0],
                                                           1,
                                                           cell_type[raster_type]['gdal'])
        destination.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        # EPSG6668: JGD2011
        res = srs.ImportFromEPSG(6668)
        if res != 0:
            raise RuntimeError(repr(res) + ': could not import from EPSG')
        destination.SetProjection(srs.ExportToWkt())
        destination.GetRasterBand(1).WriteArray(raster)
        if raster_type == 'Altitude':
            destination.GetRasterBand(1).SetNoDataValue(-9999.0)
        destination.FlushCache()
        destination = None
        return

def convert_dem (path_string, destdir, write_cell_type):
    r = DEMRasterizer(path_string)
    r.write_raster('Altitude', destdir)
    if write_cell_type:
        r.write_raster('Cell_type', destdir)
    return

argparser = ArgumentParser(description = '''
Covert DEMs published by Geospatial Information Authority of Japan to GeoTiffs
''')
argparser.add_argument('targets', metavar='target', type=str, nargs='*',
                       help='a target zip file or a directory containing zip files')
argparser.add_argument('--dest-dir', dest='dest', type=str,
                       default='.',
                       help='a destination directory to write tiff files')
argparser.add_argument('--with-cell-type', dest='with_type',
                       action='store_true',
                       default=False,
                       help='also write types of DEM cells (as separate GeoTiffs)')
argparser.add_argument('--nproc', type=int,
                       default=-1,
                       help='the number of processes (default: same as the number of CPU cores)')

if __name__ == '__main__':
    args = argparser.parse_args()
    #
    target_paths = [Path(x) for x in args.targets]
    targets_as_files = (d for d in target_paths if d.is_file())
    targets_as_dirs  = [Path(d).glob('*.zip') for d in target_paths if d.is_dir()]
    targets = chain(targets_as_files, *targets_as_dirs)
    #
    Parallel(n_jobs = args.nproc)(
        delayed(convert_dem)(p, args.dest, args.with_type) for p in targets
    )
    exit

