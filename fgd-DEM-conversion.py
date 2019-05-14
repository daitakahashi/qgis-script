
import io
from pathlib import Path
import zipfile

import lxml.etree as ET
import pandas as pd
import numpy as np
from joblib import Parallel, delayed

from osgeo import gdal
from osgeo import osr

# v.4.0のドキュメントとv4.1のschemaを参照した

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
def conv_type(s):
    return {
        '地表面': 0,
        '表層面': 1,
        '海水面': 2,
        '内水面': 3,
        'その他': 4,
        'データなし': -1
    }[s]

ns = {
    None: 'http://fgd.gsi.go.jp/spec/2008/FGD_GMLSchema',
    'gml': 'http://www.opengis.net/gml/3.2'
}

cell_type = {
    'Status': {
        'np': np.byte,
        'gdal': gdal.GDT_Byte
    },
    'Altitude': {
        'np': np.float32,
        'gdal': gdal.GDT_Float32
    }
}

class demPatch5m:
    def __init__(self, et):
        georange = et.find('./DEM/coverage/gml:boundedBy/gml:Envelope', ns)
        self.extent = BBox(read_as_num(georange.find('./gml:lowerCorner', ns).text),
                           read_as_num(georange.find('./gml:upperCorner', ns).text))
        d_index = et.find('./DEM/coverage/gml:gridDomain/gml:Grid/gml:limits/gml:GridEnvelope', ns)
        # axisNames : Sequence<CharacterString>
        # グリッドセルの座標軸の名称。
        # 基盤地図情報では，”x y” と定義する。x軸は経度の正方向，y軸は緯度の正方向を意
        # 味している。
        # セル番号 (m, n) で，mはx軸の値，nはy軸の値である。(p.17)
        #
        # low属性値: [0, 0], high属性値: [224, 149] (p.17)
        self.indexRange = np.array([read_as_num(d_index.find('./gml:low', ns).text),
                                    read_as_num(d_index.find('./gml:high', ns).text)]).astype(int)
        shape = tuple(self.indexRange[1,:] - self.indexRange[0,:] + 1)
        cell_count = np.prod(shape)
        if shape != (225, 150):
            raise ValueError('Shape of an image should be (225, 150), but got %s' % shape)
        #
        # なお，先頭部分で連続した構成点の値が存在しない場合は，valuesにおける値の指定
        # を省略することができる。その場合，実際に構成点の値が開始するグリッドセルを
        # startSequenceで指定する。(p.18) <-- startPoint?
        startPoint = read_as_num(
            et.find('./DEM/coverage/gml:coverageFunction/gml:GridFunction/gml:startPoint', ns).text
        ).astype(int)
        npaddings = startPoint[0] + startPoint[1]*shape[0]
        #
        # Read data table:
        dem_data_path = './DEM/coverage/gml:rangeSet/gml:DataBlock/gml:tupleList'
        tbl = pd.read_csv(io.StringIO(et.find(dem_data_path, ns).text),
                          header=None,
                          dtype={'Altitude': cell_type['Altitude']['np']},
                          converters={'Status': conv_type},
                          names=['Status', 'Altitude'])
        #
        # また，末尾部分で連続した構成点の値が存在しない場合は，valuesにおける値の指定
        # を省略することができる。valueの配列で設定された値の数がグリッドセルの末尾に
        # 達しない場合，その後ろは省略されている。(p.18)
        npaddings_end = cell_count - tbl.shape[0] - npaddings
        #
        # 基盤地図情報では，type属性値=”Linear”，scanDirection属性値=”+x –y” と設定する。
        # この設定値は，先頭セルは北西端にあって，配列順序がx軸の正方向（西→東の順）
        # へ順に並んでおり，東端に達すると次に，y軸の負方向（北→南の順）に進む方式で
        # 南東端に至る配列であることを示している。(p.18)
        #
        # The default value is 'nodata'
        self.state_array = np.full(cell_count, conv_type('データなし'),
                                   dtype=cell_type['Status']['np'])
        self.cell_array = np.full(cell_count, -9999.0,
                                  dtype=cell_type['Altitude']['np'])
        #
        self.state_array[npaddings:(npaddings + tbl.shape[0])] = tbl['Status']
        self.cell_array[npaddings:(npaddings + tbl.shape[0])]  = tbl['Altitude']
        #
        self.cell_array.shape = (shape[1], shape[0])
        self.state_array.shape = (shape[1], shape[0])
        return

class DEMRasterizer:
    def __init__(self, zip_file_name):
        self.zipfile_path = Path(zip_file_name)
        #
        zp = zipfile.ZipFile(str(self.zipfile_path))
        zp_filelist = [x.filename for x in zp.infolist()]
        self.patches = [demPatch5m(ET.parse(zp.open(x))) for x in zp_filelist]
        ext = self.patches[0].extent
        for x in self.patches:
            ext = ext.merge(x.extent)
        self.extent = ext
        #
        self.lat = [ext.low[0], ext.high[0]]
        self.lon = [ext.low[1], ext.high[1]]
        #
        average_bbox_size = sum([x.extent.size() for x in self.patches])/len(self.patches)
        #
        # !! [セル数: N-S, セル数: W-E] !!
        self.grid_size = np.round(self.extent.size()/average_bbox_size).astype(int)
        self.patch_locations = [(self.grid_size[0] - int(x[0]) - 1, int(x[1])) # (N-S, E-W)
                                for x in [np.round((x.extent.low - self.extent.low)/average_bbox_size)
                                          for x in self.patches]]
        self.image_dimension = (np.array([
            150, # N-S (y)
            225  # E-W (x)
        ])*self.grid_size).astype(int)
        #
        raster = np.full(tuple(self.image_dimension), -9999.0,
                              dtype=cell_type['Altitude']['np'])
        for ((loc_y, loc_x), p) in zip(self.patch_locations, self.patches):
            start_ix_y = loc_y*150 # N-S
            start_ix_x = loc_x*225 # E-W
            # Numpy stores data in row major order. And the first index is a row index.
            raster[start_ix_y:(start_ix_y + 150),:][:,start_ix_x:(start_ix_x + 225)] =  p.cell_array
        #
        status = np.full(tuple(self.image_dimension), conv_type('データなし'),
                             dtype=cell_type['Status']['np'])
        for ((loc_y, loc_x), p) in zip(self.patch_locations, self.patches):
            start_ix_y = loc_y*150 # N-S
            start_ix_x = loc_x*225 # E-W
            # Numpy stores data in row major order. And the first index is a row index.
            status[start_ix_y:(start_ix_y + 150),:][:,start_ix_x:(start_ix_x + 225)] =  p.state_array
        #
        self.raster = {
            'Altitude': raster,
            'Status': status
        }
        return
    
    def write_raster(self, raster_type='Altitude', dest_dir='.'):
        fpath = (Path(dest_dir)
                 / (Path(self.zipfile_path.name).stem + '_' + raster_type)).with_suffix('.tif')
        #
        raster = self.raster[raster_type]
        #
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
        res = srs.ImportFromEPSG(4326)
        if res != 0:
            raise RuntimeError(repr(res) + ': could not import from EPSG')
        # Qgis complains that there are no crs information... why?
        # Interestingly, there are no such problems on linux,
        # so, I guess it is (1) an OSX issue, (2) a multibyte-pathname issue,
        # or (3) something other problems in this code.
        destination.SetProjection(srs.ExportToWkt())
        destination.GetRasterBand(1).WriteArray(raster)
        destination.GetRasterBand(1).SetNoDataValue(-9999.0)
        destination.FlushCache()
        destination = None
        return

def convert_dem (path):
    r = DEMRasterizer(str(path))
    r.write_raster('Altitude')
    r.write_raster('Status')
    return

# Convert all zip files in a directory 'PackDLMap'
Parallel(n_jobs = 8)(
    delayed(convert_dem)(p)
    for p in Path('PackDLMap').glob('*.zip')
)
