# -*- coding: utf-8 -*-

"""
***************************************************************************
*   import-fgd-DEM-tile.py:                                               *
*     Download and merge DEM tiles from cyberjapandata.gsi.go.jp          *
*   AUTHOR(S): Daisuke Takahashi <dtakahshi42@gmail.com>                  *
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License version 3.       *
*                                                                         *
***************************************************************************
"""

import pathlib
import shutil
import urllib.request
import math
import itertools

from PyQt5.QtCore import QCoreApplication

from qgis.core import (QgsProcessing,
                       QgsProcessingException,
                       QgsProcessingUtils,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterExtent,
                       QgsProcessingParameterEnum,
                       QgsProcessingParameterBoolean,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterRasterDestination,
                       QgsCoordinateReferenceSystem)
import processing

from osgeo import osr, gdal

from PIL import Image
import numpy as np

class importFGDDEMTiles(QgsProcessingAlgorithm):

    INPUT        = 'INPUT'
    DEMTYPE      = 'DEMTYPE'
    CLEARCACHE   = 'CLEARCACHE'
    NODATA_VALUE = 'NODATA_VALUE'
    OUTPUT       = 'OUTPUT'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return importFGDDEMTiles()

    def name(self):
        return 'importfgddemtiles'

    def displayName(self):
        return self.tr('Import a DEM from 地理院タイル')

    def group(self):
        return self.tr('Raster')

    def groupId(self):
        return 'raster'

    def shortHelpString(self):
        return self.tr('''
国土地理院標高タイルから標高タイルを取り込む非公式ツール

このデータセットは国土地理院より地理院タイルの一つとして提供されており、使用の際には出典（「国土地理院」または「地理院タイル」等）の明記が必要です。詳細については https://maps.gsi.go.jp/development/ichiran.html (地理院タイル一覧) 及び http://www.gsi.go.jp/kikakuchousei/kikakuchousei40182.html (国土地理院コンテンツ利用規約) を参照ください。

The datasets are provided by the Geospatial Information Authority of Japan (国土地理院). When using their information, the Geospatial Information Authority of Japan asks to credit "国土地理院" or "地理院タイル" as an information source.

For the data sourses and further details, please visit https://maps.gsi.go.jp/development/ichiran.html and http://www.gsi.go.jp/kikakuchousei/kikakuchousei40182.html.
''')

    def initAlgorithm(self, config=None):
        
        self.addParameter(
            QgsProcessingParameterExtent(
                self.INPUT,
                self.tr('DEM extent')
            )
        )
        
        self.addParameter(
            QgsProcessingParameterEnum(
                self.DEMTYPE,
                self.tr('DEM type'),
                ['DEM 5A', 'DEM 5B', 'DEM 10B'],
                defaultValue=0
            )
        )
        
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.CLEARCACHE,
                self.tr('Clear cache files'),
                defaultValue=False
            )
        )
        
        self.addParameter(
            QgsProcessingParameterNumber(
                self.NODATA_VALUE,
                self.tr('Nodata value'),
                defaultValue=-9999.0
            )
        )
        
        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT,
                self.tr('Output DEM')
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        #
        # Load input parameters and choose temporary directory to write processed tiles
        #
        
        # Let (wrongly) say, the extent is rectangle even after the reprojection
        extent = self.parameterAsExtent(
            parameters,
            self.INPUT,
            context,
            QgsCoordinateReferenceSystem('EPSG:3857')
        )

        clear_cache_directory = self.parameterAsBool(
            parameters,
            self.CLEARCACHE,
            context
        )
        
        tilecache_directory = pathlib.Path(QgsProcessingUtils.tempFolder()) / self.name()
        if not tilecache_directory.exists():
            tilecache_directory.mkdir()
        elif clear_cache_directory:
            # Recreate cache directory if requested
            shutil.rmtree(tilecache_directory)
            tilecache_directory.mkdir()
        
        #
        # Setup constants
        #
        # Tile URLs are listed at here.
        # https://maps.gsi.go.jp/development/ichiran.html#dem
        baseurl = 'https://cyberjapandata.gsi.go.jp/xyz'

        demtype_index = self.parameterAsEnum(
            parameters,
            self.DEMTYPE,
            context
        )
        # See below URL for detail of DEMs.
        # https://fgd.gsi.go.jp/download/ref_dem.html
        if demtype_index == 0:
            # Use DEM 5A: 5mメッシュ (標高)
            demtype = 'dem5a_png'
            zoomlevel = 15
        elif demtype_index == 1:
            # Use DEM 5B: 5mメッシュ (数値地形)
            demtype = 'dem5b_png'
            zoomlevel = 15
        elif demtype_index == 2:
            # Use DEM 10B: 10mメッシュ (標高)
            demtype = 'dem_png'
            zoomlevel = 14
        else:
            raise QgsProcessingException('Unimplemented DEM type')
        
        # Choose a value representing nodata. Because a geotiff stores its nodata value
        # as a text-serialized number, use of non-integer values will cause roundoff
        # errors. Also, some algorithms (e.g., merge) seem to assume that the nodata
        # value is a proper number; therefore, NaN, -inf, or other special values may
        # not work well with those algorithms despite their logical correspondences
        # with 'no data'.
        # The default value is -9999.0, which is the same value as the original dataset
        # '基盤地図情報数値標高モデル'.
        nodata_value = self.parameterAsDouble(
            parameters,
            self.NODATA_VALUE,
            context
        )
        
        # PNG形式の場合
        # 24ビットカラーのPNG形式で、一つのタイルの大きさは256ピクセル×256ピクセルです。
        # http://maps.gsi.go.jp/development/demtile.html
        tile_dim = (256, 256)
        
        L = 4.00750166855784*(10**7)
        tile_size = L/(2**zoomlevel)
        dx = tile_size/tile_dim[0]
        dy = tile_size/tile_dim[1]
        
        #
        # Main processing function:
        #  Download a png tile and convert it to a geotiff DEM tile.
        #
        def get_DEM(index_x, index_y, zoomlevel, output_path):
            url = '{}/{}/{}/{}/{}.png'.format(baseurl, demtype, zoomlevel,
                                              index_x, index_y)
            
            tile_directory = output_path / demtype / str(zoomlevel) / str(index_x)
            if not tile_directory.exists():
                tile_directory.mkdir(parents=True)
            elif not tile_directory.is_dir():
                feedback.reportError(
                    'Cache inconsistency: {} exists but not a directory'.format(
                        tile_directory
                    ))
                tile_directory.unlink()
                tile_directory.mkdir()
            
            filepath = tile_directory / '{}.tif'.format(index_y)
            if filepath.exists():
                if filepath.is_file():
                    # found a processed image in cache files
                    return filepath
                else:
                    feedback.reportError(
                        'Cache inconsistency: {} exists but not a regular file'.format(
                            filepath
                        ))
                    if filepath.is_dir():
                        shutil.rmtree(filepath)
                    else:
                        filepath.unlink()
            
            # 画素値（RGB値）から標高値h（m）の計算式は下記のとおりです。
            #     x = 2^16*R + 2^8*G + B
            #     x < 2^23の場合　h = xu
            #     x = 2^23の場合　h = NA
            #     x > 2^23の場合　h = (x-2^24)u
            # uは標高分解能（0.01m）を表します。
            # http://maps.gsi.go.jp/development/demtile.html
            u = 0.01
            m = np.zeros(tile_dim, dtype=np.int32)
            
            tiff_options = []
            try:
                with urllib.request.urlopen(url) as img_data:
                    with Image.open(img_data) as img:
                        img_array = np.array(img)
                        m += img_array[:,:,0]
                        m *= 2**8
                        m += img_array[:,:,1]
                        m *= 2**8
                        m += img_array[:,:,2]
                
            except urllib.error.URLError as e:
                if e.code == 404:
                    # Quite possibly, just the tile doesn't exist (e.g., there are no ocean tiles),
                    # so generate an all-null tile as a place holder.
                    m.fill(2**23)
                    # Compression will work very well in this case
                    tiff_options.append('COMPRESS=LZW')
                else:
                    feedback.pushInfo('{}'.format(e))
                    return None
            
            m[m > 2**23] -= 2**24
            #
            H = m.astype(np.float32)*u
            H[m == 2**23] = np.nan
            #
            xmin = index_x*tile_size - L/2
            ymax = L/2 - index_y*tile_size
            dest = gdal.GetDriverByName('GTiff').Create(str(filepath),
                                                        256, 256,
                                                        1,
                                                        gdal.GDT_Float32,
                                                        options=tiff_options)
            dest.SetGeoTransform((xmin, dx, 0, ymax, 0, -dy))
            webMercator = osr.SpatialReference()
            webMercator.ImportFromEPSG(3857)
            dest.SetProjection(webMercator.ExportToWkt())
            dest.GetRasterBand(1).WriteArray(H)
            dest.FlushCache()
            dest = None
            #
            return filepath

        #
        # Actual work starts from here
        #
        xmax = extent.xMaximum()
        xmin = extent.xMinimum()
        ymax = extent.yMaximum()
        ymin = extent.yMinimum()
        margin = 0.00001

        index_x_lower = 2**(zoomlevel - 1) + int(xmin/tile_size + margin)
        index_x_upper = 2**(zoomlevel - 1) + int(math.ceil(xmax/tile_size - margin))
        index_y_lower = 2**(zoomlevel - 1) - int(math.ceil(ymax/tile_size - margin))
        index_y_upper = 2**(zoomlevel - 1) - int(ymin/tile_size + margin)
        
        index_x_range = range(index_x_lower, index_x_upper)
        index_y_range = range(index_y_lower, index_y_upper)
        
        total_tiles = (index_x_upper - index_x_lower)*(index_y_upper - index_y_lower)
        progress_step = 100.0/total_tiles if total_tiles > 0 else 0
        feedback.pushInfo('Import tiles from {} to {} (total {} tiles)'.format(
            (index_x_lower, index_y_lower), (index_x_upper - 1, index_y_upper - 1),
            total_tiles
            ))
        
        feedback.setProgressText('Importing tiles...')
        imported_files = []
        for (ix, x) in enumerate(itertools.product(index_x_range, index_y_range)):
            if feedback.isCanceled():
                raise QgsProcessingException('Canceled')

            result_tif = get_DEM(x[0], x[1], zoomlevel, tilecache_directory)
            if result_tif is not None:
                imported_files.append(str(result_tif))
            
            feedback.setProgress(int((ix + 1)*progress_step))
                
        # OGR may emit a lot of warnings like '/.../***.tif.() does not exist' due to a bug in QGIS,
        # which (should) had been fixed after the Revision 4930061b.
        # [processing] Fix incorrect OGR warnings when loading raster layer results
        # https://issues.qgis.org/projects/qgis/repository/revisions/4930061b21852c6201ddc1c9878122bfd239f1b3
        if len(imported_files) > 0:
            feedback.setProgressText('Merging tiles...')
            output = processing.run('gdal:merge', {
                'INPUT':         imported_files,
                'PICT':          False,
                'SEPARATE':      False,
                'NODATA_INPUT':  np.nan,
                'NODATA_OUTPUT': nodata_value,
                'DATA_TYPE':     5, # float32
                'OUTPUT':        parameters[self.OUTPUT]
            }, context=context, feedback=feedback, is_child_algorithm=True)['OUTPUT']
        else:
            raise QgsProcessingException('No tiles have been imported')
        
        return {self.OUTPUT: output}

