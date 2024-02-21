"""
Author: zjh
Email: zjh@cnic.cn
"""

import os
from osgeo import ogr,gdal

dem_path = '/mnt/mfs31/SRTML1/'
unzip_path = "/mnt/mfs31/SRTML1_unzip/"
LS_path = "/mnt/mfs31/China_LS/"
tmp_dir='/mnt/mfs31/tmp'
# zip_name='N19E101.SRTMGL1.hgt.zip'
# filename='N19E108.hgt'
def get_extent(filename):
    # 打开栅格数据文件
    ds = gdal.Open(filename)
    # 获取地理范围
    xmin, xres, _, ymax, _, yres = ds.GetGeoTransform()
    xmax = xmin + (ds.RasterXSize * xres)
    ymin = ymax - (ds.RasterYSize * abs(yres))
    # 关闭栅格数据文件
    ds = None
    return xmin,ymin,xmax,ymax
    
def whether_xiangjiao(n,e, border_shp ):
    # 创建格网多边形
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(e, n)
    ring.AddPoint(e+1, n)
    ring.AddPoint(e+1, n+1)
    ring.AddPoint(e, n+1)
    ring.AddPoint(e, n)  # 闭合多边形

    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    input_ds = ogr.Open(border_shp)
    if input_ds is None:
        print("无法打开输入shapefile文件")
        return

    # 获取输入shapefile的第一个图层
    input_layer = input_ds.GetLayer()

    # 获取输入shapefile的空间参考信息
    source_srs = input_layer.GetSpatialRef()

    # 判断格网是否与多边形相交或在多边形内部
    input_layer.SetSpatialFilter(poly)

    for feature in input_layer:
        if feature.geometry().Intersects(poly):
            return True
        else:
            return False
            
def ls_factor(N_min, N_max, E_min, E_max):
    for n in range(N_min, N_max):
        for e in range(E_min, E_max):
            print(n,e)
            # # 创建格网多边形
            # ring = ogr.Geometry(ogr.wkbLinearRing)
            # ring.AddPoint(e, n)
            # ring.AddPoint(e+1, n)
            # ring.AddPoint(e+1, n+1)
            # ring.AddPoint(e, n+1)
            # ring.AddPoint(e, n)  # 闭合多边形

            # poly = ogr.Geometry(ogr.wkbPolygon)
            # poly.AddGeometry(ring)

            e = str(e).zfill(3)
            
            lsfile = os.path.join(LS_path, 'LS_N%s_E%s.tif'%(n,e))
            zipfile=os.path.join(dem_path,'N%sE%s'%(n,e)+'.SRTMGL1.hgt.zip')
            if os.path.exists(lsfile):
                print('%s exists'%lsfile)
                continue
            if not os.path.exists(zipfile):
                print('%s not exists'%zipfile)
                continue
            
            c_dem = os.path.join(unzip_path,'N%sE%s.hgt'%( n,e))     
            
            #第一步合并d8的dem文件并投影到3857平面坐标系下
            srcfiles1 = ''
            for i in range(n-1,n+2):
                for j in range(int(e)-1, int(e)+2):
                    jj=str(j).zfill(3)
                    zipfile=os.path.join(dem_path,'N%sE%s'%( i,jj)+'.SRTMGL1.hgt.zip')
                    datafile = os.path.join(unzip_path,'N%sE%s.hgt'%( i,jj))                               
                    if not os.path.exists(datafile):
                        unzip_cmd = 'unzip  %s -d %s'%(zipfile,unzip_path)
                        os.system(unzip_cmd)
                    if os.path.exists(datafile):
                        srcfiles1 = srcfiles1+datafile
                        srcfiles1 = srcfiles1+' '
                
            num_neighbors = len(srcfiles1.split(' '))
            print('the %s dem file has %s neighbor files.'% (c_dem, num_neighbors))
                    
            d8merge_dem=os.path.join(tmp_dir, 'DEM_N%s_E%s.tif'%(n,e))                   

            merge_cmd1 = 'gdalwarp -wo NUM_THREADS=ALL_CPUS -t_srs EPSG:3857 %s %s'%(srcfiles1, d8merge_dem)
            os.system(merge_cmd1)   
            if not os.path.exists(d8merge_dem):
                print('the dem file %s has no enough d8 files.'% c_dem)
                continue  

            #第二步计算d8范围的LS因子
            tmp_ls=os.path.join(tmp_dir, 'tmp_LS_N%s_E%s.tif'%(n,e))
            cmd = "saga_cmd terrain_analysis 'ta_ls_factor' -DEM %s -LS_FACTOR %s -LS_METHOD 1 -PREPROCESSING 1"%(d8merge_dem, tmp_ls)
            print(cmd)
            os.system(cmd)       

            #第三步
            xmin,ymin,xmax,ymax = get_extent(c_dem)  
            
            #4 重投影8邻域的LS因子到WGS84地理坐标系            
            tmp_ls_4326=os.path.join(tmp_dir,'ls_4326_N%s_E%s.tif'%(n,e))  
            warp_cmd = 'gdalwarp -t_srs EPSG:4326 %s %s'%(tmp_ls, tmp_ls_4326)
            os.system(warp_cmd)   
            
            # 5 裁剪LS因子，只保留中心切片区域
            cutline_cmd = "gdalwarp -te %s %s %s %s -crop_to_cutline %s %s"%( xmin,ymin,xmax,ymax,tmp_ls_4326, lsfile)
            os.system(cutline_cmd)  

            #第6步删除中间临时文件
            rm_cmd='rm -rf %s %s %s'%(d8merge_dem,tmp_ls,tmp_ls_4326)
            print(rm_cmd)
            os.system(rm_cmd)

def bounds_for_tif(input_dir,output_file):
    # 创建一个新的矢量文件
    driver = ogr.GetDriverByName('GeoJSON')  # 或者 'ESRI Shapefile'
    if os.path.exists(output_file):
        driver.DeleteDataSource(output_file)
    out_ds = driver.CreateDataSource(output_file)

    # 创建图层
    output_layer = out_ds.CreateLayer('boundaries', geom_type=ogr.wkbPolygon)

    # 遍历GeoTIFF文件并提取边界
    tif_files = [os.path.join(input_dir, filename) for filename in os.listdir(input_dir) if filename.endswith('.tif')]
    for tif_file in tif_files:
        ds = gdal.Open(tif_file)

        if ds is None:
            continue

        # 获取GeoTIFF的边界
        geo_transform = ds.GetGeoTransform()
        ulx, xres, xskew, uly, yskew, yres = geo_transform
        lrx = ulx + (ds.RasterXSize * xres)
        lry = uly + (ds.RasterYSize * yres)

        # 创建一个多边形几何对象
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(ulx, uly)
        ring.AddPoint(lrx, uly)
        ring.AddPoint(lrx, lry)
        ring.AddPoint(ulx, lry)
        ring.AddPoint(ulx, uly)  # 闭合多边形

        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        # 创建一个新要素并添加到输出图层
        feature_defn = output_layer.GetLayerDefn()
        feature = ogr.Feature(feature_defn)
        feature.SetGeometry(poly)
        output_layer.CreateFeature(feature)

        ds = None  # 关闭数据集

    # 保存输出数据源
    out_ds = None

    print(f'Merged boundaries are saved to {output_file}')           

def zip_file(folder_path, zip_path):
    # 遍历文件
    files = os.listdir(folder_path)
    i=0
    for file in files:
        i=i+1
        # print(i)
        file_path = os.path.join(folder_path, file)
        zipfile = os.path.join(zip_path, file)
        if os.path.exists(zipfile):
            continue 

        try:
            # 打开栅格数据文件
            dataset = gdal.Open(file_path, gdal.GA_ReadOnly)
            
            # 如果打开成功，则 dataset 不为 None
            if dataset is None:
                print(f"无法打开文件: {file_path}")
                # 如果无法打开，可以选择删除该文件
                os.remove(file_path)
            else:
                print(f"成功打开文件: {file_path}")
                
                zip_cmd = 'gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9 -co PREDICTOR=2 -co TILED=YES %s %s'%(file_path,zipfile)
                os.system(zip_cmd)
            dataset = None  # 释放资源
        except Exception as e:
            print(f"出现错误: {e}")
    
    
if __name__ == '__main__':
    #指定空间范围，计算LS因子
    N_min=15;N_max=56;E_min=72;E_max=137
    ls_factor(N_min,N_max, E_min, E_max)
    
    #压缩LS结果，减少存储空间
    LS_zip_path='/mnt/mfs31/China_LS_zip'
    zip_file(LS_path, LS_zip_path)
    
    #输出LS因子数据矢量边界:GeoJSON或Shapefile文件
    output_file = '/mnt/mfs31/tmp/output.geojson'
    bounds_for_tif(LS_zip_path,output_file)
    




        

       
