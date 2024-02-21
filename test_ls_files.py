import os
from osgeo import gdal

# 设置文件夹路径
# folder_path = '/mnt/mfs31/SRTML1_unzip'
folder_path = '/mnt/mfs31/China_LS_zip/'

# /mnt/mfs31/SRTML1_unzip/N43E106.hgt

# 获取文件夹下所有文件
files = os.listdir(folder_path)

# 遍历文件
for file in files:
    file_path = os.path.join(folder_path, file)
    filename = file.split('.')[0]
    
    try:
        # 打开栅格数据文件
        dataset = gdal.Open(file_path, gdal.GA_ReadOnly)
        
        # 如果打开成功，则 dataset 不为 None
        if dataset is None:
            print(f"无法打开文件: {file_path}")
            # 如果无法打开，可以选择删除该文件
            os.remove(file_path)
            # zipfile = os.path.join('/mnt/mfs31/SRTML1',filename+".SRTMGL1.hgt.zip")
            # unzip_cmd = 'unzip  %s -d %s'%(zipfile,folder_path)
            # print(unzip_cmd)
            # os.system(unzip_cmd)
        else:
            print(f"成功打开文件: {file_path}")
            dataset = None  # 释放资源
    except Exception as e:
        print(f"出现错误: {e}")
