使用开源工具（SAGA和GDAL）和公开的高分辨率数字高程模型（30米SRTM）计算获得中国区域30米分辨率LS因子数据集。
在进行大范围DEM数据的LS因子计算时，涉及每个单元与其相邻单元的流量分布计算，因此存在数据依赖性。为了高效处理覆盖中国区域的大规模LS因子计算，我们采用了空间分解方法对DEM数据进行分块处理。具体操作如下：每个格网的大小为1°× 1°，在计算格网区域（i, j）的DEM数据时，同时读取周围8个邻域格网空间范围内的DEM数据。将这些数据进行投影转换到平面坐标系后，计算整个区域的LS数据。最终，我们只保留了当前格网区域（i, j）的LS值，并采用地理坐标系存储结果。对于中国海岸带，即东南沿海区域，由于邻域切片小于8个，我们在处理时进行了相应的适应性调整。