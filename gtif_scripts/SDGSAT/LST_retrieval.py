import os
import rasterio
import numpy as np
import geopandas as gpd
import xml.etree.ElementTree as ET
import subprocess
from shapely.geometry import box
from rasterio.mask import mask
from osgeo.gdal import Warp, WarpOptions
'''
LST retrieval script developed by Adriaan Keurhorst from Compass Informatics. This script computes LST from SDGSat-1
thermal imagery using ndvi and thermodynamics. The ndvi is computed from openeo_ndvi.py

For any questions, contact adriaan.keurhorst@compass.ie
'''
image_path = 'SDGSAT/Daytime-SDGSat-1/KX10_TIS_20240626_E0.57_N50.57_202400065617_L4B/KX10_TIS_20240626_E0.57_N50.57_202400065617_L4B.tiff'
NDVI_path='sdgsat/ndvi_london_2024_06.tif'
output_folder = 'sdgsat'
clipped_path = image_path[:-5]+'_CLIPPED.tif'   # DO NOT EDIT
metadata_path = image_path[:-5]+'.calib.xml'    # DO NOT EDIT

def clip_raster_with_raster(input, extent, clipped_path):
    def get_coords(input, extent):
        extent = rasterio.open(extent)
        input = rasterio.open(input)
        minx, miny = extent.bounds.left, extent.bounds.bottom
        maxx, maxy = extent.bounds.right, extent.bounds.top
        bbox = box(minx, miny, maxx, maxy)
        geo = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=extent.crs)
        geo = geo.to_crs(crs=input.crs)

        def getFeatures(gdf):
            """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
            import json
            return [json.loads(gdf.to_json())['features'][0]['geometry']]

        coords = getFeatures(geo)
        return coords

    ### step 2 = crop lcm too
    # lst_cropped = rasterio.open('landsat/l8/LST_cropped.tif')
    input_cropped, input_cropped_transform = mask(dataset=rasterio.open(input), shapes=get_coords(input=input, extent=extent),
                                              crop=True, nodata=0)

    # plt.imshow(lcm_cropped[0])
    kwargs = rasterio.open(input).meta

    kwargs.update({"driver": "GTiff",
                   "height": input_cropped.shape[1],
                   "width": input_cropped.shape[2],
                   "transform": input_cropped_transform}
                  )

    with rasterio.open(clipped_path, "w", **kwargs) as dest:
        dest.write(input_cropped)

clip_raster_with_raster(input=image_path, extent=NDVI_path, clipped_path=clipped_path)

opts = WarpOptions(dstSRS=f'{rasterio.open(image_path).crs}')
Warp(NDVI_path[:-4]+'_warp.tif', NDVI_path, options=opts)

merged_path = 'sdgsat/SDGSAT_NDVI_merged.tif'
# cmd = f'python adriaankeurhorst//opt//anaconda3//envs//GTIF//lib/python3.12//site-packages//osgeo_utils//gdal_merge.py -separate -o {merged_path} -ot Float32 {clipped_path} {NDVI_path[:-4]}_warp.tif'
cmd = f'gdal_merge -separate -o {merged_path} -ot Float32 {clipped_path} {NDVI_path[:-4]}_warp.tif'
subprocess.call(cmd)

def radiance_TOA(image_path, metadata_path):
    # This is the TOA radiance

    raw_im = rasterio.open(image_path)
    f = ET.parse(metadata_path, parser=ET.XMLParser(encoding='iso-8859-5'))

    for i in f.iter():
        # print(i.tag)
        # print(i.text)
        if i.tag == 'RADIANCE_GAIN_BAND_1':
            rgb1 = float(i.text)
        if i.tag == 'RADIANCE_GAIN_BAND_2':
            rgb2 = float(i.text)
        if i.tag == 'RADIANCE_GAIN_BAND_3':
            rgb3 = float(i.text)
        if i.tag == 'RADIANCE_BIAS_BAND_1':
            rbb1 = float(i.text)
        if i.tag == 'RADIANCE_BIAS_BAND_2':
            rbb2 = float(i.text)
        if i.tag == 'RADIANCE_BIAS_BAND_3':
            rbb3 = float(i.text)
    # L = DN * GAIN + BIAS
    L_1 = rgb1 * raw_im.read(1) + rbb1
    L_2 = rgb2 * raw_im.read(2) + rbb2
    L_3 = rgb3 * raw_im.read(3) + rbb3

    return L_1, L_2, L_3

L_1, L_2, L_3 = radiance_TOA(image_path=merged_path, metadata_path=metadata_path)
# plt.imshow(L_1)

def satellite_temp(L_1, L_2, L_3):
    # h = 6.626e-34  # Planck's constant (J·s)
    # c = 3.0e8  # Speed of light (m/s)
    # k_B = 1.381e-23  # Boltzmann's constant (J/K)

    wavelength1 = 9.35
    wavelength2 = 10.73
    wavelength3 = 11.72

    # c1 = 2 * h * c ** 2  # First radiation constant
    # c2 = (h * c) / k_B  # Second radiation constant

    # Tb1 = (c2 / wavelength1) / np.log((c1 / (wavelength1**5 * L_1)) + 1)
    # Tb2 = (h * c) / (wavelength2 * k_B) / np.log((2 * h * c ** 2) / (wavelength2 ** 5 * L_2) + 1)
    # Tb3 = (h * c) / (wavelength3 * k_B) / np.log((2 * h * c ** 2) / (wavelength3 ** 5 * L_3) + 1)

    # Tb1 = Tb1 - np.nanmin(Tb1) + 273.15
    # Tb2 = Tb2 - np.nanmin(Tb2) + 273.15
    # Tb3 = Tb3 - np.nanmin(Tb3) + 273.15

    K1 = 1.19104e8  # W μm⁴ m⁻² sr⁻¹
    K2 = 14387.7      # μm K
    Tb1 = K2 / (wavelength1 * np.log((K1 / (L_1 * (wavelength1 ** 5))) + 1))
    Tb2 = K2 / (wavelength2 * np.log((K1 / (L_2 * (wavelength2 ** 5))) + 1))
    Tb3 = K2 / (wavelength3 * np.log((K1 / (L_3 * (wavelength3 ** 5))) + 1))
    return Tb1, Tb2, Tb3

Tb1, Tb2, Tb3 = satellite_temp(L_1, L_2, L_3)
# plt.imshow(Tb1)

def NDVI_proportion(NDVI_path):
    NDVI = rasterio.open(NDVI_path).read(4).clip(0, 1)
    NDVI[NDVI==0] = np.nan
    # plt.imshow(NDVI)
    proportion = ((NDVI - np.nanmin(NDVI)) / (np.nanmax(NDVI) - np.nanmin(NDVI))) ** 0.5

    return proportion

ndvi_proportion = NDVI_proportion(NDVI_path=merged_path)
# plt.imshow(ndvi_proportion)
def emissivity(NDVI_proportion):
    # return 0.004 * NDVI_proportion + 0.986
    return 0.985 * NDVI_proportion + 0.96 * (1 - NDVI_proportion) + 0.005
e = emissivity(ndvi_proportion)
# plt.imshow(e)

def LST(sat_temp, e):
    wavelength_list = [9.35e-6, 10.73e-6, 11.72e-6]
    # row = (6.626 * 10 ** -34) * ((2.998 * 10 ** 8) / (1.38 * 10 ** -23))
    c2 = 1.4388e-2  # Planck’s second constant (m·K)
    lsts = []
    for i in range(len(wavelength_list)):
        # lsts += [(sat_temp[i]) / (1 + (wavelength_list[i] * sat_temp[i] / row) * np.log(e))[0]]
        correction = (wavelength_list[i] * sat_temp[i]) / c2 * np.log(e)
        lsts += [sat_temp[i] / (1 + correction)]
    # plt.imshow(lsts[0])
    final_lst = np.mean(np.array(lsts), axis=0) - 273.15
    return final_lst

final_lst = LST(sat_temp=[Tb1, Tb2, Tb3], e=e)

# plt.imshow(final_lst)

kwargs = rasterio.open(merged_path).meta
kwargs.update(count=1)
with rasterio.open(f'{output_folder}/LST_{image_path.split("/")[-1][:-5]}.tif', 'w',**kwargs) as dst:
    dst.write_band(1, final_lst.astype(rasterio.float32))

os.remove(NDVI_path[:-4]+'_warp.tif')
os.remove(merged_path)



