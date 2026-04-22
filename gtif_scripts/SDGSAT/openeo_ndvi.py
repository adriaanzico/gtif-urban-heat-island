import openeo

'''
NDVI retrieval script developed by Adriaan Keurhorst from Compass Informatics. This script computes the mean NDVI over
an AOI and a specified month. The output is required for LST_retrieval.py to compute LST from SDGSat-1 thermal imagery.

For any questions, contact adriaan.keurhorst@compass.ie
'''

temp_ex = [f"2024-06-01", f"2024-06-30"]
bounding_box_aoi = [-0.563, 51.261, 0.280, 51.732]
output_folder = ''

connection = openeo.connect("https://openeo.dataspace.copernicus.eu/openeo/1.2").authenticate_oidc()
connection.list_collection_ids()
bbox = {"west": bounding_box_aoi[0], "south": bounding_box_aoi[1], "east": bounding_box_aoi[2],
                        "north": bounding_box_aoi[3]}

cube_s2 = connection.load_collection(
    "SENTINEL2_L2A",
    spatial_extent = bbox,
    temporal_extent = temp_ex,
    bands = ["B04", "B08"],
    max_cloud_cover=10
)

s2_b4 = cube_s2.band('B04')
s2_b8 = cube_s2.band('B08')

#compute ndvi
s2_ndvi = (s2_b8 - s2_b4)/(s2_b8 + s2_b4)

s2_ndvi = s2_ndvi.reduce_temporal(reducer='mean')
# s2_ndvi.download('test.tif')
s2_ndvi.download(f'{output_folder}/ndvi_london_{temp_ex[0].split("-")[0]}_{temp_ex[0].split("-")[1]}.tif')
