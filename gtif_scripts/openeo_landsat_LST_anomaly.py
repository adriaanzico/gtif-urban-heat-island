import openeo
from openeo.processes import anomaly, subtract, mean, ProcessBuilder, climatological_normal
import numpy as np

# Define the full temporal range used to compute the baseline climatology
# (last 25 years in this example)
years_all = [2000, 2025]

# Connect to the openEO backend hosted by Copernicus Data Space
# and authenticate using OpenID Connect
connection = openeo.connect(
    "https://openeo.dataspace.copernicus.eu/openeo/1.2"
).authenticate_oidc()

# Define the Area of Interest (AOI) as a bounding box
# Several examples are provided; Paris is used here
# bounding_box_aoi = [-0.563, 51.261, 0.280, 51.732]  # London
bounding_box_aoi = [7.546, 44.989, 7.781, 45.145]   # Torino
# bounding_box_aoi = [1.349, 43.544, 1.56, 43.67]    # Toulouse
# bounding_box_aoi = [2.273312, 48.818733, 2.412701, 48.897340]  # Paris
# bounding_box_aoi = [4.2276, 52.0270, 4.3953, 52.1359] # Den Haag

# STAC endpoint for the Landsat Collection 2 Level-2 dataset
landsat = "https://planetarycomputer.microsoft.com/api/stac/v1/collections/landsat-c2-l2"

# Load the Landsat data as an openEO data cube
# - Spatial extent: AOI bounding box
# - Temporal extent: full baseline period (2000–2025)
# - Band: TIRS_B10 (thermal infrared band used for LST)
# - Filter scenes to those with ≤5% cloud cover
landsat_cube = connection.load_stac(
    landsat,
    spatial_extent={
        "west": bounding_box_aoi[0],
        "south": bounding_box_aoi[1],
        "east": bounding_box_aoi[2],
        "north": bounding_box_aoi[3],
    },
    temporal_extent=[f"{years_all[0]}-01-01", f"{years_all[1]}-12-31"],
    bands=["TIRS_B10"],
    properties={
        "eo:cloud_cover": lambda v: v <= 5,
    }
)

# Convert the thermal band digital numbers to land surface temperature (°C)
# Scaling and offset follow Landsat Collection 2 documentation
# 0.00341802 * DN + 149.0 converts to Kelvin
# Subtract 273.15 to convert from Kelvin to Celsius
lwir = landsat_cube.band("TIRS_B10") * 0.00341802 + 149.0 - 273.15

# Function to compute the mean land surface temperature for
# the summer season (June–July–August, JJA)
def jja_mean(datacube, years):
    # Aggregate the data into yearly JJA intervals
    datacube = datacube.aggregate_temporal(
        intervals=[
            [f"{year}-06-01", f"{year}-08-31"]
            for year in np.arange(years[0], years[1] + 1)
        ],
        reducer="mean"
    )
    # Compute the mean across all selected years
    return datacube.reduce_temporal(reducer="mean")

# Compute the baseline mean JJA LST over the full 25-year period
baseline_means = jja_mean(lwir, years_all)

# Download the baseline mean LST as a GeoTIFF
baseline_means.download(
    f"baseline_means_{years_all[0]}_{years_all[1]}.tif"
)

# Define the years of interest for which anomalies will be calculated
years_of_interest = [2020, 2025]

# Compute the mean JJA LST for the years of interest
lst_years_of_interest = jja_mean(lwir, years_of_interest)

# Download the JJA mean LST for the years of interest
lst_years_of_interest.download(
    f"torino_LST{years_of_interest[0]}_{years_of_interest[1]}.tif"
)

# Compute the land surface temperature anomaly
# Anomaly = (JJA mean LST for years of interest) − (25-year baseline JJA mean)
lst_anomaly = lst_years_of_interest - baseline_means

# Download the LST anomaly map as a GeoTIFF
lst_anomaly.download(
    f"paris_lst_anomaly_{years_of_interest[0]}_{years_of_interest[1]}.tif"
)
