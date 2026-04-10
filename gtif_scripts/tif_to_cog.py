from rio_cogeo.cogeo import cog_translate
from rio_cogeo.profiles import cog_profiles
import rasterio


def convert_to_cog(input_tif, output_tif, profile="deflate"):
    # Define COG profile
    cog_profile = cog_profiles.get(profile)

    # Open input GeoTIFF
    with rasterio.open(input_tif) as src:
        # Translate to COG
        cog_translate(
            src,
            output_tif,
            cog_profile,
            in_memory=False,
            quiet=False,
        )

    print(f"✅ Converted {input_tif} → {output_tif} (COG format)")

convert_to_cog("/Users/adriaankeurhorst/Documents/GTIF/scratch/lst_anomaly_2020_2025.tif", "/Users/adriaankeurhorst/Documents/GTIF/scratch/lst_anomaly_2020_2025_cog.tif")
