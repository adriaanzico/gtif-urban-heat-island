import openeo
from openeo.processes import anomaly, subtract, mean, ProcessBuilder, climatological_normal
from openeo.processes import ProcessBuilder as P
import numpy as np

years_all = [2022, 2025]
connection = openeo.connect("https://openeo.dataspace.copernicus.eu/openeo/1.2").authenticate_oidc()
bounding_box_aoi = [-0.563, 51.261, 0.280, 51.732] #london
# bounding_box_aoi = [7.546, 44.989, 7.781, 45.145] # torino
# bounding_box_aoi = [1.349, 43.544, 1.56, 43.67] # toulouse
# bounding_box_aoi = [2.273312, 48.818733, 2.412701, 48.897340] #paris
landsat = (
        "https://planetarycomputer.microsoft.com/api/stac/v1/collections/landsat-c2-l2"
    )

landsat_cube = connection.load_stac(
    landsat,
    spatial_extent={"west": bounding_box_aoi[0], "south": bounding_box_aoi[1], "east": bounding_box_aoi[2],
                    "north": bounding_box_aoi[3]},
    temporal_extent=[f"{years_all[0]}-01-01", f"{years_all[1]}-12-31"],
    bands=["TIRS_B10"],
    properties={
            "eo:cloud_cover": lambda v: v <= 5,
        }
    )

lwir = landsat_cube.band('TIRS_B10') * 0.00341802 + 149.0 - 273.15

def jja_mean(datacube, years):
    datacube = datacube.aggregate_temporal(
        intervals=[[f"{year}-06-01", f"{year}-08-31"] for year in np.arange(years[0], years[1] + 1)],
        reducer="mean")
    return datacube.reduce_temporal(reducer='mean')

years_of_interest=[2022, 2025]
lst_years_of_interest = jja_mean(lwir, years_of_interest)
# lst_years_of_interest.download("lwir_mean.tif", format="GTiff")

# year=2018
# def load_corine(year, bounding_box_aoi):
#     connection = openeo.connect("openeo.dataspace.copernicus.eu").authenticate_oidc()
#     connection = openeo.connect('https://openeocloud.vito.be/openeo/1.0.0').authenticate_oidc()
#     corine_cube = connection.load_collection("CORINE_LANDCOVER",
#                     spatial_extent={"west": bounding_box_aoi[0], "south": bounding_box_aoi[1], "east": bounding_box_aoi[2],"north": bounding_box_aoi[3]}
#                                              ).reduce_temporal(reducer='first')
#
#     corine_cube = connection.load_stac(
#         "https://stac.openeo.vito.be/collections/CORINE_LANDCOVER",
#         spatial_extent={"west": bounding_box_aoi[0], "south": bounding_box_aoi[1], "east": bounding_box_aoi[2],
#                     "north": bounding_box_aoi[3]}).reduce_temporal(reducer='first')
#
#     corine_cube = connection.load_stac(
#         url='https://s3.ecodatacube.eu/arco/landcover_clc.plus_f_30m_0..0cm_20170101_20191231_eu_epsg.3035_v20250327.tif',
#         spatial_extent={"west": bounding_box_aoi[0], "south": bounding_box_aoi[1], "east": bounding_box_aoi[2],"north": bounding_box_aoi[3]}
#                 ).reduce_temporal(reducer='first')
#
#     # corine_cube.download("corine.tif", format="GTiff")
#     return corine_cube

def load_cci(bounding_box_aoi, year):
    url = 'https://planetarycomputer.microsoft.com/api/stac/v1/collections/esa-cci-lc'
    cci_cube = connection.load_stac(
        url,
        spatial_extent={"west": bounding_box_aoi[0], "south": bounding_box_aoi[1], "east": bounding_box_aoi[2],
                        "north": bounding_box_aoi[3]},
        temporal_extent=[f"{year}-01-01", f"{year}-12-31"],
        bands=["lccs_class"],
    )

    # cci_cube.download('cci_test.tif')
    return cci_cube

# corine = load_corine(2018, bounding_box_aoi)
# corine.download('corine_test.tif')

cci = load_cci(bounding_box_aoi, 2018)

def run1():
    # corine_urban_mask = ~ ((corine <1) | (corine >0))
    urban_mask = ((cci == 190))
    rural_mask = ((cci != 190))
    # corine_urban_mask.download('corine_urban_mask.tif')

    # Mask LST cube
    urban_temp = lst_years_of_interest.mask(urban_mask)
    rural_temp = lst_years_of_interest.mask(rural_mask)


    # -------------------------------
    # 4. Iterative buffer T_rural
    # -------------------------------
    # We precompute rural mean LST for several window sizes and pick
    # the first non-NaN per pixel.

    def rural_mean_for_window(lst: openeo.DataCube, rural_mask: openeo.DataCube, size: int) -> openeo.DataCube:
        # Mask LST by rural area only
        masked = lst.mask(rural_mask)

        # Apply neighborhood mean over masked LST
        def reducer(data: P) -> P:
            # data is a vector over the neighborhood pixels
            # NaNs (from non-rural) are ignored in mean
            return data.mean()

        return masked.apply_neighborhood(
            process=reducer,
            size=[
                {"dimension": "x", "value": size},
                {"dimension": "y", "value": size},
            ],
            overlap=[
                {"dimension": "x", "value": size // 2},
                {"dimension": "y", "value": size // 2},
            ],
        )


    buffer_sizes = [5, 10, 15, 20, 25]

    rural_means = [rural_mean_for_window(lst_years_of_interest, rural_mask, s)
                   for s in buffer_sizes]

    # Stack along a new dimension ("window") and pick first valid (non-NaN) value
    stack = rural_means[0]
    for rm in rural_means[1:]:
        stack = stack.merge_cubes(rm)  # backend-specific: may create extra band dimension

    # If your backend supports it, make sure `stack` has a "window" or "bands" dimension.
    # Here we assume "bands".

    from openeo.processes import array_element, is_nan, ProcessBuilder as P

    def first_valid(v: P) -> P:
        # v is the 1D array over the reduced dimension
        v0 = array_element(data=v, index=0)
        v1 = array_element(data=v, index=1)
        v2 = array_element(data=v, index=2)
        v3 = array_element(data=v, index=3)
        v4 = array_element(data=v, index=4)

        r01   = P.if_(is_nan(v0), v1,   v0)
        r012  = P.if_(is_nan(r01), v2,  r01)
        r0123 = P.if_(is_nan(r012), v3, r012)
        r0124 = P.if_(is_nan(r0123), v4, r0123)

        return r0124

    # t_rural = stack.reduce_dimension(
    #     dimension="bands",
    #     reducer=first_valid,
    # )
    t_rural = stack

    t_urban = urban_temp  # LST for urban pixels (already masked)

    uhi = t_urban / (t_urban - t_rural)

    from openeo.processes import is_nan, if_ as ifproc, ProcessBuilder as P

    # 1. Build a validity mask based on t_rural
    #    -> 1 where t_rural is NOT NaN, 0 where it is NaN
    valid_rural_mask = t_rural.apply(lambda x: ifproc(is_nan(x), 0, 1))

    # 2. Combine with your CORINE urban mask: urban & valid_rural
    #    corine_urban_mask is a boolean cube already
    from openeo.processes import is_nan, if_ as ifproc

    # valid_rural_mask: 1 where t_rural is NOT NaN, 0 where it is NaN
    valid_rural_mask = t_rural.apply(lambda x: ifproc(is_nan(x), 0, 1))

    # If corine_urban_mask is boolean, convert to 0/1:
    urban_01 = urban_mask.apply(lambda x: ifproc(x, 1, 0))

    # Combined mask: urban AND valid_rural  -> multiply 0/1 masks
    combined_mask = urban_01 * valid_rural_mask

    # Apply mask to UHI
    uhi_masked = uhi.mask(combined_mask)

    uhi_masked.download("uhi_summer_torino.tif", format="GTiff")


# -------------------------------------------------
# 1. DEFINE MASKS (ESA CCI: 190 = Urban)
# -------------------------------------------------
cci_resampled = cci.reduce_temporal(reducer='first').resample_cube_spatial(lst_years_of_interest,)
rural_mask = (cci_resampled == 190)
urban_mask = (cci_resampled != 190)

# Convert boolean mask to numeric (1 rural, 0 urban)
rural_indicator = rural_mask * 1

# -------------------------------------------------
# 2. KEEP ONLY RURAL TEMPERATURE
# -------------------------------------------------
# Multiply LST by rural mask (urban pixels become 0)
rural_temp_zero = lst_years_of_interest * rural_indicator

# -------------------------------------------------
# 3. DEFINE MOVING WINDOW KERNEL
# -------------------------------------------------
kernel_size = 51
kernel = [[1] * kernel_size for _ in range(kernel_size)]

# -------------------------------------------------
# 4. COMPUTE LOCAL RURAL SUM + COUNT
# -------------------------------------------------
rural_sum = rural_temp_zero.apply_kernel(
    kernel=kernel,
    factor=1,
    border=0
)

rural_count = rural_indicator.apply_kernel(
    kernel=kernel,
    factor=1,
    border=0
)

# -------------------------------------------------
# 5. LOCAL RURAL MEAN (safe division)
# -------------------------------------------------
local_rural_mean = rural_sum / rural_count
local_rural_mean = local_rural_mean.mask(rural_count == 0)

# -------------------------------------------------
# 6. LOCAL PER-PIXEL UHI
# -------------------------------------------------
local_uhi = lst_years_of_interest - local_rural_mean

# Keep only urban pixels in final map
local_uhi_urban = local_uhi.mask(urban_mask)

# -------------------------------------------------
# 7. EXPORT
# -------------------------------------------------
job = local_uhi_urban.execute_batch(outputfile='uhi_summer_london.tif',
    out_format="GTiff",
    title="UHI_JJA_2022_2025"
)

job.get_results().download_files(target='uhi_summer_london')

# # Configurable kernel sizes (in pixels)
# kernel_sizes = [5, 15, 30, 50]
#
# # Compute local rural references for each kernel
# rural_locals = []
# for k in kernel_sizes:
#     rural_locals.append(
#         rural_temp.apply_neighborhood(
#             process="mean",
#             size=[
#                 {"dimension": "x", "unit": "px", "value": k},
#                 {"dimension": "y", "unit": "px", "value": k}
#             ],
#             overlap=[
#                 {"dimension": "x", "unit": "px", "value": k},
#                 {"dimension": "y", "unit": "px", "value": k}
#             ]
#         )
#     )
#
# from openeo.processes import if_, is_valid
# # Fallback logic: use smallest kernel first, then larger ones if needed
# rural_local = rural_locals[-1]  # start with the largest as default
# for r in reversed(rural_locals[:-1]):
#     rural_local = if_(is_valid(r), r, rural_local)
#
# # Compute UHI: only for urban pixels
# uhi = (urban_temp) / (urban_temp - rural_local)
# uhi = uhi.mask(corine_urban_mask)
#
# uhi.download("uhi_local.tif", format="GTiff")