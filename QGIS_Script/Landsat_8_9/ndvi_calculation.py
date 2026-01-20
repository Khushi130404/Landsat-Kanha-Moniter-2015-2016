import os
from osgeo import gdal
import numpy as np

# -------------------------------
# FOLDERS
# -------------------------------
red_dir  = r"D:/Landsat_Kanha_Moniter_2015_2016/USGS/Landsat_8_9/Mask/RED"
nir_dir  = r"D:/Landsat_Kanha_Moniter_2015_2016/USGS/Landsat_8_9/Mask/NIR"
ndvi_dir = r"D:/Landsat_Kanha_Moniter_2015_2016/USGS/Landsat_8_9/Mask/NDVI"

os.makedirs(ndvi_dir, exist_ok=True)

# -------------------------------
# HELPER FUNCTIONS
# -------------------------------
def read_raster(path):
    ds = gdal.Open(path)
    arr = ds.ReadAsArray().astype(np.float32)
    return ds, arr

def save_raster(path, array, ref_ds):
    driver = gdal.GetDriverByName("GTiff")
    out = driver.Create(
        path,
        ref_ds.RasterXSize,
        ref_ds.RasterYSize,
        1,
        gdal.GDT_Float32
    )
    out.SetGeoTransform(ref_ds.GetGeoTransform())
    out.SetProjection(ref_ds.GetProjection())
    band = out.GetRasterBand(1)
    band.WriteArray(array)
    band.SetNoDataValue(-9999)
    out.FlushCache()

# -------------------------------
# PROCESS EACH SCENE
# -------------------------------
for red_file in os.listdir(red_dir):
    if not red_file.endswith(".tif"):
        continue

    # Match NIR
    base = red_file.replace("_RED_MASKED.tif", "")
    red_path = os.path.join(red_dir, red_file)
    nir_path = os.path.join(nir_dir, base + "_NIR_MASKED.tif")

    if not os.path.exists(nir_path):
        print(f"❌ Missing NIR for {base}")
        continue

    print(f"Processing NDVI for {base}")

    # READ RASTERS
    red_ds, red = read_raster(red_path)
    nir_ds, nir = read_raster(nir_path)

    # CALCULATE NDVI
    ndvi = (nir - red) / (nir + red)
    ndvi[np.isnan(ndvi)] = -9999  # set NoData for invalid pixels

    # SAVE NDVI
    save_raster(os.path.join(ndvi_dir, base + "_NDVI.tif"), ndvi, red_ds)

print("✅ All NDVI rasters saved in:", ndvi_dir)