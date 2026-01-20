import os
import numpy as np
from osgeo import gdal
from qgis.core import QgsProject, QgsVectorLayer, QgsField, QgsFeature
from PyQt5.QtCore import QVariant

# ------------------------
# INPUT FOLDERS
# ------------------------
folders = [
    {
        "path": r"D:\Landsat_Kanha_Moniter_2015_2016\USGS\Landsat_8_9\Mask\NDVI",
        "sensor": "Landsat8/9 (OLI)"
    }
]

# ------------------------
# USE LOADED AOI VECTOR LAYER
# ------------------------
# Replace 'north_west_zone' with the name of your loaded shapefile layer in QGIS
aoi_layers = QgsProject.instance().mapLayersByName("kanha_whole")
if not aoi_layers:
    raise Exception("AOI layer not found in QGIS")
aoi_layer = aoi_layers[0]

# ------------------------
# CREATE MEMORY LAYER FOR RESULTS
# ------------------------
layer = QgsVectorLayer("None", "whole_kanha_table", "memory")
pr = layer.dataProvider()

pr.addAttributes([
    QgsField("date", QVariant.String),
    QgsField("year", QVariant.Int),
    QgsField("month", QVariant.Int),
    QgsField("day", QVariant.Int),
    QgsField("median_ndvi", QVariant.Double),
    QgsField("landsat", QVariant.String)
])
layer.updateFields()

# ------------------------
# PROCESS RASTERS
# ------------------------
for entry in folders:
    folder = entry["path"]
    sensor = entry["sensor"]

    files = sorted([f for f in os.listdir(folder) if f.lower().endswith(".tif")])

    for file in files:
        try:
            # ------------------------
            # EXTRACT DATE
            # ------------------------
            parts = file.split("_")
            date_part = parts[3]  # YYYYMMDD
            year = int(date_part[0:4])
            month = int(date_part[4:6])
            day = int(date_part[6:8])
            date_str = f"{day:02d}-{month:02d}-{year}"

            raster_path = os.path.join(folder, file)
            ds = gdal.Open(raster_path)
            if ds is None:
                continue

            # ------------------------
            # CLIP USING LOADED AOI
            # ------------------------
            clipped = gdal.Warp(
                "",
                ds,
                format="MEM",
                cutlineDSName=aoi_layer.source(),  # use the loaded layer's source path
                cropToCutline=True,
                dstNodata=np.nan
            )

            band = clipped.GetRasterBand(1)
            arr = band.ReadAsArray().astype(float)

            # ------------------------
            # CLEAN NDVI
            # ------------------------
            arr[arr <= 0] = np.nan

            if np.all(np.isnan(arr)):
                median_ndvi = np.nan
            else:
                median_ndvi = float(np.nanmedian(arr))

            print(f"{sensor} | {date_str} → {median_ndvi}")

            # ------------------------
            # ADD FEATURE
            # ------------------------
            feature = QgsFeature()
            feature.setAttributes([
                date_str,
                year,
                month,
                day,
                median_ndvi,
                sensor
            ])
            pr.addFeature(feature)

        except Exception as e:
            print("❌ Error:", file, e)

# ------------------------
# ADD TO QGIS
# ------------------------
QgsProject.instance().addMapLayer(layer)
print("✅ AOI-based NDVI table created successfully")