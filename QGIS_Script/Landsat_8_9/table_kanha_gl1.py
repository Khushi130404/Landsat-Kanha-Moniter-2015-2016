import os
import numpy as np
from osgeo import gdal
from qgis.core import (
    QgsProject,
    QgsVectorLayer,
    QgsField,
    QgsFeature
)
from PyQt5.QtCore import QVariant

# ============================================================
# INPUT FOLDERS
# ============================================================
folders = [
    {
        "path": r"D:\MODIS_Kanha_Moniter_2010_2011\GEE_GeoTiff\MODIS_NDVI_TIFS",
        "sensor": "MODIS"
    }
]

# ============================================================
# USE LOADED AOI VECTOR LAYER
# ============================================================
aoi_layers = QgsProject.instance().mapLayersByName("kanha_grass_1")
if not aoi_layers:
    raise Exception("❌ AOI layer 'kanha_grass_1' not found in QGIS")

aoi_layer = aoi_layers[0]

# ============================================================
# CREATE MEMORY LAYER FOR OUTPUT TABLE
# ============================================================
layer = QgsVectorLayer("None", "gl1_kanha_table", "memory")
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

# ============================================================
# PROCESS MODIS NDVI RASTERS
# ============================================================
for entry in folders:
    folder = entry["path"]
    sensor = entry["sensor"]

    files = sorted([f for f in os.listdir(folder) if f.lower().endswith(".tif")])

    for file in files:
        try:
            # ------------------------------------------------
            # EXTRACT DATE FROM MODIS FILENAME
            # MODIS_NDVI_YYYY_MM_DD.tif
            # ------------------------------------------------
            name = os.path.splitext(file)[0]
            parts = name.split("_")

            year = int(parts[2])
            month = int(parts[3])
            day = int(parts[4])

            date_str = f"{day:02d}-{month:02d}-{year}"

            raster_path = os.path.join(folder, file)
            ds = gdal.Open(raster_path)
            if ds is None:
                continue

            # ------------------------------------------------
            # CLIP USING AOI
            # ------------------------------------------------
            clipped = gdal.Warp(
                "",
                ds,
                format="MEM",
                cutlineDSName=aoi_layer.source(),
                cropToCutline=True,
                dstNodata=np.nan
            )

            band = clipped.GetRasterBand(1)
            arr = band.ReadAsArray().astype(float)

            # ------------------------------------------------
            # CLEAN NDVI VALUES (MODIS-APPROPRIATE)
            # ------------------------------------------------
            arr[(arr < -0.3) | (arr > 1.0)] = np.nan

            if np.all(np.isnan(arr)):
                median_ndvi = np.nan
            else:
                median_ndvi = float(np.nanmedian(arr))

            print(f"{sensor} | {date_str} → Median NDVI = {median_ndvi}")

            # ------------------------------------------------
            # ADD FEATURE TO MEMORY TABLE
            # ------------------------------------------------
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
            print(f"❌ Error processing {file}: {e}")

# ============================================================
# ADD RESULT LAYER TO QGIS
# ============================================================
QgsProject.instance().addMapLayer(layer)
print("✅ AOI-based MODIS Median NDVI table created successfully")
