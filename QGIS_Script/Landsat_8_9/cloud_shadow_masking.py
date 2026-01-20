import os
import processing

# -------------------------------
# PATHS
# -------------------------------
qa_dir   = r"D:/Landsat_Kanha_Moniter_2015_2016/USGS/Landsat_8_9/QA_Pixel"
red_dir  = r"D:/Landsat_Kanha_Moniter_2015_2016/USGS/Landsat_8_9/RED"
nir_dir  = r"D:/Landsat_Kanha_Moniter_2015_2016/USGS/Landsat_8_9/NIR"

out_red = r"D:/Landsat_Kanha_Moniter_2015_2016/USGS/Landsat_8_9/Mask/RED"
out_nir = r"D:/Landsat_Kanha_Moniter_2015_2016/USGS/Landsat_8_9/Mask/NIR"

os.makedirs(out_red, exist_ok=True)
os.makedirs(out_nir, exist_ok=True)

# -------------------------------
# PROCESS FILES
# -------------------------------
for qa_file in os.listdir(qa_dir):
    if not qa_file.lower().endswith(".tif"):
        continue

    base = qa_file.split("_QA_PIXEL")[0]

    qa_path  = os.path.join(qa_dir, qa_file)
    red_path = os.path.join(red_dir, f"{base}_SR_B4.tif")
    nir_path = os.path.join(nir_dir, f"{base}_SR_B5.tif")

    if not (os.path.exists(red_path) and os.path.exists(nir_path)):
        print(f"❌ Skipping {base} (missing RED/NIR)")
        continue

    # -------------------------------
    # CREATE CLOUD + SHADOW MASK
    # -------------------------------
    mask = processing.run(
        "gdal:rastercalculator",
        {
            'INPUT_A': qa_path,
            'BAND_A': 1,
            'FORMULA': '((A & 8) == 0) * ((A & 16) == 0)',
            'NO_DATA': 0,
            'RTYPE': 5,  # Byte mask is OK
            'OUTPUT': 'TEMPORARY_OUTPUT'
        }
    )['OUTPUT']

    # -------------------------------
    # APPLY MASK TO RED
    # -------------------------------
    processing.run(
        "gdal:rastercalculator",
        {
            'INPUT_A': red_path,
            'BAND_A': 1,
            'INPUT_B': mask,
            'BAND_B': 1,
            'FORMULA': 'A * B',
            'NO_DATA': -9999,
            'RTYPE': 6,  # Float32
            'OUTPUT': os.path.join(out_red, f"{base}_RED_MASKED.tif")
        }
    )

    # -------------------------------
    # APPLY MASK TO NIR
    # -------------------------------
    processing.run(
        "gdal:rastercalculator",
        {
            'INPUT_A': nir_path,
            'BAND_A': 1,
            'INPUT_B': mask,
            'BAND_B': 1,
            'FORMULA': 'A * B',
            'NO_DATA': -9999,
            'RTYPE': 6,  # Float32
            'OUTPUT': os.path.join(out_nir, f"{base}_NIR_MASKED.tif")
        }
    )

    print(f"✅ Processed: {base}")

print("🎉 Cloud masking completed successfully")