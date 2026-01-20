from qgis.core import (
    QgsProject,
    QgsRasterShader,
    QgsColorRampShader,
    QgsSingleBandPseudoColorRenderer
)
from PyQt5.QtGui import QColor

# -----------------------------
# Path to your pseudo-color file
# -----------------------------
txt_path = r"D:/Landsat_Kanha_Moniter_2015_2016/pseudo_color.txt"

# -----------------------------
# Parse TXT file
# -----------------------------
color_items = []

with open(txt_path, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("INTERPOLATION"):
            continue  # skip header/empty lines

        parts = line.split(',')
        if len(parts) < 5:
            continue

        # Handle "inf" for upper bound
        val_str = parts[0].strip()
        value = float('inf') if val_str.lower() == 'inf' else float(val_str)

        # RGBA
        r, g, b, a = map(int, parts[1:5])
        color = QColor(r, g, b, a)

        # Label
        label = parts[5].strip() if len(parts) > 5 else str(value)

        # Create ColorRampItem
        color_items.append(QgsColorRampShader.ColorRampItem(value, color, label))

# -----------------------------
# Apply pseudo-color to all raster layers
# -----------------------------
for layer in QgsProject.instance().mapLayers().values():
    if layer.type() == layer.RasterLayer:

        raster_shader = QgsRasterShader()
        color_ramp = QgsColorRampShader()
        color_ramp.setColorRampItemList(color_items)
        color_ramp.setColorRampType(QgsColorRampShader.Discrete)  # discrete like in TXT

        raster_shader.setRasterShaderFunction(color_ramp)

        renderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, raster_shader)
        layer.setRenderer(renderer)
        layer.triggerRepaint()

print("✅ Exact pseudo-color file applied to all raster layers")