import matplotlib.pyplot as plt
from qgis.core import QgsProject
from datetime import date

def plot_ndvi_chronological(
    layer_name,
    title,
    ndvi_field="median_ndvi",
    year_field="year",
    month_field="month",
    day_field="day"
):
    """
    Plot all NDVI points starting from the earliest date in the table
    """

    layer = QgsProject.instance().mapLayersByName(layer_name)[0]

    data = []

    for f in layer.getFeatures():
        try:
            y = int(f[year_field])
            m = int(f[month_field])
            d = int(f[day_field])
            ndvi = float(f[ndvi_field])
        except:
            continue

        if ndvi <= 0 or ndvi > 1:
            continue

        dt = date(y, m, d)
        data.append((dt, ndvi))

    if not data:
        print("⚠ No valid data points found")
        return

    # Sort by actual date
    data.sort(key=lambda x: x[0])

    dates = [d[0] for d in data]
    ndvi_vals = [d[1] for d in data]

    # ------------------------
    # PLOT
    # ------------------------
    plt.figure(figsize=(15,5))
    plt.scatter(dates, ndvi_vals, s=50)
    plt.plot(dates, ndvi_vals, alpha=0.6)

    plt.xlabel("Date")
    plt.ylabel("NDVI")
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

plot_ndvi_chronological(
    layer_name="whole_kanha_table",
    title="NDVI Time Series - Whole Kanha (2010-2011)"
)
