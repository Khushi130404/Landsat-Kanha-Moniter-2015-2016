# ===================================================
# PLOT NDVI CURVES FOR ALL CSV FILES IN A FOLDER
# ===================================================

library(ggplot2)
library(lubridate)

# -----------------------------
# INPUT / OUTPUT FOLDERS
# -----------------------------
input_dir <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_8_9/data_interpolated"
plot_dir  <- "D:/Landsat_Kanha_Moniter_2015_2016/image/Landsat_8_9/plot_interpolation"

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# -----------------------------
# LIST CSV FILES
# -----------------------------
csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) stop("❌ No CSV files found")

# -----------------------------
# LOOP THROUGH FILES
# -----------------------------
for (file in csv_files) {
  
  cat("Plotting:", basename(file), "\n")
  
  df <- read.csv(file, stringsAsFactors = FALSE)
  df$date <- as.Date(df$date, tryFormats = c("%d-%m-%Y", "%Y-%m-%d"))
  
  # Detect NDVI column automatically
  ndvi_col <- intersect(
    c("ndvi_sg", "ndvi_interp", "median_ndvi"),
    names(df)
  )[1]
  
  if (is.na(ndvi_col)) {
    cat("⚠️ No NDVI column found, skipping:", basename(file), "\n")
    next
  }
  
  # Plot
  p <- ggplot(df, aes(x = date, y = .data[[ndvi_col]])) +
    geom_line(color = "blue", linewidth = 1) +
    geom_point(color = "red", size = 1.2) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b (%Y)") +
    labs(
      title = "NDVI Time Series",
      subtitle = tools::file_path_sans_ext(basename(file)),
      x = "Time",
      y = "NDVI"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save plot
  out_plot <- paste0(
    tools::file_path_sans_ext(basename(file)),
    "_NDVI.png"
  )
  
  ggsave(
    filename = file.path(plot_dir, out_plot),
    plot = p,
    width = 10,
    height = 5,
    dpi = 300
  )
  
  cat("✅ Saved:", out_plot, "\n\n")
}

cat("🎉 All NDVI plots generated successfully.\n")