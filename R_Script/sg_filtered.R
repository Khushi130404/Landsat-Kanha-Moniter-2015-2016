# ============================================================
# NDVI SMOOTHING USING SAVITZKY–GOLAY FILTER (ROBUST)
# ============================================================

# -----------------------------
# Load / install packages
# -----------------------------
pkgs <- c("dplyr", "lubridate", "zoo", "signal", "ggplot2")
for(p in pkgs){
  if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# -----------------------------
# Input / Output folders
# -----------------------------
input_dir  <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_7/data_interpolated"
output_dir <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_7/data_sg_smoothed"
plot_dir   <- "D:/Landsat_Kanha_Moniter_2015_2016/image/Landsat_7/plot_sg_smoothed"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# List CSV files
# -----------------------------
csv_files <- list.files(input_dir, pattern="\\.csv$", full.names=TRUE)
if(length(csv_files) == 0) stop("❌ No CSV files found")

# -----------------------------
# Loop through files
# -----------------------------
for(file in csv_files){
  
  cat("Processing:", basename(file), "\n")
  
  ndvi <- read.csv(file, stringsAsFactors = FALSE)
  
  # ---- Date parsing
  ndvi$date <- as.Date(ndvi$date, tryFormats = c("%Y-%m-%d", "%d-%m-%Y"))
  ndvi <- ndvi[order(ndvi$date), ]
  
  # ---- Detect NDVI column automatically
  ndvi_col <- intersect(c("ndvi_interp", "median_ndvi"), names(ndvi))[1]
  if(is.na(ndvi_col)){
    cat("⚠️ No NDVI column found, skipping\n\n")
    next
  }
  
  # ---- Ensure numeric
  ndvi[[ndvi_col]] <- as.numeric(ndvi[[ndvi_col]])
  
  # ---- Fill gaps BEFORE SG
  ndvi$ndvi_filled <- zoo::na.approx(
    ndvi[[ndvi_col]],
    na.rm = FALSE
  )
  
  # ---- Edge filling (important for SG)
  ndvi$ndvi_filled <- zoo::na.locf(ndvi$ndvi_filled, na.rm=FALSE)
  ndvi$ndvi_filled <- zoo::na.locf(ndvi$ndvi_filled, fromLast=TRUE)
  
  # ---- Apply Savitzky–Golay filter
  window_size <- 9   
  poly_order  <- 2
  
  ndvi$ndvi_sg <- signal::sgolayfilt(
    ndvi$ndvi_filled,
    p = poly_order,
    n = window_size
  )
  
  # ---- Save CSV
  out_csv <- paste0(tools::file_path_sans_ext(basename(file)), "_SG.csv")
  write.csv(ndvi, file.path(output_dir, out_csv), row.names = FALSE)
  cat("✅ Saved CSV:", out_csv, "\n")
  
  # ---- Plot
  p <- ggplot(ndvi, aes(x=date)) +
    geom_point(aes(y=ndvi[[ndvi_col]]), color="red", size=1.4) +
    geom_line(aes(y=ndvi_filled), color="grey60", linewidth=0.6) +
    geom_line(aes(y=ndvi_sg), color="blue", linewidth=1.2) +
    scale_x_date(date_breaks="1 month", date_labels="%b (%Y)") +
    labs(
      title = "NDVI Savitzky–Golay Smoothed Curve",
      subtitle = tools::file_path_sans_ext(basename(file)),
      x = "Time",
      y = "NDVI"
    ) +
    theme_minimal(base_size=13) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  out_plot <- paste0(tools::file_path_sans_ext(basename(file)), "_SG.png")
  ggsave(file.path(plot_dir, out_plot), p, width=10, height=5, dpi=300)
  
  cat("✅ Saved plot:", out_plot, "\n\n")
}

cat("🎉 SG smoothing completed successfully for all regions.\n")