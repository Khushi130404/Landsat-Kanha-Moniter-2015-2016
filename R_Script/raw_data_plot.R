# ============================================================
# NDVI CONNECTED LINE PLOT (NO INTERPOLATION)
# ============================================================

# -----------------------------
# 0. Install missing packages
# -----------------------------
packages <- c("dplyr","lubridate","ggplot2")
for(p in packages){
  if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
}

# -----------------------------
# 1. Load libraries
# -----------------------------
library(dplyr)
library(lubridate)
library(ggplot2)

# -----------------------------
# 2. Input / Output folders
# -----------------------------
input_dir <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_7/data_raw"
plot_dir  <- "D:/Landsat_Kanha_Moniter_2015_2016/image/Landsat_7/plot_raw"

if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive=TRUE)

# -----------------------------
# 3. List CSV files
# -----------------------------
csv_files <- list.files(input_dir, pattern="\\.csv$", full.names=TRUE)
if(length(csv_files)==0) stop("❌ No CSV files found!")

# -----------------------------
# 4. Process each CSV
# -----------------------------
for(file in csv_files){
  
  cat("Processing:", basename(file), "\n")
  
  ndvi <- read.csv(file, stringsAsFactors=FALSE)
  
  # Convert date
  ndvi$date <- as.Date(ndvi$date, format="%d-%m-%Y")
  
  # Remove NA and sort
  ndvi <- ndvi %>%
    filter(!is.na(median_ndvi)) %>%
    arrange(date)
  
  if(nrow(ndvi) < 2){
    cat("⚠️ Not enough points, skipping:", basename(file), "\n")
    next
  }
  
  # -----------------------------
  # Plot: ONLY CONNECTED POINTS
  # -----------------------------
  p <- ggplot(ndvi, aes(x = date, y = median_ndvi)) +
    geom_line(color = "blue", linewidth = 0.8) +
    geom_point(color = "red", size = 2) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b-%Y") +
    labs(
      title = "NDVI Time Series (Observed Data Only)",
      subtitle = tools::file_path_sans_ext(basename(file)),
      x = "Date",
      y = "NDVI"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save plot
  out_plot <- paste0(
    tools::file_path_sans_ext(basename(file)),
    "_NDVI_connected.png"
  )
  
  ggsave(
    filename = file.path(plot_dir, out_plot),
    plot = p,
    width = 10,
    height = 5,
    dpi = 300
  )
  
  cat("✅ Saved plot:", out_plot, "\n\n")
}

cat("🎉 All CSV files plotted with connected lines only.\n")