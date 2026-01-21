# ============================================================
# NDVI-SG TIME SERIES COMPARISON (DAY PRESERVED)
# Landsat-7 vs Landsat-8/9
# ============================================================

# -----------------------------
# 0. Install missing packages
# -----------------------------
packages <- c("dplyr","lubridate","ggplot2","readr")
for(p in packages){
  if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
}

# -----------------------------
# 1. Load libraries
# -----------------------------
library(dplyr)
library(lubridate)
library(ggplot2)
library(readr)

# -----------------------------
# 2. Input / Output folders
# -----------------------------
folder_1 <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_7/data_sg_smoothed"
folder_2 <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_8_9/data_sg_smoothed"

plot_dir <- "D:/Landsat_Kanha_Moniter_2015_2016/image/Comparision/plot_sg_smoothed"
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# -----------------------------
# 3. List CSV files
# -----------------------------
csv_files <- list.files(folder_1, pattern="\\.csv$", full.names = FALSE)
if(length(csv_files) == 0) stop("❌ No CSV files found!")

# -----------------------------
# 4. Process each matching file
# -----------------------------
for(fname in csv_files){

  file1 <- file.path(folder_1, fname)
  file2 <- file.path(folder_2, fname)

  if(!file.exists(file2)){
    cat("⚠️ Missing in Landsat 8/9:", fname, "\n")
    next
  }

  cat("Processing:", fname, "\n")

  # -----------------------------
  # Read CSVs
  # -----------------------------
  l7  <- read_csv(file1, show_col_types = FALSE)
  l89 <- read_csv(file2, show_col_types = FALSE)

  # -----------------------------
  # Parse date (KEEP day info)
  # -----------------------------
  l7$date  <- as.Date(l7$date)
  l89$date <- as.Date(l89$date)

  # -----------------------------
  # Clean & prepare (USE ndvi_sg)
  # -----------------------------
  l7 <- l7 %>%
    filter(!is.na(ndvi_sg)) %>%
    arrange(date) %>%
    mutate(sensor = "Landsat 7")

  l89 <- l89 %>%
    filter(!is.na(ndvi_sg)) %>%
    arrange(date) %>%
    mutate(sensor = "Landsat 8/9")

  data_all <- bind_rows(l7, l89)

  if(nrow(data_all) < 4){
    cat("⚠️ Not enough data, skipping:", fname, "\n")
    next
  }

  # -----------------------------
  # Plot: NDVI-SG vs TIME
  # -----------------------------
  p <- ggplot(
    data_all,
    aes(x = date, y = ndvi_sg,
        color = sensor, group = sensor)
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%b"
    ) +
    scale_color_manual(
      values = c("Landsat 7" = "blue",
                 "Landsat 8/9" = "red")
    ) +
    labs(
      title = "SG-Filtered NDVI Time Series Comparison",
      subtitle = tools::file_path_sans_ext(fname),
      x = "Month",
      y = "NDVI (Savitzky–Golay)",
      color = "Sensor"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # -----------------------------
  # Save plot
  # -----------------------------
  out_plot <- paste0(
    tools::file_path_sans_ext(fname),
    "_NDVI_SG_L7_vs_L89.png"
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

cat("🎉 All SG-filtered NDVI comparison plots generated.\n")
