# ============================================================
# NDVI SHAPE-PRESERVING INTERPOLATION (5-DAY INTERVAL)
# ============================================================

# -----------------------------
# 0. Install missing packages
# -----------------------------
packages <- c("dplyr","lubridate","ggplot2","zoo")
for(p in packages){
  if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
}

# -----------------------------
# 1. Load libraries
# -----------------------------
library(dplyr)
library(lubridate)
library(ggplot2)
library(zoo)

# -----------------------------
# 2. Input / Output folders
# -----------------------------
input_dir  <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_8_9/data_raw"
output_dir <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_8_9/data_interpolated"
plot_dir   <- "D:/Landsat_Kanha_Moniter_2015_2016/image/Landsat_8_9/plot_double_logistic"

if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)
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
  ndvi$date <- as.Date(ndvi$date, format="%d-%m-%Y")
  
  # Remove NA and sort
  ndvi <- ndvi %>% dplyr::filter(!is.na(median_ndvi)) %>% dplyr::arrange(date)
  if(nrow(ndvi)<2){
    cat("⚠️ Not enough points, skipping:", basename(file), "\n")
    next
  }
  
  # Create 5-day sequence for full period
  start_date <- min(ndvi$date)
  end_date   <- max(ndvi$date)
  date_seq   <- seq(start_date, end_date, by="5 days")
  
  # Shape-preserving interpolation
  ndvi_interp <- zoo::na.spline(ndvi$median_ndvi, x=as.numeric(ndvi$date),
                                xout=as.numeric(date_seq), method="monoH.FC")
  
  # Fill landsat by last observation carried forward
  landsat_interp <- zoo::na.locf(ndvi$landsat[match(ndvi$date, ndvi$date)], na.rm=FALSE)
  landsat_full <- rep(landsat_interp[length(landsat_interp)], length(date_seq))
  
  # Create output dataframe
  df_out <- data.frame(
    date = date_seq,
    median_ndvi = ndvi_interp,
    landsat = landsat_full
  )
  df_out$year  <- year(df_out$date)
  df_out$month <- month(df_out$date)
  df_out$day   <- day(df_out$date)
  
  # Save CSV
  out_csv <- paste0(tools::file_path_sans_ext(basename(file)), "_5day_interpolated.csv")
  write.csv(df_out, file.path(output_dir, out_csv), row.names=FALSE)
  cat("✅ Saved CSV:", out_csv, "\n")
  
  # Plot
  p <- ggplot() +
    geom_point(data=ndvi, aes(x=date, y=median_ndvi), color="red", size=2) +
    geom_line(data=df_out, aes(x=date, y=median_ndvi), color="blue", size=1) +
    scale_x_date(date_breaks="1 month", date_labels="%b (%Y)") +
    labs(title="NDVI Shape-Preserving Interpolation",
         subtitle=tools::file_path_sans_ext(basename(file)),
         x="Date (Apr 2015 → Apr 2016)",
         y="NDVI") +
    theme_minimal(base_size=13) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  
  # Save plot
  out_plot <- paste0(tools::file_path_sans_ext(basename(file)), "_NDVI_curve.png")
  ggsave(filename=file.path(plot_dir, out_plot), plot=p, width=10, height=5, dpi=300)
  cat("✅ Saved plot:", out_plot, "\n\n")
}

cat("🎉 All regions processed successfully.\n")