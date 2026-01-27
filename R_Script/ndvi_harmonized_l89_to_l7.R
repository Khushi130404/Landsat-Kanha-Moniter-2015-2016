# ===============================
# Landsat 8/9 → Landsat 7 NDVI Harmonization
# + Monthly NDVI Plot Saving
# ===============================

library(dplyr)
library(lubridate)
library(readr)
library(purrr)
library(ggplot2)

# -------- PATHS --------
l7_dir  <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_7/data_sg_smoothed"
l89_dir <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_8_9/data_sg_smoothed"

out_data_dir <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_8_9/data_sg_scaled_to_L7"
out_plot_dir <- "D:/Landsat_Kanha_Moniter_2015_2016/image/Landsat_8_9/plot_harmonization"

dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_plot_dir, showWarnings = FALSE, recursive = TRUE)

# -------- PARAMETERS --------
max_day_diff <- 5

# -------- FILE LIST --------
l7_files  <- list.files(l7_dir,  pattern = "\\.csv$", full.names = TRUE)
l89_files <- list.files(l89_dir, pattern = "\\.csv$", full.names = TRUE)

common_names <- intersect(basename(l7_files), basename(l89_files))

# -------- NEAREST DATE MATCH FUNCTION --------
nearest_match <- function(d1, d2) {
  which.min(abs(difftime(d2, d1, units = "days")))
}

# -------- PROCESS FILES --------
results <- map(common_names, function(fname) {
  
  cat("Processing:", fname, "\n")
  
  # ---------- REGION NAME ----------
  # Example:
  # se_kanha_table_5day_interpolated_SG.csv
  # → se_kanha
  region_name <- paste(strsplit(fname, "_")[[1]][1:2], collapse = "_")
  
  # ---------- READ DATA ----------
  l7  <- read_csv(file.path(l7_dir, fname), show_col_types = FALSE)
  l89 <- read_csv(file.path(l89_dir, fname), show_col_types = FALSE)
  
  l7$date  <- as.Date(l7$date)
  l89$date <- as.Date(l89$date)
  
  l7  <- l7  %>% select(date, ndvi_sg)
  l89 <- l89 %>% select(date, ndvi_sg)
  
  # ---------- DATE MATCHING ----------
  matched <- l7 %>%
    rowwise() %>%
    mutate(
      idx_l89 = nearest_match(date, l89$date),
      date_l89 = l89$date[idx_l89],
      day_diff = abs(as.numeric(date - date_l89)),
      ndvi_l89 = l89$ndvi_sg[idx_l89]
    ) %>%
    ungroup() %>%
    filter(day_diff <= max_day_diff) %>%
    rename(ndvi_l7 = ndvi_sg)
  
  if (nrow(matched) < 5) {
    warning(paste("Skipping (insufficient overlap):", fname))
    return(NULL)
  }
  
  # ---------- OLS REGRESSION ----------
  model <- lm(ndvi_l7 ~ ndvi_l89, data = matched)
  
  # ---------- APPLY MODEL ----------
  l89_scaled <- l89 %>%
    mutate(
      ndvi_sg_scaled = coef(model)[1] + coef(model)[2] * ndvi_sg,
      month = floor_date(date, "month")
    )
  
  l7_plot <- l7 %>%
    mutate(
      month = floor_date(date, "month")
    )
  
  # ---------- SAVE SCALED DATA ----------
  write_csv(
    l89_scaled,
    file.path(out_data_dir, fname)
  )
  
  # ---------- PLOT NDVI vs MONTH ----------
  p <- ggplot() +
    geom_line(data = l7_plot,
              aes(x = month, y = ndvi_sg, color = "Landsat 7"),
              linewidth = 0.9) +
    geom_line(data = l89_scaled,
              aes(x = month, y = ndvi_sg_scaled, color = "Landsat 8/9 (scaled)"),
              linewidth = 0.9) +
    scale_color_manual(values = c(
      "Landsat 7" = "blue",
      "Landsat 8/9 (scaled)" = "darkgreen"
    )) +
    labs(
      title = paste("NDVI Time Series –", region_name),
      x = "Month",
      y = "NDVI",
      color = "Sensor"
    ) +
    theme_minimal(base_size = 13)
  
  ggsave(
    filename = paste0(region_name, "_ndvi_harmonized.png"),
    plot = p,
    path = out_plot_dir,
    width = 8,
    height = 4.5,
    dpi = 300
  )
  
  # ---------- RETURN SUMMARY ----------
  tibble(
    region = region_name,
    intercept = coef(model)[1],
    slope = coef(model)[2],
    r2 = summary(model)$r.squared,
    n_pairs = nrow(matched)
  )
})

# -------- SAVE REGRESSION SUMMARY --------
reg_summary <- bind_rows(results)

write_csv(
  reg_summary,
  file.path(out_plot_dir, "landsat_harmonization_summary.csv")
)

cat("ALL REGIONS PROCESSED SUCCESSFULLY ✅\n")
