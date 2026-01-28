# ===============================
# Landsat 8/9 → Landsat 7 NDVI Harmonization
# Robust Version (No filter() error)
# ===============================

library(dplyr)
library(lubridate)
library(readr)
library(purrr)
library(ggplot2)
library(tools)

# -------- PATHS --------
l7_dir  <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_7/data_sg_smoothed"
l89_dir <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_8_9/data_sg_smoothed"

out_data_dir <- "D:/Landsat_Kanha_Moniter_2015_2016/Data_Table/Landsat_8_9/data_sg_scaled_to_L7"

plot_1curve_dir <- "D:/Landsat_Kanha_Moniter_2015_2016/image/Landsat_8_9/plot_harmonization"
plot_3curve_dir <- "D:/Landsat_Kanha_Moniter_2015_2016/image/Comparision/plot_harmonization"

summary_csv <- "D:/Landsat_Kanha_Moniter_2015_2016/landsat_harmonization_summary.csv"

dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_1curve_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_3curve_dir, recursive = TRUE, showWarnings = FALSE)

# -------- PARAMETERS --------
max_day_diff <- 5

# -------- FILE LIST --------
l7_files  <- list.files(l7_dir,  pattern = "\\.csv$", full.names = TRUE)
l89_files <- list.files(l89_dir, pattern = "\\.csv$", full.names = TRUE)

common_names <- intersect(basename(l7_files), basename(l89_files))

# -------- SAFE NEAREST MATCH FUNCTION --------
nearest_match <- function(d1, d2) {
  if (length(d2) == 0 || all(is.na(d2))) return(NA_integer_)
  which.min(abs(as.numeric(difftime(d2, d1, units = "days"))))
}

# -------- PROCESS FILES --------
results <- map(common_names, function(fname) {

  cat("Processing:", fname, "\n")

  # ---------- REGION NAME ----------
  region_name <- paste(strsplit(fname, "_")[[1]][1:2], collapse = "_")

  # ---------- READ DATA ----------
  l7  <- read_csv(file.path(l7_dir, fname), show_col_types = FALSE)
  l89 <- read_csv(file.path(l89_dir, fname), show_col_types = FALSE)

  if (nrow(l7) == 0 || nrow(l89) == 0) {
    warning("Skipping (empty file): ", fname)
    return(NULL)
  }

  l7$date  <- as.Date(l7$date)
  l89$date <- as.Date(l89$date)

  l7  <- l7  %>% select(date, ndvi_sg)
  l89 <- l89 %>% select(date, ndvi_sg)

  # ---------- DATE MATCHING ----------
  matched <- l7 %>%
    rowwise() %>%
    mutate(
      idx_l89  = nearest_match(date, l89$date),
      date_l89 = if (!is.na(idx_l89)) l89$date[idx_l89] else NA_Date_,
      ndvi_l89 = if (!is.na(idx_l89)) l89$ndvi_sg[idx_l89] else NA_real_,
      day_diff = abs(as.numeric(date - date_l89))
    ) %>%
    ungroup() %>%
    dplyr::filter(!is.na(ndvi_l89), day_diff <= max_day_diff) %>%
    rename(ndvi_l7 = ndvi_sg)

  if (nrow(matched) < 5) {
    warning("Skipping (insufficient overlap): ", fname)
    return(NULL)
  }

  # ---------- OLS REGRESSION ----------
  model <- lm(ndvi_l7 ~ ndvi_l89, data = matched)

  # ---------- APPLY MODEL ----------
  l89_scaled <- l89 %>%
    mutate(
      ndvi_sg_scaled = coef(model)[1] + coef(model)[2] * ndvi_sg
    )

  # ---------- SAVE SCALED DATA ----------
  write_csv(
    l89_scaled,
    file.path(out_data_dir, fname)
  )

  # ============================================================
  # 3-CURVE PLOT (L7 + L8/9 RAW + HARMONIZED)
  # ============================================================
  plot_3curve <- bind_rows(
    l7 %>% mutate(sensor = "Landsat 7", ndvi = ndvi_sg),
    l89 %>% mutate(sensor = "Landsat 8/9 (raw)", ndvi = ndvi_sg),
    l89_scaled %>% mutate(sensor = "Landsat 8/9 (harmonized)", ndvi = ndvi_sg_scaled)
  )

  p3 <- ggplot(plot_3curve, aes(date, ndvi, color = sensor)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b (%Y)") +
    labs(
      title = "NDVI Time Series (Landsat Harmonization)",
      subtitle = file_path_sans_ext(fname),
      x = "Time",
      y = "NDVI",
      color = "Sensor"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    filename = paste0(file_path_sans_ext(fname), "_3curve.png"),
    plot     = p3,
    path     = plot_3curve_dir,
    width    = 9,
    height   = 4.8,
    dpi      = 300
  )
  

  # ============================================================
  # 1-CURVE PLOT (L8/9 RAW vs HARMONIZED)
  # ============================================================
  plot_1curve <- bind_rows(
    l89_scaled %>% mutate(sensor = "Landsat 8/9 (harmonized)", ndvi = ndvi_sg_scaled)
  )

  p1 <- ggplot(plot_1curve, aes(date, ndvi, color = sensor)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b (%Y)") +
    labs(
      title = "NDVI Time Series (Landsat 8/9)",
      subtitle = file_path_sans_ext(fname),
      x = "Time",
      y = "NDVI",
      color = "Sensor"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    filename = paste0(file_path_sans_ext(fname), "_L89_only.png"),
    plot     = p1,
    path     = plot_1curve_dir,
    width    = 9,
    height   = 4.8,
    dpi      = 300
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

write_csv(reg_summary, summary_csv)

cat("ALL REGIONS DONE ✅\n")
