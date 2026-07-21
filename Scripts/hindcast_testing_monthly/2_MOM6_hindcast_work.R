library(ncdf4)
library(terra)
library(tidync)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(zoo)
library(matrixStats)

# ══════════════════════════════════════════════════════════════════════════════
# MHW DETECTION — SST AND BOTTOM TEMPERATURE
# Following Freeman et al. (2026) methodology
# DATA: MOM6 Northeast Pacific Hindcast 1993-2025
#
# IMPROVEMENTS OVER PREVIOUS VERSION:
#   - Bathymetry now available (deptho), good for depth mask
#   - 3D oxygen available — HLD/SHX now possible
#   - Bottom oxygen available
#   - Longer record: (1993-2025 vs 1995-2020)
#
# REMAINING LIMITATIONS:
#   - Monthly resolution (Freeman uses daily)
#   - No chlorophyll (LCX)
# ══════════════════════════════════════════════════════════════════════════════

data_dir <- "D:/LorenzoRFM_ Hollings/Downloads"

# File paths
sst_file    <- file.path(data_dir, "tos.nep.full.hcast.monthly.raw.r20250912.199301-202506.nc")
bT_file     <- file.path(data_dir, "tob.nep.full.hcast.monthly.raw.r20250912.199301-202506.nc")
bO2_file    <- file.path(data_dir, "btm_o2.nep.full.hcast.monthly.raw.r20250912.199301-202506.nc")
o2_3d_file  <- file.path(data_dir, "o2.nep.full.hcast.monthly.raw.r20250912.199301-202506.nc")
grid_file   <- file.path(data_dir, "ocean_static.nc")

# explore grid -----------------------------------------------------------------

nc_grid <- nc_open(grid_file)
names(nc_grid$var)
# "deptho" geolat"       "geolat_c"     "geolat_u"     "geolat_v"     "geolon"       "geolon_c"   "geolon_u"     "geolon_v" 
names(nc_grid$dim)
# jq , iq, jh , ih


#load in coords and bathymetry
lon2d <- ncvar_get(nc_grid, "geolon")   # 2D lon array
lat2d <- ncvar_get(nc_grid, "geolat")   # 2D lat array
depth <- ncvar_get(nc_grid, "deptho")   # bottom depth

nc_close(nc_grid)

cat("Grid size:", dim(lon2d), "\n")
# 342 816
cat("Lon range:", range(lon2d, na.rm = TRUE), "\n")
# ~ 156.9248 254.9719
cat("Lat range:", range(lat2d, na.rm = TRUE), "\n")
# 10.80904 - 80.71795
cat("Depth range:", range(depth, na.rm = TRUE), "\n")
# 5 - 6500

# explore one data file --------------------------------------------------------

nc_sst <- nc_open(sst_file)
names(nc_sst$var)
#[1] "average_DT" "average_T1" "average_T2" "time_bnds"  "tos"       
names(nc_sst$dim)
# time, ih, jh, nv

# check time
time_raw <- ncvar_get(nc_sst, "time")
time_units <- nc_sst$dim$time$units
cat("Time units:", time_units, "\n")
# days since 1993-01-01 00:00:00 
cat("N timesteps:", length(time_raw), "\n")   # = 390
cat("SST dims:", dim(ncvar_get(nc_sst, "tos", start=c(1,1,1), count=c(-1,-1,1))), "\n")
# 342 816

nc_close(nc_sst)

# parse time -------------------------------------------------------------------

# "days since 1993-01-01 00:00:00"
parse_time_mom6 <- function(file, timevar = "time") {
  nc         <- nc_open(file)
  raw        <- ncvar_get(nc, timevar)
  time_units <- nc$dim[[timevar]]$units
  nc_close(nc)
  
  origin <- as.Date(sub("days since ", "", time_units))
  return(origin + as.numeric(raw))
}

time_vec <- parse_time_mom6(sst_file)
cat("Time range:", as.character(range(time_vec)), "\n")  # 1993-01-16 to 2025-06-16
cat("Total timesteps:", length(time_vec), "\n")           # 390

# find test point --------------------------------------------------------------

# MOM6 uses 0-360 longitude convention
# Convert target from -180:180 to 0:360
target_lon <- -122 + 360   # 231
target_lat <- 36.7

ocean_mask <- !is.na(depth) & depth > 0
dist <- sqrt((lon2d - target_lon)^2 + (lat2d - target_lat)^2)
dist[!ocean_mask] <- NA

idx <- which(dist == min(dist, na.rm = TRUE), arr.ind = TRUE)[1, ]
ilon <- as.integer(idx[1])
ilat <- as.integer(idx[2])
stopifnot(ilon >= 1, ilon <= dim(lon2d)[1])
stopifnot(ilat >= 1, ilat <= dim(lon2d)[2])

cat("Test point: ih =", ilon, ", jh =", ilat, "\n") # ih = 232 , jh = 235
cat("Actual lon:", lon2d[ilon, ilat], "(=", lon2d[ilon, ilat] - 360, "in -180:180)\n")
# actual lon: 237.9789 (= -122.0211)
cat("Actual lat:", lat2d[ilon, ilat], "\n")
# actual lat: 36.67264
cat("Bottom depth:", depth[ilon, ilat], "m\n")
# bottom depth: 312.7541  m
cat("Qualifies for bottom MHW (< 1000m)?", depth[ilon, ilat] < 1000, "\n")
# TRUE!

# ── extract time series ───────────────────────────────────────────────────────

# Confirm dimension order before extracting
# So start = c(ilon, ilat, 1), count = c(1, 1, -1)

extract_point_mom6 <- function(file, varname, ilon, ilat) {
  nc  <- nc_open(file)
  val <- ncvar_get(nc, varname,
                   start = c(ilon, ilat, 1),
                   count = c(1, 1, -1))
  nc_close(nc)
  return(as.numeric(val))
}

sst_ts  <- extract_point_mom6(sst_file, "tos", ilon, ilat)
bT_ts   <- extract_point_mom6(bT_file, "tob", ilon, ilat)
bO2_ts  <- extract_point_mom6(bO2_file, "btm_o2", ilon, ilat)

cat("SST  — length:", length(sst_ts),  "| NAs:", sum(is.na(sst_ts)),
    "| range:", round(range(sst_ts,  na.rm=TRUE), 2), "\n")
# SST: length = 390 | NAs: 0 | range: 11.05 - 18.32
cat("bT   — length:", length(bT_ts),   "| NAs:", sum(is.na(bT_ts)),
    "| range:", round(range(bT_ts,   na.rm=TRUE), 2), "\n")
# bT   — length: 390 | NAs: 0 | range: 6.43 8.03
cat("bO2  — length:", length(bO2_ts),  "| NAs:", sum(is.na(bO2_ts)),
    "| range:", round(range(bO2_ts,  na.rm=TRUE), 2), "\n")
# bO2  — length: 390 | NAs: 0 | range: 0 0 

par(mfrow = c(3, 1), mar = c(2, 4, 2, 1))
plot(time_vec, sst_ts,  type = "l", main = "SST",        ylab = "°C")
plot(time_vec, bT_ts,   type = "l", main = "Bottom Temp", ylab = "°C")
plot(time_vec, bO2_ts,  type = "l", main = "Bottom O2",   ylab = "mmol/m3")
par(mfrow = c(1, 1))

# Shared functions -------------------------------------------------------------

# Seasonally-varying percentila threshold (Freeman Section 2.3)
# Freeman: 11 day winddow centered on each day of year x 24 years = 264 values
# ADAPTATION: 3-month window centered on each month x 26 years = 78 values
# Reason: monthly resolution would make sub-monthly windows meaningless

# 3 month window centered on each month, pooled across all 26 yrs
compute_threshold <- function(ts, months, percentile = 0.90, window = 1) {
  thresh <- numeric(12)
  for (m in 1:12) {
    window_months <- ((m - 1 - window):(m - 1 + window)) %% 12 + 1
    idx <- which(months %in% window_months)
    thresh[m] <- quantile(ts[idx], probs = percentile, na.rm = TRUE)
  }
  return(thresh)
}

compute_threshold_matrix <- function(ts_mat, months, percentile = 0.90, window = 1) {
  n_cells <- nrow(ts_mat)
  thresh_mat <- matrix(NA_real_, n_cells, 12)
  for (m in 1:12) {
    window_months <- ((m - 1 - window):(m - 1 + window)) %% 12 + 1
    idx <- which(months %in% window_months)
    thresh_mat[, m] <- rowQuantiles(ts_mat[, idx, drop = FALSE],
                             probs = percentile, na.rm = TRUE)
  }
  thresh_mat
}


# Merge consecutive events within max_gap timesteps (Freeman Secion 2.3)
# Freeman: merge if gap <= 3 days
# ADAPTATION: merge if gap <= 1 month

merge_events <- function(exceeds, max_gap = 1) {
  extreme_idx <- which(exceeds == 1)
  if (length(extreme_idx) == 0) return(exceeds)
  
  gaps           <- diff(extreme_idx)
  fill_positions <- which(gaps > 1 & gaps <= max_gap + 1)
  
  for (pos in fill_positions) {
    start_fill <- extreme_idx[pos] + 1
    end_fill   <- extreme_idx[pos + 1] - 1
    if (start_fill <= end_fill) {  # protects against zero-length fills
      exceeds[start_fill:end_fill] <- 1
    }
  }
  return(exceeds)
}


# Remove events shorter than minimum duration (Freeman et al. Section 2.3)
# Freeman: minimum 5 days
# Adaptation: Minimum 2 months

remove_short_events <- function(exceeds, min_duration = 2) {
  r <- rle(exceeds)
  r$values[r$values == 1 & r$lengths < min_duration] <- 0
  return(inverse.rle(r))
}


# Climatological z-score anomaly (Freeman et al. Section 2.3)
# Anomaly computed relative to monthly climatological mean and sd
# NOT the overall time series mean. removes seasonal cycle
compute_anomaly <- function(ts, months) {
  clim_mean <- tapply(ts, months, mean, na.rm = TRUE)
  clim_sd   <- tapply(ts, months, sd,   na.rm = TRUE)
  (ts - clim_mean[months]) / clim_sd[months]
}


# Duration, intensity, frequency (Freeman Section 2.4)
characterize_events <- function(events_flag, anomaly, time_vec) {
  r      <- rle(events_flag)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  
  
  
  out <- data.frame()
  for (i in seq_along(r$values)) {
    if (r$values[i] == 1) {
      idx <- starts[i]:ends[i]
      out <- rbind(out, data.frame(
        start           = time_vec[starts[i]],
        end             = time_vec[ends[i]],
        duration_months = length(idx),           # Freeman: number of days
        mean_intensity  = round(mean(anomaly[idx], na.rm = TRUE), 3)  # Freeman: mean z-score anomaly
      ))
    }
  }
  # Frequency = number of events over study period (Freeman Section 2.4)
  # Reported as events per year at this location
  n_years <- as.numeric(diff(range(time_vec))) / 365.25
  location_freq <- round(nrow(out) / n_years, 2)
  cat("Event frequency at this location:", location_freq, "events/year\n")
  return(out)
  
}


# ══════════════════════════════════════════════════════════════════════════════
# PART 1a: SURFACE MHW (SMHW) — SST ---------------------------------------------
# Freeman Section 2.3: >= 5 days exceeding P90
# No depth restriction for surface MHW
# ══════════════════════════════════════════════════════════════════════════════

df <- data.frame(
  time  = time_vec,
  sst   = sst_ts,
  month = as.integer(format(time_vec, "%m"))
)

## Step 1: P90 threshold --------------------------------------------------------
thresh_90 <- compute_threshold(df$sst, df$month, percentile = 0.90, window = 1)
df$thresh_90 <- thresh_90[df$month]
print(data.frame(month = 1:12, SST_P90 = round(thresh_90, 2)))
#    month SST_P90
# 1      1   13.55
# 2      2   13.74
# 3      3   14.26
# 4      4   14.64
# 5      5   15.30
# 6      6   16.31
# 7      7   16.94
# 8      8   17.16
# 9      9   17.16
# 10    10   16.74
# 11    11   15.68
# 12    12   14.28

## Step 2: Flag exceedances ----------------------------------------------------
df$exceeds <- ifelse(df$sst > df$thresh_90, 1, 0)
cat("SST months exceeding P90:", sum(df$exceeds), "of", nrow(df), "\n")
# 28 of 390

## Step 3: Merge close events ----------------------------------------------------
df$exceeds_merged <- merge_events(df$exceeds, max_gap = 1)

## Step 4: Remove short events ----------------------------------------------------
df$events <- remove_short_events(df$exceeds_merged, min_duration = 2)
cat("SST months flagged as SMHW:", sum(df$events), "\n")
# 24 months 

## Step 5: Anomaly --------------------------------------------------------------
df$anomaly <- compute_anomaly(df$sst, df$month)

## Step 6: Characterize --------------------------------------------------------
sst_events <- characterize_events(df$events, df$anomaly, df$time)

# event freq at this location: 0.12 events/yr
print(sst_events)
#        start        end duration_months mean_intensity
# 1 1997-08-16 1998-03-16               8          2.017
# 2 2005-03-16 2005-05-16               3          1.247
# 3 2014-07-16 2015-03-16               9          2.211
# 4 2015-07-16 2015-10-16               4          2.444

## Step 7: Plot ----------------------------------------------------------------
ggplot(df, aes(x = time)) +
  geom_rect(
    data = df %>% filter(events == 1),
    aes(xmin = time - 15, xmax = time + 15, ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.5
  ) +
  geom_line(aes(y = sst),      color = "black", linewidth = 0.75) +
  geom_line(aes(y = thresh_90), color = "black", linetype = "dashed", linewidth = 0.5) +
  labs(
    title    = paste("SMHW | lon: -122.02",
                     "lat:", round(lat2d[ilon, ilat], 2)),
    subtitle = paste("Method from Freeman et al. | P90 threshold (dashed) | Red = SMHW\n",
                     "ADAPTATIONS: monthly res; 3-month window; 2-month min duration"),
    y = "SST (°C)", x = NULL
  ) 

ggplot(df, aes(x = time)) +
  geom_rect(
    data = df %>% filter(events == 1),
    aes(xmin = time - 15, xmax = time + 15, ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.5
  ) +
  geom_line(aes(y = sst),      color = "black", linewidth = 0.75) +
  geom_line(aes(y = thresh_90), color = "black", linetype = "dashed", linewidth = 0.5) +
  labs(
    title    = paste("SMHW | lon: -122.02",
                     "lat:", round(lat2d[ilon, ilat], 2)),
    subtitle = paste("Method from Freeman et al. | P90 threshold (dashed) | Red = SMHW\n"),
    y = "SST (°C)", x = NULL
  ) + theme_light()


# PART 1b: SMHW Full Grid ------------------------------------------------------
# Going chunk by chunk
# All functions already defined

# Grid dimensions (already known)
n_ih    <- 342
n_jh    <- 816
n_times <- 390
months  <- as.integer(format(time_vec, "%m"))  
n_years <- as.numeric(diff(range(time_vec))) / 365.25

## Step 1: Output arrays -------------------------------------------------------

# Summary stats per cell
out_n_events       <- matrix(NA_real_, n_ih, n_jh)
out_freq           <- matrix(NA_real_, n_ih, n_jh)
out_mean_duration  <- matrix(NA_real_, n_ih, n_jh)
out_mean_intensity <- matrix(NA_real_, n_ih, n_jh)

# Full event time series [ih x jh x time]
out_events_ts <- array(NA_integer_, dim = c(n_ih, n_jh, n_times))

## Step 2: Chunk loop ----------------------------------------------------------

nc_sst     <- nc_open(sst_file)
start_time <- Sys.time()

for (j in 1:n_jh) {
  
  # Progress every 50 rows
  if (j %% 50 == 0) {
    elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
    cat("Row", j, "of", n_jh, "| Elapsed:", elapsed, "mins\n")
  }
  
  # Skip if no ocean cells in this row
  if (!any(ocean_mask[, j], na.rm = TRUE)) next
  
  # Load entire row at once: [342 x 390]
  sst_row <- ncvar_get(nc_sst, "tos",
                       start = c(1, j, 1),
                       count = c(-1, 1, -1))
  
  for (i in 1:n_ih) {
    
    # Skip land
    if (!ocean_mask[i, j]) next
    
    # EXACT SAME STEPS AS SINGLE POINT CODE
    
    ts <- sst_row[i, ]
    if (all(is.na(ts))) next
    
    ### 1: P90 threshold --------------------------------------------------
    thresh_90 <- compute_threshold(ts, months, percentile = 0.90, window = 1)
    thresh_ts <- thresh_90[months]
    
    ### 2: Flag exceedances ----------------------------------------------------
    exceeds <- ifelse(ts > thresh_ts, 1, 0)
    exceeds[is.na(exceeds)] <- 0
    
    ### 3: Merge close events --------------------------------------------------
    exceeds_merged <- merge_events(exceeds, max_gap = 1)
    
    ### 4: Remove short events --------------------------------------------------------
    events <- remove_short_events(exceeds_merged, min_duration = 2)
    
    ### Step 5: Anomaly --------------------------------------------------------
    anomaly <- compute_anomaly(ts, months)
    
    ### 6: Characterize (inline version of characterize_events) -----------------
    r      <- rle(events)
    ends   <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    
    durations  <- c()
    intensities <- c()
    for (k in seq_along(r$values)) {
      if (r$values[k] == 1) {
        idx         <- starts[k]:ends[k]
        durations   <- c(durations,   length(idx))
        intensities <- c(intensities, mean(anomaly[idx], na.rm = TRUE))
      }
    }
    
    # Store results
    out_n_events[i, j]       <- length(durations)
    out_freq[i, j]           <- length(durations) / n_years
    out_mean_duration[i, j]  <- ifelse(length(durations) > 0,
                                       mean(durations), 0)
    out_mean_intensity[i, j] <- ifelse(length(intensities) > 0,
                                       mean(intensities), 0)
    out_events_ts[i, j, ]   <- events
  }
}

nc_close(nc_sst)
cat("Done! Total time:",
    round(difftime(Sys.time(), start_time, units = "mins"), 1), "mins\n")
# first round went 17 min

## Part 3. Save ----------------------------------------------------------------

# only run if the loop changes
# setwd("Outputs")
# 
# saveRDS(list(
#   n_events       = out_n_events,
#   freq           = out_freq,
#   mean_duration  = out_mean_duration,
#   mean_intensity = out_mean_intensity,
#   lon2d          = lon2d,
#   lat2d          = lat2d,
#   ocean_mask     = ocean_mask
# ), "sst_mhw_summary.rds")
# 
# saveRDS(list(
#   events_ts = out_events_ts,
#   time_vec  = time_vec,
#   lon2d     = lon2d,
#   lat2d     = lat2d
), "sst_mhw_events_ts.rds")

setwd("C:/Users/mcclo/OneDrive/Documents/Hollings_Internship")
getwd()

sst_summary <- readRDS(file.path(output_dir, "sst_mhw_summary.rds"))
out_n_events       <- sst_summary$n_events
out_freq           <- sst_summary$freq
out_mean_duration  <- sst_summary$mean_duration
out_mean_intensity <- sst_summary$mean_intensity
ocean_mask         <- sst_summary$ocean_mask

sst_ts_data   <- readRDS(file.path(output_dir, "sst_mhw_events_ts.rds"))
out_events_ts <- sst_ts_data$events_ts
time_vec      <- sst_ts_data$time_vec

## Part 4: Validate against single point ---------------------------------------

# Should exactly match your single point results
cat("Spatial run — freq at test point:", round(out_freq[ilon, ilat], 2), "\n")
cat("Single point run — freq:          0.12\n")
cat("N events at test point:", out_n_events[ilon, ilat], "(expect 4)\n")
# passed all tests!


# Part 1c: Domain-averaged SST with MHW periods ------------------------------------------

# Load full SST array [ih x jh x time]
nc_sst  <- nc_open(sst_file)
sst_all <- ncvar_get(nc_sst, "tos")   # [342 x 816 x 390]
nc_close(nc_sst)

# Mask land
for (t in 1:n_times) {
  sst_all[, , t][!ocean_mask] <- NA
}

# Domain mean SST at each timestep
sst_mean_ts <- apply(sst_all, 3, mean, na.rm = TRUE)

# Fraction of ocean cells in MHW at each timestep
# out_events_ts is [ih x jh x time] from spatial run
mhw_fraction <- apply(out_events_ts, 3, function(x) {
  mean(x[ocean_mask] == 1, na.rm = TRUE)
})

cat("MHW fraction range:", round(range(mhw_fraction, na.rm = TRUE), 3), "\n")
# 0 - 0.436

## Build DF -------------------------------------------------

df_spatial <- data.frame(
  time         = time_vec,
  sst_mean     = sst_mean_ts,
  mhw_fraction = mhw_fraction,
  month        = as.integer(format(time_vec, "%m"))
)

# Domain-mean P90 threshold
thresh_90_mean <- compute_threshold(df_spatial$sst_mean,
                                    df_spatial$month,
                                    percentile = 0.90,
                                    window = 1)
df_spatial$thresh_90 <- thresh_90_mean[df_spatial$month]

# Flag periods where >= 15% of ocean cells are in MHW
# (threshold can be adjusted needed)
df_spatial$widespread_mhw <- ifelse(df_spatial$mhw_fraction >= 0.15, 1, 0)

## Plot 1: Domain mean SST w widespread MHW shading ----------------------------

ggplot(df_spatial, aes(x = time)) +
  geom_rect(
    data = df_spatial %>% filter(widespread_mhw == 1),
    aes(xmin = time - 15, xmax = time + 15,
        ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.5
  ) +
  geom_line(aes(y = sst_mean),  color = "black", linewidth = 0.75) +
  geom_line(aes(y = thresh_90), color = "black",
            linetype = "dashed", linewidth = 0.5) +
  labs(
    title    = "Domain-Mean SST MHW | MOM6 Hindcast 1993-2025",
    subtitle = paste("Dashed = P90 | Red = >= 15% of ocean cells in a MHW\n",
                     "CHANGES: monthly res; 3-month window; 2-month min duration"),
    y = "SST (°C)", x = NULL
  ) 

## Plot 2: Fraction of cells in MHW over time ----------------------------------
# Shows how spatially extensive each event was

ggplot(df_spatial, aes(x = time, y = mhw_fraction * 100)) +
  geom_area(fill = "red", alpha = 0.4) +
  geom_line(color = "red", linewidth = 0.75) +
  labs(
    title    = "Spatial Extent of SST MHWs | MOM6 Hindcast 1993-2025",
    subtitle = "% of ocean cells experiencing MHW conditions each month",
    y = "% of domain in MHW", x = NULL
  )

## Plot 3: Both together -------------------------------------------------------

library(patchwork)   

p1 <- ggplot(df_spatial, aes(x = time)) +
  geom_rect(
    data = df_spatial %>% filter(widespread_mhw == 1),
    aes(xmin = time - 15, xmax = time + 15,
        ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.4
  ) +
  geom_line(aes(y = sst_mean),  color = "black", linewidth = 0.75) +
  geom_line(aes(y = thresh_90), color = "black",
            linetype = "dashed", linewidth = 0.5) +
  labs(title    = "Domain-Mean SST | Red = widespread MHW (>= 15% cells)",
       y = "SST (°C)", x = NULL) 

p2 <- ggplot(df_spatial, aes(x = time, y = mhw_fraction * 100)) +
  geom_area(fill = "red", alpha = 0.4) +
  geom_line(color = "darkred", linewidth = 0.75) +
  geom_hline(yintercept = 10, linetype = "dashed", linewidth = 0.5) +
  labs(y = "% domain in MHW", x = NULL) 

p1 / p2 



# ══════════════════════════════════════════════════════════════════════════════
# PART 2: BOTTOM MHW (BMHW) — bT ---------------------------------------------
# Freeman Section 2.3: P90, restricted to depth < 1000m
# Depth mask via ocean_static.nc
# ══════════════════════════════════════════════════════════════════════════════

# Confirm test point qualifies for bottom MHW (depth < 1000m)
cat("Bottom depth at test point:", depth[ilon, ilat], "m\n")
cat("Qualifies for bottom MHW:", depth[ilon, ilat] < 1000, "\n")
# TRUE

df_bT <- data.frame(
  time  = time_vec,
  bT    = bT_ts,
  month = as.integer(format(time_vec, "%m"))
)

## Step 1: P90 threshold --------------------------------------------------------
thresh_90_bT    <- compute_threshold(df_bT$bT, df_bT$month, percentile = 0.90, window = 1)
df_bT$thresh_90 <- thresh_90_bT[df_bT$month]
print(data.frame(month = 1:12, bT_P90 = round(thresh_90_bT, 2)))
#    month bT_P90
# 1      1   7.67
# 2      2   7.58
# 3      3   7.47
# 4      4   7.41
# 5      5   7.42
# 6      6   7.58
# 7      7   7.73
# 8      8   7.76
# 9      9   7.77
# 10    10   7.75
# 11    11   7.73
# 12    12   7.69

## Step 2: Flag exceedances  --------------------------------------------------------
df_bT$exceeds <- ifelse(df_bT$bT > df_bT$thresh_90, 1, 0)
cat("bT months exceeding P90:", sum(df_bT$exceeds), "of", nrow(df_bT), "\n")
# 34 out of 390

## Step 3: Merge close events   --------------------------------------------------------
df_bT$exceeds_merged <- merge_events(df_bT$exceeds, max_gap = 1)

## Step 4: Remove short events   --------------------------------------------------------
df_bT$events_bT <- remove_short_events(df_bT$exceeds_merged, min_duration = 2)
cat("bT months flagged as BMHW:", sum(df_bT$events_bT), "\n")
# 24 months flagged 

## Step 5: Anomaly   --------------------------------------------------------
df_bT$anomaly <- compute_anomaly(df_bT$bT, df_bT$month)

## Step 6: Characterize   --------------------------------------------------------
bT_events <- characterize_events(df_bT$events_bT, df_bT$anomaly, df_bT$time)
# freq: 0.19 events/yr
print(bT_events)
#        start        end duration_months mean_intensity
# 1 2008-08-16 2008-09-16               2          1.216
# 2 2009-12-16 2010-02-15               3          1.436
# 3 2014-05-16 2014-07-16               3          1.315
# 4 2015-08-16 2015-09-16               2          1.239
# 5 2019-11-16 2020-03-16               5          1.409
# 6 2023-04-16 2023-12-16               9          1.703


## Step 7: Plot --------------------------------------------------------
ggplot(df_bT, aes(x = time)) +
  geom_rect(
    data = df_bT %>% filter(events_bT == 1),
    aes(xmin = time - 15, xmax = time + 15, ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.5
  ) +
  geom_line(aes(y = bT),       color = "black", linewidth = 0.75) +
  geom_line(aes(y = thresh_90), color = "black", linetype = "dashed", linewidth = 0.5) +
  labs(
    title    = paste("BMHW | lon: -122.02",
                     "lat:", round(lat2d[ilon, ilat], 2),
                     "| depth:", round(depth[ilon, ilat], 0), "m"),
    subtitle = paste("Method from Freeman et al. | P90 threshold (dashed) | Red = BMHW\n",
                     "Depth mask applied (< 1000m) | monthly res; 3-month window; 2-month min duration"),
    y = "Bottom Temperature (°C)", x = NULL
  )

ggplot(df_bT, aes(x = time)) +
  geom_rect(
    data = df_bT %>% filter(events_bT == 1),
    aes(xmin = time - 15, xmax = time + 15, ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.5
  ) +
  geom_line(aes(y = bT),       color = "black", linewidth = 0.75) +
  geom_line(aes(y = thresh_90), color = "black", linetype = "dashed", linewidth = 0.5) +
  labs(
    title    = paste("BMHW | lon: -122.02",
                     "lat:", round(lat2d[ilon, ilat], 2),
                     "| depth:", round(depth[ilon, ilat], 0), "m"),
    subtitle = paste("Method from Freeman et al. | P90 threshold (dashed) | Red = BMHW\n"),
    y = "Bottom Temperature (°C)", x = NULL
  ) + theme_light()

# Grid dimensions (already known)
n_ih    <- 342
n_jh    <- 816
n_times <- 390
months  <- as.integer(format(time_vec, "%m"))  
n_years <- as.numeric(diff(range(time_vec))) / 365.25

# PART 2b: BMHW Full Grid ------------------------------------------------------

## Step 1: Output arrays -------------------------------------------------------

# Summary stats per cell
out_n_events_bT       <- matrix(NA_real_, n_ih, n_jh)
out_freq_bT           <- matrix(NA_real_, n_ih, n_jh)
out_mean_duration_bT  <- matrix(NA_real_, n_ih, n_jh)
out_mean_intensity_bT <- matrix(NA_real_, n_ih, n_jh)

# Full event time series [ih x jh x time]
out_events_ts_bT <- array(NA_integer_, dim = c(n_ih, n_jh, n_times))

bottom_mask <- ocean_mask & depth < 1000

## Step 2: Chunk loop ----------------------------------------------------------

# nc_bT     <- nc_open(bT_file)
# start_time <- Sys.time()
# 
# for (j in 1:n_jh) {
#   
#   # Progress every 50 rows
#   if (j %% 50 == 0) {
#     elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
#     cat("Row", j, "of", n_jh, "| Elapsed:", elapsed, "mins\n")
#   }
#   
#   # Skip if no ocean cells in this row
#   if (!any(bottom_mask[, j], na.rm = TRUE)) next
#   
#   # Load entire row at once: [342 x 390]
#   bT_row <- ncvar_get(nc_bT, "tob",
#                        start = c(1, j, 1),
#                        count = c(-1, 1, -1))
#   
#   for (i in 1:n_ih) {
#     
#     # Skip land
#     if (!bottom_mask[i, j]) next
#     
#     # EXACT SAME STEPS AS SINGLE POINT CODE
#     
#     ts <- bT_row[i, ]
#     if (all(is.na(ts))) next
#     
#     ### 1: P90 threshold --------------------------------------------------
#     thresh_90 <- compute_threshold(ts, months, percentile = 0.90, window = 1)
#     thresh_ts <- thresh_90[months]
#     
#     ### 2: Flag exceedances ----------------------------------------------------
#     exceeds <- ifelse(ts > thresh_ts, 1, 0)
#     exceeds[is.na(exceeds)] <- 0
#     
#     ### 3: Merge close events --------------------------------------------------
#     exceeds_merged <- merge_events(exceeds, max_gap = 1)
#     
#     ### 4: Remove short events --------------------------------------------------------
#     events <- remove_short_events(exceeds_merged, min_duration = 2)
#     
#     ### Step 5: Anomaly --------------------------------------------------------
#     anomaly <- compute_anomaly(ts, months)
#     
#     ### 6: Characterize (inline version of characterize_events) -----------------
#     r      <- rle(events)
#     ends   <- cumsum(r$lengths)
#     starts <- ends - r$lengths + 1
#     
#     durations  <- c()
#     intensities <- c()
#     for (k in seq_along(r$values)) {
#       if (r$values[k] == 1) {
#         idx         <- starts[k]:ends[k]
#         durations   <- c(durations,   length(idx))
#         intensities <- c(intensities, mean(anomaly[idx], na.rm = TRUE))
#       }
#     }
#     
#     # Store results
#     out_n_events_bT[i, j]       <- length(durations)
#     out_freq_bT[i, j]           <- length(durations) / n_years
#     out_mean_duration_bT[i, j]  <- ifelse(length(durations) > 0,
#                                        mean(durations), 0)
#     out_mean_intensity_bT[i, j] <- ifelse(length(intensities) > 0,
#                                        mean(intensities), 0)
#     out_events_ts_bT[i, j, ]   <- events
#   }
# }
# 
# nc_close(nc_bT)
# cat("Done! Total time:",
#     round(difftime(Sys.time(), start_time, units = "mins"), 1), "mins\n")
# first round went 17 min

## Part 3. Save ----------------------------------------------------------------

# setwd("Outputs")
# 
# saveRDS(list(
#   n_events       = out_n_events_bT,
#   freq           = out_freq_bT,
#   mean_duration  = out_mean_duration_bT,
#   mean_intensity = out_mean_intensity_bT,
#   lon2d          = lon2d,
#   lat2d          = lat2d,
#   bottom_mask     = bottom_mask
# ), "bT_mhw_summary.rds")
# 
# saveRDS(list(
#   events_ts = out_events_ts_bT,
#   time_vec  = time_vec,
#   lon2d     = lon2d,
#   lat2d     = lat2d
# ), "bT_mhw_events_ts.rds")
# 
# setwd("C:/Users/mcclo/OneDrive/Documents/Hollings_Internship")
# getwd()

bT_summary <- readRDS(file.path(output_dir, "bT_mhw_summary.rds"))
out_n_events_bT       <- bT_summary$n_events
out_freq_bT           <- bT_summary$freq
out_mean_duration_bT  <- bT_summary$mean_duration
out_mean_intensity_bT <- bT_summary$mean_intensity
bottom_mask           <- bT_summary$bottom_mask

bT_ts_data       <- readRDS(file.path(output_dir, "bT_mhw_events_ts.rds"))
out_events_ts_bT <- bT_ts_data$events_ts

## Part 4: Validate against single point ---------------------------------------

# Should exactly match your single point results
cat("Spatial run — freq at test point:", round(out_freq_bT[ilon, ilat], 2), "\n")
cat("Single point run — freq:          0.19\n")
cat("N events at test point:", out_n_events_bT[ilon, ilat], "\n")

# Check if results actually exist
cat("Non-NA cells in out_freq:", sum(!is.na(out_freq_bT)), "\n")
cat("Range of freq values:", range(out_freq_bT, na.rm = TRUE), "\n")

# Check test point specifically
cat("Freq at test point:", out_freq_bT[ilon, ilat], "\n")  

# Check lon conversion
cat("lon2d_180 range:", range(lon2d, na.rm = TRUE), "\n")
# -180 to 180

# PART 2c: Domain-averaged bT with MHW periods ------------------------------------------

# Load full bT array [ih x jh x time]
nc_bT  <- nc_open(bT_file)
bT_all <- ncvar_get(nc_bT, "tob") 
nc_close(nc_bT)

# Mask land
for (t in 1:n_times) {
  bT_all[, , t][!bottom_mask] <- NA
}

# Domain mean SST at each timestep
bT_mean_ts <- apply(bT_all, 3, mean, na.rm = TRUE)


# Fraction of ocean cells in MHW at each timestep
# out_events_ts is [ih x jh x time] from spatial run
mhw_fraction_bT <- apply(out_events_ts_bT, 3, function(x) {
  mean(x[bottom_mask] == 1, na.rm = TRUE)
})

cat("BMHW fraction range:", round(range(mhw_fraction_bT, na.rm = TRUE), 3), "\n")
# 0.005 - 0.341

## Build DF -------------------------------------------------

df_spatial_bT <- data.frame(
  time         = time_vec,
  bT_mean     = bT_mean_ts,
  mhw_fraction_bT = mhw_fraction_bT,
  month        = as.integer(format(time_vec, "%m"))
)

# Domain-mean P90 threshold
thresh_90_mean <- compute_threshold(df_spatial_bT$bT_mean,
                                    df_spatial_bT$month,
                                    percentile = 0.90,
                                    window = 1)
df_spatial_bT$thresh_90 <- thresh_90_mean[df_spatial_bT$month]

# Flag periods where >= 15% of ocean cells are in MHW
# (threshold can be adjusted needed)
df_spatial_bT$widespread_mhw_bT <- ifelse(df_spatial_bT$mhw_fraction_bT >= 0.15, 1, 0)

## Plot 1: Domain mean SST w widespread MHW shading ----------------------------

ggplot(df_spatial_bT, aes(x = time)) +
  geom_rect(
    data = df_spatial_bT %>% filter(widespread_mhw_bT == 1),
    aes(xmin = time - 15, xmax = time + 15,
        ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.5
  ) +
  geom_line(aes(y = bT_mean),  color = "black", linewidth = 0.75) +
  geom_line(aes(y = thresh_90), color = "black",
            linetype = "dashed", linewidth = 0.5) +
  labs(
    title    = "Domain-Mean bT MHW | MOM6 Hindcast 1993-2025",
    subtitle = paste("Dashed = P90 | Red = >= 15% of ocean cells in a MHW\n",
                     "CHANGES: monthly res; 3-month window; 2-month min duration"),
    y = "SST (°C)", x = NULL
  ) 

## Plot 2: Fraction of cells in MHW over time ----------------------------------
# Shows how spatially extensive each event was

ggplot(df_spatial_bT, aes(x = time, y = mhw_fraction_bT * 100)) +
  geom_area(fill = "red", alpha = 0.4) +
  geom_line(color = "red", linewidth = 0.75) +
  geom_hline(yintercept = 10, linetype = "dashed",
             color = "black", linewidth = 0.5) +
  labs(
    title    = "Spatial Extent of bT MHWs | MOM6 Hindcast 1993-2025",
    subtitle = "% of ocean cells experiencing BMHW conditions each month",
    y = "% of domain in BMHW", x = NULL
  )

## Plot 3: Both together -------------------------------------------------------

p1 <- ggplot(df_spatial_bT, aes(x = time)) +
  geom_rect(
    data = df_spatial_bT %>% filter(widespread_mhw_bT == 1),
    aes(xmin = time - 15, xmax = time + 15,
        ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.4
  ) +
  geom_line(aes(y = bT_mean_ts),  color = "black", linewidth = 0.75) +
  geom_line(aes(y = thresh_90), color = "black",
            linetype = "dashed", linewidth = 0.5) +
  labs(title    = "Domain-Mean bT | Red = widespread BMHW (>= 15% cells)",
       y = "bT (°C)", x = NULL) 

p2 <- ggplot(df_spatial_bT, aes(x = time, y = mhw_fraction_bT * 100)) +
  geom_area(fill = "red", alpha = 0.4) +
  geom_line(color = "darkred", linewidth = 0.75) +
  geom_hline(yintercept = 10, linetype = "dashed", linewidth = 0.5) +
  labs(y = "% domain in BMHW", x = NULL) 

p1 / p2 

# ══════════════════════════════════════════════════════════════════════════════
# PART 3:  COMPOUND SMHW + BMHW ---------------------------------------------
# Freeman Section 2.3: co-occurrence >= 5 days
# Freeman Section 2.4: intensity = Z_SST * Z_bT (standardized)
# CHANGE: co-ocurrence >= 2 months
# ══════════════════════════════════════════════════════════════════════════════

df_compare <- data.frame(
  time = df$time,
  sst = df$sst,
  bT = df_bT$bT,
  events_SST = df$events,
  events_bT  = df_bT$events_bT
)

## Step 1: Flag raw co-occurrence ----------------------------------------------
compound_raw <- ifelse(
  df_compare$events_SST == 1 & df_compare$events_bT == 1, 1, 0
)

## Step 2: Apply minimum duration to compound events (Freeman 2.3) ------
df_compare$compound <- remove_short_events(compound_raw, min_duration = 2)
cat("Compound SMHW + BMHW months:", sum(df_compare$compound), "\n")
# compound months: 2
 
## Step 3: Compound intensity = Z_SST * Z_bT (Freeman 2.4) -------------
z_SST <- compute_anomaly(df_compare$sst, as.integer(format(df_compare$time, "%m")))
z_bT  <- compute_anomaly(df_compare$bT,  as.integer(format(df_compare$time, "%m")))

raw_intensity <- z_SST * z_bT
df_compare$compound_intensity <- (raw_intensity - mean(raw_intensity, na.rm = TRUE)) /
  sd(raw_intensity, na.rm = TRUE)

## Step 4: Characterize compound events ----------------------------------------
compound_events <- characterize_events(
  df_compare$compound,
  df_compare$compound_intensity,
  df_compare$time
)
# event frequency: 0.03 events/yr
print(compound_events)
#        start        end duration_months mean_intensity
# 1 2015-08-16 2015-09-16               2          2.359

## Step 5: Plot ----------------------------------------------------------------
ggplot(df_compare, aes(x = time)) +
  geom_rect(
    data = df_compare %>% filter(compound == 1),
    aes(xmin = time - 15, xmax = time + 15, ymin = -Inf, ymax = Inf),
    fill = NA, color = "gold", linewidth = 0.9
  ) +
  geom_line(aes(y = scale(z_SST)[, 1], color = "SST"),         linewidth = 0.75) +
  geom_line(aes(y = scale(z_bT)[, 1],  color = "Bottom Temp"), linewidth = 0.75) +
  scale_color_manual(values = c("SST" = "#D55E00", "Bottom Temp" = "#0072B2")) +
  labs(
    title    = paste("Compound SMHW + BMHW | lon: -122.02",
                     "lat:", round(lat2d[ilon, ilat], 2)),
    subtitle = paste(
      "Freeman et al. method | Gold = compound\n",
      "Compound intensity = Z_SST x Z_bT (standardized)"
    ),
    y = "Standardized Anomaly (z-score)", x = NULL, color = NULL
  ) + theme_light()


print(sst_events)
print(bT_events)
print(compound_events)

# PART 3b: MHW Spatial Summaries -----------------------------------------------

# SST anomaly grid
nc_sst <- nc_open(sst_file)
sst_all <- ncvar_get(nc_sst, "tos")
nc_close(nc_sst)
sst_all[!array(ocean_mask, dim = dim(sst_all))] <- NA

anomaly_sst_grid <- array(NA_real_, dim = dim(sst_all))
for (i in 1:n_ih) for (j in 1:n_jh) {
  if (!ocean_mask[i, j]) next
  anomaly_sst_grid[i, j, ] <- compute_anomaly(sst_all[i, j, ], months)
}

# bT anomaly grid
nc_bT <- nc_open(bT_file)
bT_all <- ncvar_get(nc_bT, "tob")
nc_close(nc_bT)
bT_all[!array(bottom_mask, dim = dim(bT_all))] <- NA

anomaly_bT_grid <- array(NA_real_, dim = dim(bT_all))
for (i in 1:n_ih) for (j in 1:n_jh) {
  if (!bottom_mask[i, j]) next
  anomaly_bT_grid[i, j, ] <- compute_anomaly(bT_all[i, j, ], months)
}

saveRDS(anomaly_sst_grid, file.path(output_dir, "anomaly_sst_grid.rds"))
saveRDS(anomaly_bT_grid,  file.path(output_dir, "anomaly_bT_grid.rds"))

fraction_over_time <- function(events_array, mask) {
  apply(events_array, 3, function(x) {
    mean(x[mask] == 1, na.rm = TRUE) 
  })
}

## Step 1: compute fractions ---------------------------------------------------

# SST MHW fraction
frac_sst <- fraction_over_time(out_events_ts, ocean_mask)

# BHMW
frac_bT <- fraction_over_time(out_events_ts_bT, bottom_mask)

# SHX
# frac_SHX <- fraction_over_time(out_events_ts_shx, ocean_mask)

## Step 2: build combined df ---------------------------------------------------

df_fractions <- data.frame(
  time = time_vec,
  sst = frac_sst * 100,
  bT = frac_bT * 100
  # shx = frac_shx * 100
)

## Part 3: Plots  --------------------------------------------------------------

### Plot 1 - SST MHW Fraction --------------------------------------------------
p_sst <- ggplot(df_fractions, aes(x = time, y = sst)) +
  geom_area(fill = "red", alpha = 0.4) +
  geom_line(color = "darkred", linewidth = 0.75) +
  labs(
    title    = "Spatial Extent of Surface MHWs | MOM6 1993-2025",
    y        = "% domain in SMHW",
    x        = NULL
  )

### Plot 2 - bT MHW Fraction ---------------------------------------------------

p_bT <- ggplot(df_fractions, aes(x = time, y = bT)) +
  geom_area(fill = "blue", alpha = 0.4) +
  geom_line(color = "darkblue", linewidth = 0.75) +
  labs(
    title    = "Spatial Extent of Bottom MHWs | MOM6 1993-2025",
    subtitle = "depth < 1000m only",
    y        = "% domain in BMHW",
    x        = NULL
  )


### Plot 3 - Compound Fraction -------------------------------------------------

# Co-occurrence: cell must be in both SST and bottom MHW simultaneously
compound_spatial <- out_events_ts * out_events_ts_bT

# Use bottom_mask so only cells valid for BOTH are included, bottom is stricter than ocean
frac_compound <- fraction_over_time(compound_spatial, bottom_mask) * 100

df_fractions$compound <- frac_compound

p_compound <- ggplot(df_fractions, aes(x = time, y = compound)) +
  geom_area(fill = "purple", alpha = 0.4) +
  geom_line(color = "purple4", linewidth = 0.75) +
  labs(
    title    = "Spatial Extent of Compound SMHW + BMHW | MOM6 1993-2025",
    y        = "% domain in compound",
    x        = NULL
  )

### Plot 4 - All 3 stacked -----------------------------------------------------

p_sst / p_bT / p_compound +
  plot_annotation(
    title    = "Marine Heatwave Spatial Extent | MOM6 Hindcast 1993-2025",
    subtitle = "Monthly resolution"
  )

# ══════════════════════════════════════════════════════════════════════════════
# PART 4: SHX DETECTION (SHALLOW HYPOXIA EXTREME) ------------------------------
# Following Freeman et al. (2026) Section 2.3
# ADAPTATIONS FROM FREEMAN:
#   - Monthly resolution (Freeman: daily)
#   - 3-month window for threshold (Freeman: 11-day)
#   - 2-month minimum duration (Freeman: 5 days)
#   - Merge gap 1 month (Freeman: 3 days)
# ══════════════════════════════════════════════════════════════════════════════

## Step 1: Explore 3D oxygen file ----------------------------------------------

nc_o2 <- nc_open(o2_3d_file)
names(nc_o2$var)# oxygen var names
# "average_DT" "average_T1" "average_T2" "o2"         "time_bnds"  "volcello" 
names(nc_o2$dim) # dim names
# "time" "ih"   "jh"   "nv"   "z_l"  "z_i"

# check depth dim
nc_o2$dim$z_l$len # 52
nc_o2$dim$z_l$vals 

nc_o2$var$o2$dim     # [ih, jh, z_l, time] 

sapply(nc_o2$var$o2$dim, \(x) x$name)

nc_close(nc_o2)

# check dim orderfor o2
# ncvar_get start/count has to match this order
sapply(nc_o2$var$o2$dim, function(d) d$name)
# output: ih, jh, z_l, time

# depth level check
z_levels <- nc_o2$dim$z_l$vals
n_depths  <- length(z_levels)
cat("Number of depth levels:", n_depths, "\n")
# 52
cat("Depth range:", range(z_levels), "m\n")
# 2.5 - 6250 m
cat("First few depths:", head(z_levels), "\n")  
# 2.5, 7.5, 12.5, 12.5, 17.5, 22.5, 27.5

# Check n timesteps
cat("N timesteps:", nc_o2$dim$time$len, "\n")  #  timesteps: 390

nc_grid <- nc_open(grid_file)

lon_grid <- ncvar_get(nc_grid, "geolon")
lat_grid <- ncvar_get(nc_grid, "geolat")

nc_close(nc_grid)

# Extract one profile at test point, first timestep to verify structure
extract_o2 <- function(file, ilon, ilat) {
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  
  ncvar_get(nc, "o2",
            start = c(ilon, ilat, 1, 1),
            count = c(1, 1, -1, 1))
}
o2_test <- extract_o2(o2_3d_file, ilon, ilat)
cat("Profile length:", length(o2_test), "\n")   # should = n_depths
cat("O2 profile range:", range(o2_test, na.rm = TRUE), "\n")  # check units

o2_test_umol <- o2_test * 1e6
range(o2_test_umol, na.rm = TRUE)
# 59.1699 258.4101 µmol/kg

nc_o2 <- nc_open(o2_3d_file)

nc_o2$var$o2$units

nc_close(nc_o2)

## Step 2: Load depth levels ---------------------------------------------------

nc_o2 <- nc_open(o2_3d_file)
z_levels <- nc_o2$dim$z_l$vals   # depth of each layer in meters (positive = down)
n_depths <- length(z_levels)
cat("Number of depth levels:", n_depths, "\n") # 52
cat("Depth range:", range(z_levels), "m\n") # 2.5 - 6250
nc_close(nc_o2)


## Step 3: Extract full O2 profile time series at test point ----------------------------

# Dimensions: [ih, jh, z_l, time] = [342, 816, n_depths, 390]
# Extract all depths at test point across all timesteps
# Result: matrix of [n_depths x n_times]

nc_o2 <- nc_open(o2_3d_file)
o2_profile_ts <- ncvar_get(nc_o2, "o2",
                           start = c(ilon, ilat, 1, 1),
                           count = c(1, 1, -1, -1))
# o2_profile_ts is now [n_depths x n_times]
nc_close(nc_o2)

cat("O2 profile dimensions:", dim(o2_profile_ts), "\n")
# [52, 390]
cat("O2 range:", range(o2_profile_ts, na.rm = TRUE))
# Converted from MOM6 mol/kg using rho=1025 kg/m3

## Step 4: Unit Conversion -----------------------------------------------------

# Freeman threshold: 2 mg/L = 62/5 mmol/m3
# MOM6 oxygen units: mol/kg
# Convert to mmol/m3:
# mol/kg × seawater density (kg/m3) × 1000 mmol/mol

rho_sw <- 1025  # kg/m3

o2_profile_ts <- o2_profile_ts * rho_sw * 1000

cat("Converted O2 range:",
    range(o2_profile_ts, na.rm = TRUE),
    "mmol/m3\n")

# Freeman threshold:
# 2 mg/L = 62.5 mmol/m3
O2_threshold <- 62.5


## Step 5: Compute HLD at each timestep ----------------------------------------

# Freeman: linearly interpolate vert O2 profile to find shallowest 
# depth where O2 crosses 62.5 mmol/m3 (hypoxic threshold)

# Why linear interpolation:O2 doesnt jump between model layers
# interpolating gives more precise depth than just taking the nearest layer

# Vectorized HLD computation for an entire [depth x time] matrix
# Assumes z_levels is already sorted ascending (surface to bottom)
compute_HLD_matrix <- function(o2_mat, z_levels, threshold = 62.5) {
  n_depths <- nrow(o2_mat)
  n_times  <- ncol(o2_mat)
  
  below <- o2_mat < threshold          # logical matrix [depth x time]
  below[is.na(below)] <- FALSE
  
  # first depth index where hypoxic, per column (timestep)
  # apply(..., 2, which.max) on a logical matrix returns first TRUE index
  has_hypoxia <- colSums(below) > 0
  first_idx   <- apply(below, 2, function(x) which(x)[1])
  
  HLD <- rep(NA_real_, n_times)
  
  # surface already hypoxic -> HLD = 0
  surf_hypoxic <- has_hypoxia & (first_idx == 1)
  HLD[surf_hypoxic] <- 0
  
  # interior crossing -> linear interpolation
  interior <- has_hypoxia & (first_idx > 1) & !is.na(first_idx)
  idx2 <- first_idx[interior]
  idx1 <- idx2 - 1
  cols <- which(interior)
  
  z1 <- z_levels[idx1]
  z2 <- z_levels[idx2]
  o1 <- o2_mat[cbind(idx1, cols)]
  o2 <- o2_mat[cbind(idx2, cols)]
  
  HLD[interior] <- z1 + (threshold - o1) * (z2 - z1) / (o2 - o1)
  
  # fully oxygenated columns stay NA (has_hypoxia == FALSE)
  return(HLD)
}


z_levels <- z_levels[order(z_levels)]
HLD_ts <- compute_HLD_matrix(o2_profile_ts, z_levels, threshold = O2_threshold)

cat("HLD range:", range(HLD_ts, na.rm = TRUE), "m\n")
# HLD range: 184.5534 324.3789 m
cat("Timesteps with hypoxia (defined HLD):", sum(!is.na(HLD_ts)), "\n")
# Timesteps w hypoxia (defined HLD ): 385
cat("Timesteps fully oxygenated (undefined HLD):", sum(is.na(HLD_ts)), "\n")
# timesteps fully oxygenated (undefined HLD): 5


## Step 6: Handle Undefined HLDs -----------------------------------------------

# Freeman 2.3: two separate treatments for undefined HLD

# FOR THRESHOLD COMPUTATION (HLD_for_SHX):
#   Replace undefined HLD with bottom depth at that location
#   WHY: represents a fully oxygenated column. HLD effectively at the seafloor
#   Keeps the time series continuous so percentile can be computed properly
#   Without this, NA-heavy series would bias the P10 threshold
#
# FOR EVENT DETECTION:
#   Undefined HLD timesteps are NOT valid extreme days
#   WHY: a fully oxygenated column cannot be a shallow hypoxia extreme
#   Even if HLD_for_SHX falls below P10, it's excluded from event flagging

bottom_depth_here <- depth[ilon, ilat]
cat("Bottom depth at test point:", bottom_depth_here, "m\n")
# 312.7541 m

# HLD for threshold: undefined -> bottom depth
HLD_for_SHX <- ifelse(is.na(HLD_ts), bottom_depth_here, HLD_ts)

# true only where hypoxia was acc present
valid_extreme <- !is.na(HLD_ts)

cat("HLD_for_SHX range:", range(HLD_for_SHX, na.rm = TRUE), "m\n")
# HLD_for_SHX range: 184.5534 324.3789 m

# Quick plot of raw HLD time series
plot(time_vec, HLD_ts, type = "l",
     main = "Hypoxic Layer Depth (HLD) at test point",
     ylab = "Depth (m)", xlab = NULL,
     ylim = c(0, bottom_depth_here))
abline(h = bottom_depth_here, lty = 2, col = "gray")


## Step 7: Build O2 DF ----------------------------------------------------------

df_o2 <- data.frame(
  time          = time_vec,
  HLD           = HLD_ts, # raw HLD (NA where fully oxygenated)
  HLD_for_SHX   = HLD_for_SHX, # NA replaced with bottom depth
  valid_extreme  = valid_extreme, # TRUE = hypoxia present, valid for event ID
  month          = as.integer(format(time_vec, "%m"))
)

## Step 8: Compute P10 Threshold for SHX ---------------------------------------

# Freeman: SHX = HLD shallower than P10 (lower depth value = shallower = worse)
# WHY P10 not P90: shallower HLD = more of the water column is hypoxic
# A smaller depth number means hypoxia reaches higher up. = ecological stress
# ADAPTATION: 3-month window (Freeman: 11-day)
# Computed on HLD_for_SHX (with NAs replaced) not raw HLD

thresh_10_HLD <- compute_threshold(df_o2$HLD_for_SHX,
                                   df_o2$month,
                                   percentile = 0.10,
                                   window = 1)
df_o2$thresh_10 <- thresh_10_HLD[df_o2$month]

print(data.frame(month = 1:12, HLD_P10 = round(thresh_10_HLD, 1)))
#    month HLD_P10
# 1      1   238.0
# 2      2   236.1
# 3      3   221.3
# 4      4   214.5
# 5      5   223.5
# 6      6   230.7
# 7      7   244.1
# 8      8   249.4
# 9      9   243.8
# 10    10   236.5
# 11    11   234.7
# 12    12   234.0
# seasonal variation in HLD threshold

## Step 9:Flag SHX exceedances  ------------------------------------------------

# SHX: HLD_for_SHX <= P10 AND it was a valid hypoxic timestep
# WHY both conditions: threshold computed on continuous series including NAs
# replaced with bottom depth, but those replaced timesteps cannot be extremes
# even if they fall below P10 (a fully oxygenated column is never an SHX)

df_o2$exceeds_SHX <- ifelse(
  df_o2$HLD_for_SHX <= df_o2$thresh_10 & df_o2$valid_extreme,
  1, 0
)

cat("Months flagged as SHX before filtering:", sum(df_o2$exceeds_SHX), "\n")
# 36

## Step 10: Merge close events and remove short events -------------------------

# same as SST/bT MHW
# ADAPTATION: merge gap 1 month, min duration 2 months
df_o2$exceeds_merged <- merge_events(df_o2$exceeds_SHX, max_gap = 1)
df_o2$events_SHX <- remove_short_events(df_o2$exceeds_merged, min_duration = 2)

cat("Months flagged as SHX after filtering", sum(df_o2$events_SHX), "\n")
# 34 after filtering. 2 were removed

## Step 11: Anomaly and characterization ---------------------------------------

# Anomaly based on HLD_for_SHX z-score
# WHY HLD_for_SHX and not raw HLD: keeps series continuous for anomaly computation

df_o2$anomaly_HLD <- compute_anomaly(df_o2$HLD_for_SHX, df_o2$month)

shx_events <- characterize_events(df_o2$events_SHX,
                                  df_o2$anomaly_HLD,
                                  df_o2$time)
# event freq: 0.19/yr

print(shx_events)
#        start        end duration_months mean_intensity
# 1 2002-04-16 2002-05-16               2         -1.463
# 2 2010-07-16 2011-11-16              17         -1.721
# 3 2012-02-15 2012-03-16               2         -1.881
# 4 2017-10-16 2017-11-16               2         -1.057
# 5 2020-11-16 2021-04-16               6         -1.721
# 6 2021-08-16 2021-12-16               5         -1.403

## Step 12: Plot SHX -----------------------------------------------------------

ggplot(df_o2, aes(x = time)) +
  geom_rect(
    data = df_o2 %>% filter(events_SHX == 1),
    aes(xmin = time - 15, xmax = time + 15,
        ymin = -Inf, ymax = Inf),
    fill = "steelblue", alpha = 0.5
  ) +
  geom_line(aes(y = HLD_for_SHX), color = "black", linewidth = 0.75) +
  geom_line(aes(y = thresh_10),   color = "black",
            linetype = "dashed",  linewidth = 0.5) +
  # Mark timesteps where HLD was undefined (fully oxygenated)
  geom_point(data = df_o2 %>% filter(!valid_extreme),
             aes(y = bottom_depth_here),
             color = "gold", size = 3, alpha = 0.7) +
  scale_y_reverse() +   # depth increases downward
  labs(
    title    = paste("SHX Detection | lon:",
                     round(lon2d[ilon, ilat] - 360, 2),
                     "lat:", round(lat2d[ilon, ilat], 2)),
    subtitle = paste(
      "P10 threshold (dashed) | Blue = SHX event\n",
      "Gold dots = fully oxygenated (undefined HLD, excluded from events)\n"
    ),
    y = "HLD Depth (m, positive down)", x = NULL
  )

## Step 13: Compound MHW + SHX -------------------------------------------------

# Freeman 2.3: compound = co-occurrence >= 5 days
# Freeman Section 2.4: intensity = Z_MHW * Z_SHX, standardized
# ADAPTATION: co-occurrence >= 2 months

df_compound_shx <- df %>%
  select(time, events, anomaly) %>%
  rename(events_MHW = events, z_MHW = anomaly) %>%
  inner_join(
    df_o2 %>% select(time, events_SHX, anomaly_HLD) %>%
      rename(z_SHX = anomaly_HLD),
    by = "time"
  )

# Flag co-occurrence
compound_raw_shx <- ifelse(
  df_compound_shx$events_MHW == 1 & df_compound_shx$events_SHX == 1, 1, 0
)

# Apply minimum duration filter
df_compound_shx$compound <- remove_short_events(compound_raw_shx, min_duration = 2)

cat("Compound MHW + SHX months:", sum(df_compound_shx$compound), "\n")
# ZERO compound months

# Compound intensity = Z_MHW * Z_SHX, standardized (Freeman Section 2.4)
raw_int <- df_compound_shx$z_MHW * df_compound_shx$z_SHX
df_compound_shx$compound_intensity <- (raw_int - mean(raw_int, na.rm = TRUE)) /
  sd(raw_int, na.rm = TRUE)

compound_shx_events <- characterize_events(
  df_compound_shx$compound,
  df_compound_shx$compound_intensity,
  df_compound_shx$time
)
print(compound_shx_events)
# nothing 

# PART 4b: SHX Full Grid -------------------------------------------------------

shx_mask <- ocean_mask

## Step 1: Output arrays -------------------------------------------------------

out_n_events_shx       <- matrix(NA_real_, n_ih, n_jh)
out_freq_shx           <- matrix(NA_real_, n_ih, n_jh)
out_mean_duration_shx  <- matrix(NA_real_, n_ih, n_jh)
out_mean_intensity_shx <- matrix(NA_real_, n_ih, n_jh)

out_events_ts_shx <- array(NA_integer_, dim = c(n_ih, n_jh, n_times))

## Step 2: Chunk Loop ------------------------------------------------------------

# nc_o2      <- nc_open(o2_3d_file)
# start_time <- Sys.time()
# 
# z_levels <- z_levels[order(z_levels)]
# 
# for (j in 1:n_jh) {
#   
#   if (j %% 50 == 0) {
#     elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
#     cat("Row", j, "of", n_jh, "| Elapsed:", elapsed, "mins\n")
#   }
#   
#   if (!any(shx_mask[, j], na.rm = TRUE)) next
#   
#   o2_row <- ncvar_get(nc_o2, "o2", start = c(1, j, 1, 1), count = c(-1, 1, -1, -1))
#   o2_row <- o2_row * rho_sw * 1000
#   
#   # Step A: compute HLD_for_SHX for every cell in the row in one vectorized call
#   o2_row_perm <- aperm(o2_row, c(2, 1, 3))              # [depth, ih, time]
#   dim(o2_row_perm) <- c(n_depths, n_ih * n_times)       # [depth, ih*time]
#   
#   HLD_flat <- compute_HLD_matrix(o2_row_perm, z_levels, threshold = O2_threshold)
#   HLD_mat_row <- matrix(HLD_flat, nrow = n_ih, ncol = n_times)   # back to [ih, time]
#   
#   bottom_depth_row <- depth[, j]  # length n_ih
#   HLD_for_SHX_mat   <- HLD_mat_row
#   HLD_for_SHX_mat[is.na(HLD_mat_row)] <- bottom_depth_row[row(HLD_mat_row)[is.na(HLD_mat_row)]]
#   valid_extreme_mat <- !is.na(HLD_mat_row)
#  
#    # Mask out land / non-shx cells so they don't pollute the threshold calc below
#   land_rows <- !shx_mask[, j]
#   HLD_for_SHX_mat[land_rows, ]   <- NA
#   valid_extreme_mat[land_rows, ] <- FALSE
#   
#   
#   # Step B: vectorized threshold across the row
#   thresh_mat_row <- compute_threshold_matrix(HLD_for_SHX_mat, months, percentile = 0.10, window = 1)
#   thresh_ts_mat  <- thresh_mat_row[, months]
#   
#   exceeds_mat <- ifelse(HLD_for_SHX_mat <= thresh_ts_mat & valid_extreme_mat, 1, 0)
#   exceeds_mat[is.na(exceeds_mat)] <- 0
#   
#   # Step C: per-cell event characterization (cheap)
#   for (i in 1:n_ih) {
#     
#     if (!shx_mask[i, j]) next
#     if (all(is.na(HLD_for_SHX_mat[i, ]))) next
#     
#     exceeds <- exceeds_mat[i, ]
#     
#     exceeds_merged <- merge_events(exceeds, max_gap = 1)
#     events         <- remove_short_events(exceeds_merged, min_duration = 2)
#     
#     anomaly <- compute_anomaly(HLD_for_SHX_mat[i, ], months)
#     
#     r      <- rle(events)
#     ends   <- cumsum(r$lengths)
#     starts <- ends - r$lengths + 1
#     
#     durations   <- c()
#     intensities <- c()
#     for (k in seq_along(r$values)) {
#       if (r$values[k] == 1) {
#         idx         <- starts[k]:ends[k]
#         durations   <- c(durations,   length(idx))
#         intensities <- c(intensities, mean(anomaly[idx], na.rm = TRUE))
#       }
#     }
#     
#     out_n_events_shx[i, j]       <- length(durations)
#     out_freq_shx[i, j]           <- length(durations) / n_years
#     out_mean_duration_shx[i, j]  <- ifelse(length(durations) > 0, mean(durations), 0)
#     out_mean_intensity_shx[i, j] <- ifelse(length(intensities) > 0, mean(intensities), 0)
#     out_events_ts_shx[i, j, ]    <- events
#   }
# }
# 
# nc_close(nc_o2)
# cat("Done! Total time:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "mins\n")


## Step 3: Save -----------------------------------------------------------------

# output_dir <- 'C:/Users/mcclo/OneDrive/Documents/Hollings_Internship/Outputs'
# 
# saveRDS(list(
#   n_events       = out_n_events_shx,
#   freq           = out_freq_shx,
#   mean_duration  = out_mean_duration_shx,
#   mean_intensity = out_mean_intensity_shx,
#   lon2d          = lon2d,
#   lat2d          = lat2d,
#   shx_mask       = shx_mask
# ), file.path(output_dir, "shx_summary.rds"))
# 
# saveRDS(list(
#   events_ts = out_events_ts_shx,
#   time_vec  = time_vec,
#   lon2d     = lon2d,
#   lat2d     = lat2d
# ), file.path(output_dir, "shx_events_ts.rds"))

file.exists(file.path(output_dir, "shx_summary.rds"))
file.exists(file.path(output_dir, "shx_events_ts.rds"))
file.info(file.path(output_dir, "shx_summary.rds"))$size    
file.info(file.path(output_dir, "shx_events_ts.rds"))$size  

# saveRDS(list(
#   df_o2               = df_o2,
#   shx_events          = shx_events,
#   compound_shx_events = compound_shx_events,
#   df_compound_shx     = df_compound_shx,
#   o2_profile_ts       = o2_profile_ts,
#   HLD_ts              = HLD_ts,
#   depth               = depth,
#   z_levels            = z_levels,
#   O2_threshold        = O2_threshold,
#   rho_sw              = rho_sw,
#   n_years             = n_years,
#   months              = months,
#   ilon                = ilon,
#   ilat                = ilat
# ), file.path(output_dir, "shx_point_and_params.rds"))
# 
# file.exists(file.path(output_dir, "shx_point_and_params.rds"))
# file.info(file.path(output_dir, "shx_point_and_params.rds"))$size

shx_point <- readRDS(file.path(output_dir, "shx_point_and_params.rds"))
df_o2               <- shx_point$df_o2
shx_events          <- shx_point$shx_events
compound_shx_events <- shx_point$compound_shx_events
df_compound_shx     <- shx_point$df_compound_shx
o2_profile_ts       <- shx_point$o2_profile_ts
HLD_ts              <- shx_point$HLD_ts
depth               <- shx_point$depth
z_levels            <- shx_point$z_levels
O2_threshold        <- shx_point$O2_threshold
rho_sw              <- shx_point$rho_sw
n_years             <- shx_point$n_years
months              <- shx_point$months
ilon                <- shx_point$ilon
ilat                <- shx_point$ilat

shx_summary <- readRDS(file.path(output_dir, "shx_summary.rds"))
out_n_events_shx       <- shx_summary$n_events
out_freq_shx           <- shx_summary$freq
out_mean_duration_shx  <- shx_summary$mean_duration
out_mean_intensity_shx <- shx_summary$mean_intensity
shx_mask               <- shx_summary$shx_mask

shx_ts_data       <- readRDS(file.path(output_dir, "shx_events_ts.rds"))
out_events_ts_shx <- shx_ts_data$events_ts

## Step 4: Validate against single point ------------------------------------------
cat("Spatial run — freq at test point:", round(out_freq_shx[ilon, ilat], 2), "\n")
cat("Single point run — freq:          0.19\n")
cat("N events at test point:", out_n_events_shx[ilon, ilat], "(expect 6)\n")


## Step 5: Spatial extent plot -------------------------------------------------


frac_shx <- fraction_over_time(out_events_ts_shx, shx_mask) * 100

df_fractions$shx <- frac_shx

p_shx <- ggplot(df_fractions, aes(x = time, y = shx)) +
  geom_area(fill = "steelblue", alpha = 0.4) +
  geom_line(color = "steelblue4", linewidth = 0.75) +
  geom_hline(yintercept = 5, linetype = "dashed", linewidth = 0.5) +
  labs(
    title    = "Spatial Extent of SHX | MOM6 1993-2025",
    subtitle = "Dashed = 5% threshold",
    y        = "% domain in SHX",
    x        = NULL
  )


# Stack all four
p_sst / p_bT / p_shx / p_compound +
  plot_annotation(
    title    = "Marine Extreme Spatial Extent | MOM6 Northeast Pacific Hindcast 1993-2025"
  )



# PART 5: 2D Mapping  -----------------------------------------------------------
nc_grid  <- nc_open(grid_file)
lon_c <- ncvar_get(nc_grid, "geolon_c")  # corner longitudes
lat_c <- ncvar_get(nc_grid, "geolat_c")  # corner latitudes
nc_close(nc_grid)

cat("Center grid dims:", dim(lon2d), "\n")    #  342 x 816
cat("Corner grid dims:", dim(lon_c), "\n")    #  343 x 817

lon_c_180 <- ifelse(lon_c > 180, lon_c - 360, lon_c)

## Convert to sf and plot -----------------------------------------------------

build_sf_grid <- function(mat, mask, every_nth = 1, max_lon_span = 10) {
  
  mat[!mask] <- NA
  
  i_idx <- seq(1, n_ih, by = every_nth)
  j_idx <- seq(1, n_jh, by = every_nth)
  
  polys  <- vector("list", sum(!is.na(mat[i_idx, j_idx])))
  values <- numeric(length(polys))
  cell_id <- 0
  
  for (j in j_idx) {
    for (i in i_idx) {
      
      if (is.na(mat[i, j])) next
      
      lons <- c(lon_c_180[i,   j], lon_c_180[i+1, j],
                lon_c_180[i+1, j+1], lon_c_180[i, j+1],
                lon_c_180[i,   j])
      lats <- c(lat_c[i,   j], lat_c[i+1, j],
                lat_c[i+1, j+1], lat_c[i, j+1],
                lat_c[i,   j])
      
      if (diff(range(lons)) > max_lon_span) next
      
      cell_id <- cell_id + 1
      polys[[cell_id]]  <- st_polygon(list(cbind(lons, lats)))
      values[cell_id]   <- mat[i, j]
    }
  }
  
  if (cell_id == 0) {
    warning("No valid cells found for this mask. returning empty sf object")
    return(st_sf(value = numeric(0), geometry = st_sfc(crs = 4326)))
  }
  
  polys  <- polys[1:cell_id]
  values <- values[1:cell_id]
  
  st_sf(value = values, geometry = st_sfc(polys, crs = 4326))
}

plot_spatial_sf <- function(mat, mask, title, legend_label, every_nth = 1, max_lon_span = 10) {
  
  sf_grid <- build_sf_grid(mat, mask, every_nth, max_lon_span)
  
  ggplot(sf_grid) +
    geom_sf(aes(fill = value), color = NA) +
    scale_fill_viridis_c(name = legend_label, na.value = "gray90") +
    coord_sf(xlim = c(-180, -100), ylim = c(10, 80)) +
    labs(title = title, subtitle = "MOM6 1993-2025 | Monthly",
         x = "Longitude", y = "Latitude")
}

# every_nth = 1 for full resolution
shx_plot_mask <- ocean_mask & !is.na(out_freq_shx)

# SST
plot_spatial_sf(out_freq, ocean_mask,  "SST MHW Frequency",      "Events/yr",  every_nth = 1)
plot_spatial_sf(out_mean_duration, ocean_mask,  "SST MHW Mean Duration",  "Months",     every_nth = 1)
plot_spatial_sf(out_mean_intensity, ocean_mask,  "SST MHW Mean Intensity", "Z-score",    every_nth = 1)

# bT
plot_spatial_sf(out_freq_bT, bottom_mask, "Bottom MHW Frequency",   "Events/yr",  every_nth = 1)
plot_spatial_sf(out_mean_duration_bT, bottom_mask,  "Bottom MHW Mean Duration",  "Months",     every_nth = 1)
plot_spatial_sf(out_mean_intensity_bT, bottom_mask,  "Bottom MHW Mean Intensity", "Z-score",    every_nth = 1)

# SHX
out_mean_intensity_shx <- -out_mean_intensity_shx
plot_spatial_sf(out_freq_shx, shx_plot_mask,  "SHX Frequency",          "Events/yr",  every_nth = 1)
plot_spatial_sf(out_mean_duration_shx, shx_plot_mask,  "SHX Mean Duration",  "Months",     every_nth = 1)
plot_spatial_sf(out_mean_intensity_shx, shx_plot_mask,  "SHX Mean Intensity", "Z-score",    every_nth = 1)


# How many cells had hypoxia at any point?
cells_with_any_hypoxia <- apply(out_events_ts_shx, c(1,2), function(x) any(x == 1, na.rm = TRUE))
cat("Cells with at least one SHX event:", sum(cells_with_any_hypoxia, na.rm = TRUE), "\n")
cat("Total ocean cells:", sum(ocean_mask), "\n")

# Where are they?
plot_spatial_sf(
  ifelse(cells_with_any_hypoxia, 1, NA),
  ocean_mask, 
  "Cells with any SHX event", 
  "1 = yes"
)


# PART 4c: SHX P10 Depth Threshold Map - lean loop -----------------------------
# Skips flagging/merging/duration-filtering/anomaly/characterization
# entirely. Still has to re-read the oxygen file and rerun HLD
# interpolation, since that intermediate array was never saved — but
# this cuts out the heaviest per-cell work from the original loop.
# Only needs to run once; saves its own output at the end.

shx_mask <- ocean_mask
n_ih <- 342; n_jh <- 816; n_depths <- length(z_levels)

out_thresh_shx <- matrix(NA_real_, n_ih, n_jh)
out_anomaly_shx <- array(NA_real_, dim = c(n_ih, n_jh, n_times))

rm(sst_all, bT_all)
gc()

checkpoint_file <- file.path(output_dir, "shx_anomaly_checkpoint.rds")

if (file.exists(checkpoint_file)) {
  cp <- readRDS(checkpoint_file)
  out_thresh_shx  <- cp$out_thresh_shx
  out_anomaly_shx <- cp$out_anomaly_shx
  start_row <- cp$last_row + 1
  cat("Resuming from row", start_row, "\n")
} else {
  out_thresh_shx  <- matrix(NA_real_, n_ih, n_jh)
  out_anomaly_shx <- array(NA_real_, dim = c(n_ih, n_jh, n_times))
  start_row <- 1
}

# nc_o2 <- nc_open(o2_3d_file)
# start_time <- Sys.time()
# 
# for (j in start_row:n_jh) {
#   
#   if (j %% 20 == 0) {
#     elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
#     cat("Row", j, "of", n_jh, "| Elapsed:", elapsed, "mins\n")
#     
#     # checkpoint every 20 rows
#     saveRDS(list(out_thresh_shx = out_thresh_shx,
#                  out_anomaly_shx = out_anomaly_shx,
#                  last_row = j),
#             checkpoint_file)
#   }
#   
#   if (!any(shx_mask[, j], na.rm = TRUE)) next
#   
#   # retry up to 3 times on a transient I/O failure before giving up
#   o2_row <- NULL
#   for (attempt in 1:3) {
#     o2_row <- tryCatch(
#       ncvar_get(nc_o2, "o2", start = c(1, j, 1, 1), count = c(-1, 1, -1, -1)),
#       error = function(e) {
#         cat("Read failed on row", j, "attempt", attempt, ":", conditionMessage(e), "\n")
#         nc_close(nc_o2)
#         Sys.sleep(5)
#         nc_o2 <<- nc_open(o2_3d_file)   # reopen the connection
#         NULL
#       }
#     )
#     if (!is.null(o2_row)) break
#   }
#   if (is.null(o2_row)) {
#     cat("Row", j, "failed after 3 attempts — skipping, will need manual rerun\n")
#     next
#   }
#   
#   o2_row <- o2_row * rho_sw * 1000
#   o2_row_perm <- aperm(o2_row, c(2, 1, 3))
#   dim(o2_row_perm) <- c(n_depths, n_ih * n_times)
#   HLD_flat <- compute_HLD_matrix(o2_row_perm, z_levels, threshold = O2_threshold)
#   HLD_mat_row <- matrix(HLD_flat, nrow = n_ih, ncol = n_times)
#   
#   bottom_depth_row <- depth[, j]
#   HLD_for_SHX_mat <- HLD_mat_row
#   HLD_for_SHX_mat[is.na(HLD_mat_row)] <- bottom_depth_row[row(HLD_mat_row)[is.na(HLD_mat_row)]]
#   
#   land_rows <- !shx_mask[, j]
#   HLD_for_SHX_mat[land_rows, ] <- NA
#   
#   for (i in 1:n_ih) {
#     if (!shx_mask[i, j]) next
#     out_anomaly_shx[i, j, ] <- compute_anomaly(HLD_for_SHX_mat[i, ], months)
#   }
#   
#   thresh_mat_row <- compute_threshold_matrix(HLD_for_SHX_mat, months, percentile = 0.10, window = 1)
#   out_thresh_shx[, j] <- rowMeans(thresh_mat_row, na.rm = TRUE)
# }

# nc_close(nc_o2)
# cat("Done:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "mins\n")
# 
# saveRDS(out_thresh_shx,  file.path(output_dir, "shx_thresh_depth_map.rds"))
# saveRDS(out_anomaly_shx, file.path(output_dir, "anomaly_shx_grid.rds"))
# file.remove(checkpoint_file)
# 
# out_thresh_shx <- readRDS(file.path(output_dir, "shx_thresh_depth_map.rds"))
# out_anomalygrid_shx <- readRDS(file.path(output_dir, "anomaly_shx_grid.rds"))



# ══════════════════════════════════════════════════════════════════════════════
# PART 6: Full-Grid Compound Extremes (SMHW / BMHW / SHX) ------------------------------
# Freeman Section 2.3-2.4 adapted — LCX omitted
# ══════════════════════════════════════════════════════════════════════════════

# Reload if starting fresh
# out_events_ts<- readRDS(file.path(output_dir, "sst_mhw_events_ts.rds"))$events_ts
# out_events_ts_bT  <- readRDS(file.path(output_dir, "bT_mhw_events_ts.rds"))$events_ts
# out_events_ts_shx <- readRDS(file.path(output_dir, "shx_events_ts.rds"))$events_ts

# per-cell frequency/duration/n_events for ANY compound
# defined as element-wise product of 2 or 3 boolean event arrays
compound_grid_full <- function(events_a, events_b, anomaly_a, anomaly_b,
                               valid_mask, n_years, min_duration = 2) {
  n_ih <- dim(events_a)[1]; n_jh <- dim(events_a)[2]
  out_n <- matrix(NA_real_, n_ih, n_jh)
  out_freq <- matrix(NA_real_, n_ih, n_jh)
  out_dur  <- matrix(NA_real_, n_ih, n_jh)
  out_int  <- matrix(NA_real_, n_ih, n_jh)
  
  for (j in 1:n_jh) {
    if (!any(valid_mask[, j], na.rm = TRUE)) next
    for (i in 1:n_ih) {
      if (!valid_mask[i, j]) next
      
      exceeds <- ifelse(events_a[i, j, ] == 1 & events_b[i, j, ] == 1, 1, 0)
      exceeds[is.na(exceeds)] <- 0
      events <- remove_short_events(exceeds, min_duration)
      if (sum(events) == 0) next
      
      # Freeman 2.4: intensity = Z_a x Z_b, standardized
      raw_intensity <- anomaly_a[i, j, ] * anomaly_b[i, j, ]
      z_intensity <- (raw_intensity - mean(raw_intensity, na.rm = TRUE)) /
        sd(raw_intensity, na.rm = TRUE)
      
      r <- rle(events)
      ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
      durations <- c(); intensities <- c()
      for (k in seq_along(r$values)) {
        if (r$values[k] == 1) {
          idx <- starts[k]:ends[k]
          durations <- c(durations, length(idx))
          intensities <- c(intensities, mean(z_intensity[idx], na.rm = TRUE))
        }
      }
      out_n[i, j] <- length(durations)
      out_freq[i, j] <- length(durations) / n_years
      out_dur[i, j] <- mean(durations)
      out_int[i, j] <- mean(intensities)
    }
  }
  list(n_events = out_n, freq = out_freq, mean_duration = out_dur, mean_intensity = out_int)
}

res_sst_bT  <- compound_grid_full(out_events_ts, out_events_ts_bT,  anomaly_sst_grid, anomaly_bT_grid,  bottom_mask, n_years)
res_sst_shx <- compound_grid_full(out_events_ts, out_events_ts_shx, anomaly_sst_grid, out_anomaly_shx,  ocean_mask,  n_years)
res_bT_shx  <- compound_grid_full(out_events_ts_bT, out_events_ts_shx, anomaly_bT_grid, out_anomaly_shx, bottom_mask, n_years)

## SMHW × BMHW (full grid)
compound_sst_bT <- out_events_ts * out_events_ts_bT
res_sst_bT <- compound_grid_freq_duration(compound_sst_bT, bottom_mask, n_years)

## SMHW × SHX (full grid)
compound_sst_shx <- out_events_ts * out_events_ts_shx
res_sst_shx <- compound_grid_freq_duration(compound_sst_shx, ocean_mask, n_years)

## BMHW × SHX (full grid)
compound_bT_shx <- out_events_ts_bT * out_events_ts_shx
res_bT_shx <- compound_grid_freq_duration(compound_bT_shx, bottom_mask, n_years)

## Triple: SMHW × BMHW × SHX (full grid)
compound_triple <- out_events_ts * out_events_ts_bT * out_events_ts_shx
res_triple <- compound_grid_freq_duration(compound_triple, bottom_mask, n_years)

## Save  -----------------------------------------------------------------
saveRDS(list(
   sst_bT  = res_sst_bT,
   sst_shx = res_sst_shx,
   bT_shx  = res_bT_shx,
   triple  = res_triple,
   lon2d = lon2d, lat2d = lat2d
), file.path(output_dir, "compound_grid_summary.rds"))

compound_summary <- readRDS(file.path(output_dir, "compound_grid_summary.rds"))
res_sst_bT  <- compound_summary$sst_bT
res_sst_shx <- compound_summary$sst_shx
res_bT_shx  <- compound_summary$bT_shx
res_triple  <- compound_summary$triple


## Likelihood Multiplication Factor (grid version) -----------------------------

# CHANGE: freeman works from filtered marginal probability of an extreme on a given day
# here I use straight events/year - approx of the idea of the metric instead of reproduction
# Future: compute joint/marginal probabilities from daily-boolean events_ts arrays(fraction of timesteps flagged)

# Freeman 2.8: LMF = observed joint freq / product of marginal freq
compute_LMF <- function(joint_freq, freq_a, freq_b, n_years) {
  # Convert frequencies (events/yr) back to approximate probabilities
  # by comparing joint occurrence rate to expected-by-chance rate
  expected <- (freq_a / n_years) * (freq_b / n_years)  
  ifelse(expected > 0, joint_freq / n_years / expected, NA)
}

LMF_sst_bT  <- compute_LMF(res_sst_bT$freq, out_freq, out_freq_bT, n_years)
LMF_sst_shx <- compute_LMF(res_sst_shx$freq, out_freq, out_freq_shx, n_years)
LMF_bT_shx  <- compute_LMF(res_bT_shx$freq, out_freq_bT, out_freq_shx, n_years)

saveRDS(list(LMF_sst_bT = LMF_sst_bT, LMF_sst_shx = LMF_sst_shx, LMF_bT_shx = LMF_bT_shx),
        file.path(output_dir, "LMF_grid.rds"))

LMF_data    <- readRDS(file.path(output_dir, "LMF_grid.rds"))
LMF_sst_bT  <- LMF_data$LMF_sst_bT
LMF_sst_shx <- LMF_data$LMF_sst_shx
LMF_bT_shx  <- LMF_data$LMF_bT_shx


## Plot --------------------------------------------------------------------------
frac_sst_bT   <- fraction_over_time(compound_sst_bT,  bottom_mask) * 100
frac_sst_shx  <- fraction_over_time(compound_sst_shx, ocean_mask)  * 100
frac_bT_shx   <- fraction_over_time(compound_bT_shx,  bottom_mask) * 100
frac_triple   <- fraction_over_time(compound_triple,  bottom_mask) * 100

df_fractions$compound_sst_shx <- frac_sst_shx
df_fractions$compound_bT_shx <- frac_bT_shx
df_fractions$triple <- frac_triple

p_triple <- ggplot(df_fractions, aes(x = time, y = triple)) +
  geom_area(fill = "black", alpha = 0.4) +
  geom_line(color = "black", linewidth = 0.75) +
  labs(title = "Spatial Extent of Triple Compound (SMHW+BMHW+SHX) | MOM6 1993-2025",
       y = "% domain", x = NULL)

p_sst / p_bT / p_shx / p_compound / p_triple

## Maps ------------------------------------------------------------------------

plot_spatial_sf(res_sst_shx$freq, ocean_mask,  "SMHW x SHX Frequency", "Events/yr")
plot_spatial_sf(res_bT_shx$freq, bottom_mask, "BMHW x SHX Frequency", "Events/yr")
plot_spatial_sf(res_triple$freq, bottom_mask, "Triple Compound Frequency", "Events/yr")


# PART 7: Regional Analyses ----------------------------------------------------

# ══════════════════════════════════════════════════════════════════════════════
## Step 1: Regional Masks: Gulf of Alaska vs CCLME, with north/south splits --------------
# CCLME bounding box: 20°N-50°N, 135°W-105°W
# Gulf of Alaska bounding box: 55°N-57°N, 144°W-130°W
# ══════════════════════════════════════════════════════════════════════════════

# Convert lon2d to -180:180 for readable bounding box logic
lon2d_180 <- ifelse(lon2d > 180, lon2d - 360, lon2d)

# CCLME box
CCLME_lat_min <- 20
CCLME_lat_max <- 47.9
CCLME_lon_min <- -135
CCLME_lon_max <- -105

CCLME_mask <- ocean_mask &
  lat2d >= CCLME_lat_min & lat2d <= CCLME_lat_max &
  lon2d_180 >= CCLME_lon_min & lon2d_180 <= CCLME_lon_max

# Gulf of Alaska box
GoA_lat_min <- 53.5
GoA_lat_max <- 61
GoA_lon_min <- -163
GoA_lon_max <- -136

GoA_mask <- ocean_mask &
  lat2d >= GoA_lat_min & lat2d <= GoA_lat_max &
  lon2d_180 >= GoA_lon_min & lon2d_180 <= GoA_lon_max

cat("CCLME cells:", sum(CCLME_mask), "\n") # 57,195
cat("GoA cells:  ", sum(GoA_mask), "\n") # 10,081


## Step 2: North/South splits within each region --------------------------------
# Using midpoint of each box's latitude range as the split line
CCLME_split_lat <- (CCLME_lat_min + CCLME_lat_max) / 2   # 35
GoA_split_lat   <- (GoA_lat_min + GoA_lat_max) / 2        # 56

CCLME_north_mask <- CCLME_mask & lat2d > CCLME_split_lat
CCLME_south_mask <- CCLME_mask & lat2d <= CCLME_split_lat

GoA_north_mask <- GoA_mask & lat2d > GoA_split_lat
GoA_south_mask <- GoA_mask & lat2d <= GoA_split_lat

cat("CCLME north (>", CCLME_split_lat, "N) cells:", sum(CCLME_north_mask), "\n") # 15,390
cat("CCLME south (<=", CCLME_split_lat, "N) cells:", sum(CCLME_south_mask), "\n") # 41,805
cat("GoA north (>", GoA_split_lat, "N) cells:", sum(GoA_north_mask), "\n") # 3454
cat("GoA south (<=", GoA_split_lat, "N) cells:", sum(GoA_south_mask), "\n") # 6627

cat("Test point in CCLME?", CCLME_mask[ilon, ilat], "\n")   # should be TRUE
cat("Test point in GoA?  ", GoA_mask[ilon, ilat], "\n")     # should be FALSE

# visual check
region_check <- ifelse(CCLME_mask, 1, ifelse(GoA_mask, 2, NA))
plot_spatial_sf(region_check, ocean_mask, "Region Mask: 1=CCLME, 2=GoA", "Region")

region_check_split <- ifelse(CCLME_north_mask, 1,
                             ifelse(CCLME_south_mask, 2,
                                    ifelse(GoA_north_mask, 3,
                                           ifelse(GoA_south_mask, 4, NA))))
plot_spatial_sf(region_check_split, ocean_mask,
                "Region Splits: 1=CCLME-N, 2=CCLME-S, 3=GoA-N, 4=GoA-S", "Region")


## Step 3: Regional Summary Stats ----------------------------------------------

summarize_region <- function(mat, mask) {
  vals <- mat[mask]
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(c(mean = NA, median = NA, sd = NA, n_cells = 0))
  c(mean = mean(vals), median = median(vals), sd = sd(vals), n_cells = length(vals))
}

regions <- list(
  CCLME       = CCLME_mask,
  CCLME_north = CCLME_north_mask,
  CCLME_south = CCLME_south_mask,
  GoA         = GoA_mask,
  GoA_north   = GoA_north_mask,
  GoA_south   = GoA_south_mask
)

# SMHW
region_summary_sst_freq     <- sapply(regions, function(m) summarize_region(out_freq, m & ocean_mask))
region_summary_sst_duration <- sapply(regions, function(m) summarize_region(out_mean_duration, m & ocean_mask))
region_summary_sst_intensity<- sapply(regions, function(m) summarize_region(out_mean_intensity, m & ocean_mask))

# BMHW
region_summary_bT_freq      <- sapply(regions, function(m) summarize_region(out_freq_bT, m & bottom_mask))
region_summary_bT_duration  <- sapply(regions, function(m) summarize_region(out_mean_duration_bT, m & bottom_mask))
region_summary_bT_intensity <- sapply(regions, function(m) summarize_region(out_mean_intensity_bT, m & bottom_mask))

# SHX
region_summary_shx_freq      <- sapply(regions, function(m) summarize_region(out_freq_shx, m & shx_mask))
region_summary_shx_duration  <- sapply(regions, function(m) summarize_region(out_mean_duration_shx, m & shx_mask))
region_summary_shx_intensity <- sapply(regions, function(m) summarize_region(out_mean_intensity_shx, m & shx_mask))

# Print everything
cat("\n===== SMHW Frequency by Region =====\n");      print(t(round(region_summary_sst_freq, 3)))
#              mean median    sd n_cells
# CCLME       0.168  0.154 0.046   57195
# CCLME_north 0.147  0.154 0.042   15390
# CCLME_south 0.176  0.185 0.045   41805
# GoA         0.130  0.123 0.037   10081
# GoA_north   0.138  0.154 0.039    3454
# GoA_south   0.126  0.123 0.035    6627
cat("\n===== SMHW Duration by Region =====\n");        print(t(round(region_summary_sst_duration, 3)))
#              mean median    sd n_cells
# CCLME       3.068  2.857 0.911   57195
# CCLME_north 2.677  2.286 1.003   15390
# CCLME_south 3.212  3.000 0.830   41805
# GoA         2.232  2.200 0.264   10081
# GoA_north   2.148  2.000 0.247    3454
# GoA_south   2.276  2.250 0.261    6627
cat("\n===== SMHW Intensity by Region =====\n");       print(t(round(region_summary_sst_intensity, 3)))
#              mean median    sd n_cells
# CCLME       1.918  1.908 0.175   57195
# CCLME_north 1.950  1.933 0.147   15390
# CCLME_south 1.907  1.898 0.183   41805
# GoA         1.923  1.917 0.153   10081
# GoA_north   1.927  1.911 0.144    3454
# GoA_south   1.921  1.921 0.158    6627

cat("\n===== BMHW Frequency by Region =====\n");       print(t(round(region_summary_bT_freq, 3)))
#              mean median    sd n_cells
# CCLME       0.190  0.185 0.066    3903
# CCLME_north 0.199  0.185 0.057     912
# CCLME_south 0.187  0.185 0.069    2991
# GoA         0.178  0.185 0.057    3624
# GoA_north   0.175  0.185 0.061    2178
# GoA_south   0.184  0.185 0.048    1446
cat("\n===== BMHW Duration by Region =====\n");        print(t(round(region_summary_bT_duration, 3)))
#              mean median    sd n_cells
# CCLME       5.355  4.714 2.622    3903
# CCLME_north 4.995  4.667 1.653     912
# CCLME_south 5.464  4.750 2.844    2991
# GoA         3.846  3.667 1.685    3624
# GoA_north   3.876  3.646 1.942    2178
# GoA_south   3.800  3.714 1.198    1446
cat("\n===== BMHW Intensity by Region =====\n");       print(t(round(region_summary_bT_intensity, 3)))
#              mean median    sd n_cells
# CCLME       1.778  1.766 0.274    3903
# CCLME_north 1.799  1.837 0.287     912
# CCLME_south 1.771  1.757 0.270    2991
# GoA         1.780  1.757 0.207    3624
# GoA_north   1.806  1.784 0.240    2178
# GoA_south   1.740  1.739 0.132    1446

cat("\n===== SHX Frequency by Region =====\n");        print(t(round(region_summary_shx_freq, 3)))
#              mean median    sd n_cells
# CCLME       0.271  0.278 0.073   57195
# CCLME_north 0.268  0.278 0.081   15390
# CCLME_south 0.272  0.278 0.070   41805
# GoA         0.095  0.093 0.079   10081
# GoA_north   0.072  0.000 0.094    3454
# GoA_south   0.106  0.123 0.067    6627
cat("\n===== SHX Duration by Region =====\n");         print(t(round(region_summary_shx_duration, 3)))
#              mean median    sd n_cells
# CCLME       3.900  3.600 1.473   57195
# CCLME_north 4.240  3.778 2.022   15390
# CCLME_south 3.774  3.556 1.185   41805
# GoA         6.475  6.667 6.166   10081
# GoA_north   2.827  0.000 3.812    3454
# GoA_south   8.376  8.200 6.302    6627

cat("\n===== SHX Intensity by Region =====\n");        print(t(round(region_summary_shx_intensity, 3)))
#              mean median    sd n_cells
# CCLME       1.622  1.617 0.292   57195
# CCLME_north 1.600  1.607 0.335   15390
# CCLME_south 1.630  1.619 0.275   41805
# GoA         1.111  1.532 0.822   10081
# GoA_north   0.705  0.000 0.883    3454
# GoA_south   1.323  1.607 0.700    6627


## Step 4: regional time series ------------------------------------------------
# domain mean SST/bT + % of region in MHW over tie 

# Fraction of region in SMHW/BMHW/SHX at each timestep
frac_in_region <- function(events_array, mask) {
  apply(events_array, 3, function(x) mean(x[mask] == 1, na.rm = TRUE))
}

df_regional <- data.frame(time = time_vec)

df_regional$sst_CCLME       <- frac_in_region(out_events_ts, CCLME_mask) * 100
df_regional$sst_CCLME_north <- frac_in_region(out_events_ts, CCLME_north_mask) * 100
df_regional$sst_CCLME_south <- frac_in_region(out_events_ts, CCLME_south_mask) * 100
df_regional$sst_GoA         <- frac_in_region(out_events_ts, GoA_mask) * 100
df_regional$sst_GoA_north   <- frac_in_region(out_events_ts, GoA_north_mask) * 100
df_regional$sst_GoA_south   <- frac_in_region(out_events_ts, GoA_south_mask) * 100

df_regional$bT_CCLME  <- frac_in_region(out_events_ts_bT, CCLME_mask & bottom_mask) * 100
df_regional$bT_GoA    <- frac_in_region(out_events_ts_bT, GoA_mask & bottom_mask) * 100

df_regional$shx_CCLME <- frac_in_region(out_events_ts_shx, CCLME_mask & shx_mask) * 100
df_regional$shx_GoA   <- frac_in_region(out_events_ts_shx, GoA_mask & shx_mask) * 100

## Plot: CCLME vs GoA SMHW extent over time
df_regional_long <- df_regional %>%
  select(time, sst_CCLME, sst_GoA) %>%
  pivot_longer(-time, names_to = "region", values_to = "pct_in_MHW")

ggplot(df_regional_long, aes(x = time, y = pct_in_MHW, color = region)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = c("sst_CCLME" = "darkred", "sst_GoA" = "steelblue"),
                     labels = c("CCLME", "Gulf of Alaska")) +
  labs(
    title = "SMHW Spatial Extent: CCLME vs Gulf of Alaska",
    y = "% of region in SMHW", x = NULL, color = "Region"
  )

## Plot: CCLME north vs south
df_CCLME_ns_long <- df_regional %>%
  select(time, sst_CCLME_north, sst_CCLME_south) %>%
  pivot_longer(-time, names_to = "region", values_to = "pct_in_MHW")

ggplot(df_CCLME_ns_long, aes(x = time, y = pct_in_MHW, color = region)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = c("sst_CCLME_north" = "orange", "sst_CCLME_south" = "purple"),
                     labels = c(paste0("CCLME North (>", CCLME_split_lat, "°N)"),
                                paste0("CCLME South (<=", CCLME_split_lat, "°N)"))) +
  labs(
    title = "SMHW Spatial Extent: CCLME North vs South",
    y = "% of sub-region in SMHW", x = NULL, color = "Region"
  )


# SAVE REGIONAL RESULTS

saveRDS(list(
  masks = list(
    CCLME = CCLME_mask, CCLME_north = CCLME_north_mask, CCLME_south = CCLME_south_mask,
    GoA = GoA_mask, GoA_north = GoA_north_mask, GoA_south = GoA_south_mask
  ),
  summaries = list(
    sst_freq = region_summary_sst_freq, sst_duration = region_summary_sst_duration,
    sst_intensity = region_summary_sst_intensity,
    bT_freq = region_summary_bT_freq, bT_duration = region_summary_bT_duration,
    bT_intensity = region_summary_bT_intensity,
    shx_freq = region_summary_shx_freq, shx_duration = region_summary_shx_duration,
    shx_intensity = region_summary_shx_intensity
  ),
  df_regional = df_regional
), file.path(output_dir, "regional_analysis.rds"))

regional <- readRDS(file.path(output_dir, "regional_analysis.rds"))
CCLME_mask       <- regional$masks$CCLME
CCLME_north_mask <- regional$masks$CCLME_north
CCLME_south_mask <- regional$masks$CCLME_south
GoA_mask         <- regional$masks$GoA
GoA_north_mask   <- regional$masks$GoA_north
GoA_south_mask   <- regional$masks$GoA_south

region_summary_sst_freq      <- regional$summaries$sst_freq
region_summary_sst_duration  <- regional$summaries$sst_duration
region_summary_sst_intensity <- regional$summaries$sst_intensity
region_summary_bT_freq       <- regional$summaries$bT_freq
region_summary_bT_duration   <- regional$summaries$bT_duration
region_summary_bT_intensity  <- regional$summaries$bT_intensity
region_summary_shx_freq      <- regional$summaries$shx_freq
region_summary_shx_duration  <- regional$summaries$shx_duration
region_summary_shx_intensity <- regional$summaries$shx_intensity

df_regional <- regional$df_regional

## Step 5: North/South of GoA and CCLME

# North of Gulf of Alaska (lat > 61°N)
north_of_GoA_mask <- ocean_mask &
  lat2d > GoA_lat_max &
  lon2d_180 >= GoA_lon_min & lon2d_180 <= GoA_lon_max

# South of Gulf of Alaska (lat < 54°N)
south_of_GoA_mask <- ocean_mask &
  lat2d < GoA_lat_min &   lon2d_180 >= GoA_lon_min & lon2d_180 <= GoA_lon_max

## North of CCLME (lat > 50°N)
north_of_CCLME_mask <- ocean_mask &
  lat2d > CCLME_lat_max &   lon2d_180 >= CCLME_lon_min & lon2d_180 <= CCLME_lon_max

# South of CCLME (lat < 20°N)
south_of_CCLME_mask <- ocean_mask &
  lat2d < CCLME_lat_min &   lon2d_180 >= CCLME_lon_min & lon2d_180 <= CCLME_lon_max

cat("North of GoA cells:  ", sum(north_of_GoA_mask), "\n") # 6,345
cat("South of GoA cells:  ", sum(south_of_GoA_mask), "\n") # 53,229
cat("North of CCLME cells:", sum(north_of_CCLME_mask), "\n") # 4,010
cat("South of CCLME cells:", sum(south_of_CCLME_mask), "\n") # 18,293

cat("Overlap N.GoA & N.CCLME:", sum(north_of_GoA_mask & north_of_CCLME_mask), "\n")
cat("Overlap S.GoA & CCLME box:", sum(south_of_GoA_mask & CCLME_mask), "\n")
cat("Overlap S.GoA & S.CCLME:", sum(south_of_GoA_mask & south_of_CCLME_mask), "\n")


gap_mask <- ocean_mask &
  lat2d > CCLME_lat_max & lat2d < GoA_lat_min &
  lon2d_180 >= min(CCLME_lon_min, GoA_lon_min) & lon2d_180 <= max(CCLME_lon_max, GoA_lon_max)

cat("Gap zone (50-54°N) cells:", sum(gap_mask), "\n")


# Visual check
region_check_beyond <- ifelse(north_of_GoA_mask, 1,
                              ifelse(south_of_GoA_mask, 2,
                                     ifelse(north_of_CCLME_mask, 3,
                                            ifelse(south_of_CCLME_mask, 4,
                                                   ifelse(gap_mask, 5,
                                                          ifelse(GoA_mask, 6,
                                                                 ifelse(CCLME_mask, 7, NA)))))))

plot_spatial_sf(region_check_beyond, ocean_mask,
                "Beyond-Box Regions: 1=N.GoA, 2=S.GoA, 3=N.CCLME, 4=S.CCLME, 5=Gap",
                "Region")

regions_beyond <- list(
  North_of_GoA   = north_of_GoA_mask,
  South_of_GoA   = south_of_GoA_mask,
  North_of_CCLME = north_of_CCLME_mask,
  South_of_CCLME = south_of_CCLME_mask
)

region_summary_sst_beyond <- sapply(regions_beyond, function(m) summarize_region(out_freq, m & ocean_mask))
region_summary_sst_dur_beyond <- sapply(regions_beyond, function(m) summarize_region(out_mean_duration, m & ocean_mask))
region_summary_sst_int_beyond <- sapply(regions_beyond, function(m) summarize_region(out_mean_intensity, m & ocean_mask))

region_summary_bT_beyond <- sapply(regions_beyond, function(m) summarize_region(out_freq_bT, m & bottom_mask))
region_summary_shx_beyond <- sapply(regions_beyond, function(m) summarize_region(out_freq_shx, m & shx_mask))

cat("\n===== SMHW Frequency: Beyond-Box Regions =====\n"); print(t(round(region_summary_sst_beyond, 3)))
cat("\n===== SMHW Duration: Beyond-Box Regions =====\n");  print(t(round(region_summary_sst_dur_beyond, 3)))
cat("\n===== SMHW Intensity: Beyond-Box Regions =====\n"); print(t(round(region_summary_sst_int_beyond, 3)))
cat("\n===== BMHW Frequency: Beyond-Box Regions =====\n"); print(t(round(region_summary_bT_beyond, 3)))
cat("\n===== SHX Frequency: Beyond-Box Regions =====\n");  print(t(round(region_summary_shx_beyond, 3)))


df_regional$sst_N_of_GoA <- frac_in_region(out_events_ts, north_of_GoA_mask) * 100
df_regional$sst_S_of_GoA <- frac_in_region(out_events_ts, south_of_GoA_mask) * 100
df_regional$sst_N_of_CCLME <- frac_in_region(out_events_ts, north_of_CCLME_mask) * 100
df_regional$sst_S_of_CCLME <- frac_in_region(out_events_ts, south_of_CCLME_mask) * 100

df_beyond_long <- df_regional %>%
  select(time, sst_N_of_GoA, sst_S_of_GoA, sst_N_of_CCLME, sst_S_of_CCLME) %>%
  pivot_longer(-time, names_to = "region", values_to = "pct_in_MHW")

ggplot(df_beyond_long, aes(x = time, y = pct_in_MHW, color = region)) +
  geom_line(linewidth = 0.75) +
  labs(
    title = "SMHW Spatial Extent: Regions Beyond GoA and CCLME Boxes",
    y = "% of sub-region in SMHW", x = NULL, color = "Region"
  )


output_dir <- 'C:/Users/mcclo/OneDrive/Documents/Hollings_Internship/Outputs'

expected_files <- c(
  "sst_mhw_summary.rds", "sst_mhw_events_ts.rds",
  "bT_mhw_summary.rds", "bT_mhw_events_ts.rds",
  "shx_summary.rds", "shx_events_ts.rds",
  "shx_point_and_params.rds", "shx_thresh_depth_map.rds",
  "compound_grid_summary.rds", "LMF_grid.rds",
  "regional_analysis.rds"
)

for (f in expected_files) {
  path <- file.path(output_dir, f)
  cat(f, ":", ifelse(file.exists(path), "EXISTS", "MISSING"), "\n")
}


