install.packages('zoo')

library(ncdf4)
library(terra)
library(tidync)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(zoo)

# ══════════════════════════════════════════════════════════════════════════════
# MHW DETECTION — SST AND BOTTOM TEMPERATURE
# Following Freeman et al. (2026) methodology
#
# CURRENT LIMITATIONS (data pending):
#   - Monthly resolution (Freeman uses daily; 11-day window + 5-day minimum duration
#     not achievable, adapted to 3-month window + 2-month min duration)
#   - No bathymetry (Freeman masks bottom MHWs to depth < 1000m — pending)
#   - No 3D oxygen (HLD/SHX not yet possible — pending)
#   - No chlorophyll (LCX not yet possible — pending)
#   - No IPSL projection files; incomplete UKESM files (hindcast only for now)
#   - Single test point only (spatial application pending)
# ══════════════════════════════════════════════════════════════════════════════


### Load in Data ----------------------------------------------------------------------------
list.files('D:/LorenzoRFM_ Hollings/Downloads/ROMs')


data_dir <- "D:/LorenzoRFM_ Hollings/Downloads/ROMs"

# extract yrs present
check_years <- function(prefix, dir = data_dir) {
  files <- list.files(dir, pattern = paste0("^", prefix), full.names = FALSE)
  years <- as.integer(gsub(".*_(\\d{4})\\.nc$", "\\1", files))
  years <- sort(years)
  expected <- seq(min(years), max(years))
  missing  <- setdiff(expected, years)
  cat("\n---", prefix, "---\n")
  cat("Years found:  ", min(years), "-", max(years), "(n =", length(years), ")\n")
  if (length(missing) > 0) {
    cat("MISSING YEARS:", missing, "\n")
  } else {
    cat("No missing years\n")
  }
}

check_years("sst_hindcast")
check_years("sst_GFDL")
check_years("sst_UKESM")
check_years("bO2_atl_gfdl_monavg")
check_years("bT_atl_gfdl_monavg")
# No more missing years!


### Explore the grid

sst_files <- sort(list.files(data_dir, pattern = "^sst_hindcast", full.names = TRUE))

nc <- nc_open(sst_files[1])

#get 2d lon lat arrays
lon2d <- ncvar_get(nc, 'longitude')
lat2d <- ncvar_get(nc, 'latitude')
mask <- ncvar_get(nc, 'mask')

# check dimensions
dim(lon2d)
dim(mask)


#coordinate range
cat("Lon range:", range(lon2d, na.rm = TRUE), "\n")
cat("Lat range:", range(lat2d, na.rm = TRUE), "\n")

# check time
time_raw <- ncvar_get(nc, "time")
time_units <- nc$dim$time$units
cat("Time units:", time_units, "\n")
cat("Time values:", time_raw, "\n")

# fix for time

parse_time <- function(file) {
  nc         <- nc_open(file)
  raw        <- ncvar_get(nc, "time")
  time_units <- nc$dim$time$units  # "seconds since 1900-01-01 00:00:00"
  nc_close(nc)
  
  origin <- as.POSIXct("1900-01-01 00:00:00", tz = "UTC")
  times  <- origin + as.numeric(raw)  # add seconds
  as.Date(times)
}

# Build full time vector across all 26 files
time_vec <- do.call(c, lapply(sst_files, parse_time))
cat("Time range:", as.character(range(time_vec)), "\n")
cat("Total timesteps:", length(time_vec), "\n") 

# getting SST for a file
sst_raw <- ncvar_get(nc, "sst")
dim(sst_raw)

nc_close(nc)

### finding a test point -------------------------------------------------------
target_lon <- -125
target_lat <- 40

# find nearest grid index
dist <- sqrt((lon2d - target_lon)^2 + (lat2d - target_lat)^2)
idx <- which(dist == min(dist, na.rm = TRUE), arr.ind = TRUE)
ilon <- idx[1] # xi_rho index
ilat <- idx[2] ## eta_rho index

cat("Nearest grid point found at xi =", ilon, ",eta =", ilat, "\n")
cat("Actual lon:", lon2d[ilon, ilat], "lat:", lat2d[ilon, ilat], "\n")
cat("Mask value:", mask[ilon, ilat], " (1 = ocean, 0 = land\n")


### extract time series at test point ------------------------------------------

extract_point_ts <- function(file_list, varname, ilon, ilat) {
  ts <- c()
  for (f in file_list) {
    nc  <- nc_open(f)
    # start = [xi, eta, time], count = [1, 1, -1] means all time steps
    val <- ncvar_get(nc, varname, start = c(ilon, ilat, 1), count = c(1, 1, -1))
    ts  <- c(ts, val)
    nc_close(nc)
  }
  return(ts)
}


sst_ts <- extract_point_ts(sst_files,varname = "sst", ilon = ilon, ilat = ilat)

cat("Length:", length(sst_ts), "\n") # 312
cat("NAs:", sum(is.na(sst_ts)), "\n") # 0
cat("SST Range:", range(sst_ts, na.rm = TRUE), "\n") # ~ 9.5 - 17.5


# plot, see if it looks real
plot(time_vec, sst_ts, type = "l",
     main = "SST at test point", ylab = "SST (deg C)", xlab = NULL)



### Build df -------------------------------------------------------------------
 df <- data.frame(
   time = time_vec,
   sst = sst_ts,
   month = as.integer(format(time_vec,"%m"))
 )


### Compute 90th percentile threshold ----------------------------------------------------------------------------


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


thresh_90 <- compute_threshold(df$sst, df$month, percentile = 0.90)

# should show seasonal patterns
print(data.frame(month = 1:12, threshold = round(thresh_90, 2)))

# back onto full time series
df$thresh_90 <- thresh_90[df$month]

### Flag extremes --------------------------------------------------------------

df$exceeds <- ifelse(df$sst > df$thresh_90, 1, 0)

cat("Months exceeding threshold:", sum(df$exceeds), "out of", nrow(df), "\n")
# 29/312 months 

### Merge events close together ------------------------------------------------

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

df$exceeds_merged <- merge_events(df$exceeds, max_gap = 1)

### remove short events --------------------------------------------------------

# Remove events shorter than minimum duration (Freeman et al. Section 2.3)
# Freeman: minimum 5 days
# Adaptation: Minimum 2 months

remove_short_events <- function(exceeds, min_duration = 2) {
  r <- rle(exceeds)
  r$values[r$values == 1 & r$lengths < min_duration] <- 0
  return(inverse.rle(r))
}

df$events <- remove_short_events(df$exceeds_merged, min_duration = 2)

cat("Months flagged as MHW after filtering:", sum(df$events), "\n")
## 22

### compute anomaly ------------------------------------------------------------

# Climatological z-score anomaly (Freeman et al. Section 2.3)
# Anomaly computed relative to monthly climatological mean and sd
# NOT the overall time series mean. removes seasonal cycle

compute_anomaly <- function(ts, months) {
  clim_mean <- tapply(ts, months, mean, na.rm = TRUE)
  clim_sd   <- tapply(ts, months, sd,   na.rm = TRUE)
  (ts - clim_mean[months]) / clim_sd[months]
}

df$anomaly <- compute_anomaly(df$sst, df$month)

### characterize events --------------------------------------------------------

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

event_summary <- characterize_events(
  events_flag = df$events,
  anomaly = df$anomaly,
  time_vec = df$time
)

print(event_summary)
table(df$events)

r      <- rle(df$events)
ends   <- cumsum(r$lengths)
starts <- ends - r$lengths + 1

events_table <- data.frame()
for (i in seq_along(r$values)) {
  if (r$values[i] == 1) {
    idx      <- starts[i]:ends[i]
    events_table <- rbind(events_table, data.frame(
      start           = df$time[starts[i]],
      end             = df$time[ends[i]],
      duration_months = length(idx),
      mean_intensity  = round(mean(df$anomaly[idx], na.rm = TRUE), 3)
    ))
  }
}

n_years <- as.numeric(diff(range(df$time))) / 365.25

location_freq <- round(nrow(events_table) / n_years, 2)
cat("Event frequency at this location:", location_freq, "events/year\n")
# 0.23 freq per yr

print(events_table)
#        start        end duration_months mean_intensity
# 1 1997-08-16 1997-09-16               2          2.862
# 2 1997-12-16 1998-03-17               4          0.593
# 3 2014-09-16 2015-03-18               7          1.417
# 4 2015-07-17 2015-11-16               5          2.171
# 5 2016-03-17 2016-04-16               2          0.033
# 6 2020-05-17 2020-06-16               2          0.912


# ══════════════════════════════════════════════════════════════════════════════
# SURFACE MHW (SMHW) — SST
# Freeman Section 2.3: extreme temp event >= 5 days exceeding P90
# Applied across full domain (no depth restriction for surface MHW)
# ══════════════════════════════════════════════════════════════════════════════

# Build df
df <- data.frame(
  time = time_vec,
  sst = sst_ts,
  month = as.integer(format(time_vec, "%m"))
)

### Step 1: Compute seasonally-varying 90th percentile threshold -----------------
# Freeman: P90 for MHW (Table 1). trends not removed before detection (Section 2.3)
thresh_90 <- compute_threshold(df$sst, df$month, percentile = 0.90, window = 1)
df$thresh_90 <- thresh_90[df$month]

### Step 2: Flag months exceeding threshold ------------------------------------
# Freeman: SST > P90
df$exceeds <- ifelse(df$sst > df$thresh_90, 1, 0)
cat("SST — months exceeding P90:", sum(df$exceeds), "out of", nrow(df), "\n")
# 29/312 months exceeding

### Step 3: Merge events close in time -----------------------------------------
# Freeman: merge if gap <= 3 days | ADAPTATION: <= 1 month
df$exceeds_merged <- merge_events(df$exceeds, max_gap = 1)

### Step 4: Remove short events --------------------------------------------------
# Freeman: remove if < 5 days | ADAPTATION: remove if < 2 months
df$events <- remove_short_events(df$exceeds_merged, min_duration = 2)
cat("SST — months flagged as SMHW after filtering:", sum(df$events), "\n")
# 22 months flagged

# Step 5: Compute climatological z-score anomaly for intensity -----------------
# Freeman Section 2.4: intensity = mean anomaly over event period
df$anomaly <- compute_anomaly(df$sst, df$month)

# Step 6: Characterize events --------------------------------------------------
sst_events <- characterize_events(df$events, df$anomaly, df$time)
print(sst_events)
#    1  32  34  36  40 237 244 247 252 255 257 305 307

### plot sst -----------------------------------------------------------------------
#df$mhw_shade <- ifelse(df$events == 1, df$sst, NA)

ggplot(df, aes(x = time)) +
  geom_rect(
    data = df %>% filter(events == 1),
    aes(xmin = time - 15, xmax = time + 15,  # ~half month either side
        ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.5
  ) +
  # SST line on top
  geom_line(aes(y = sst), color = "black", linewidth = 0.75) +
  # Seasonally-varying P90 threshold (Freeman Section 2.3)
  geom_line(aes(y = thresh_90), color = "black",
            linetype = "dashed", linewidth = 0.5) +
  labs(
    title    = paste("MHW Detection — lon:", round(lon2d[ilon, ilat], 2),
                     "lat:", round(lat2d[ilon, ilat], 2)),
    subtitle = "Dashed = 90th percentile | Red shading = MHW event",
    y        = "SST (°C)",
    x        = NULL
  )

# ══════════════════════════════════════════════════════════════════════════════
# PART 2: BOTTOM MHW (BMHW) — BOTTOM TEMPERATURE
# Freeman Section 2.3: same P90 threshold logic as SMHW
# Freeman: restricted to locations where bottom depth < 1000m
# LIMITATION: depth mask not yet applied, because no bathymetry
# ══════════════════════════════════════════════════════════════════════════════

### Bottom temp MHW (NO DEPTH MASK FOR NOW) -----------------------------------------------------------

#### Check file contents ----------------------------------------------------------

bT_files <- sort(list.files(data_dir, pattern = "^bT_atl_gfdl_monavg_", full.names = TRUE))

# keep in hindcast yrs only
bT_files <- bT_files[grepl("(199[5-9]|200[0-9]|201[0-9]|2020)\\.nc$", bT_files)]
cat("# of bT files:", length(bT_files), "\n")

# check contents
nc_bT <- nc_open(bT_files[1]) 
names(nc_bT$var) # lon_rho, lat_rho, temp
names(nc_bT$dim) # ocean_time xi_rho eta_rho s_rho
nc_close(nc_bT)

#### Load grid and bathymetry

nc_bT <- nc_open(bT_files[1])
lon2d_bT <- ncvar_get(nc_bT, "lon_rho") 
lat2d_bT <- ncvar_get(nc_bT, "lat_rho")   

# what is s_rho?
nc_bT$dim$s_rho$len
# 1 depth level?
nc_bT$dim$s_rho$vals
# -0.99

# check time dimension name
time_raw_bT <- ncvar_get(nc_bT, "ocean_time")
time_units_bT <- nc_bT$dim$ocean_time$units
cat("Time units:", time_units_bT, "\n")
# Time units: seconds since 1900-01-01 00:00:00 

#check dimensions of temp
temp_test <- ncvar_get(nc_bT, "temp", start = c(1, 1, 1, 1), count = c(1, 1, -1, 1))
cat("Number of depth layers:", length(temp_test), "\n")
# 1 depth layer, bottom only, 2D data

nc_close(nc_bT)

# Update parse_time for ocean_time
parse_time_bT <- function(file) {
  nc <- nc_open(file)
  raw <- ncvar_get(nc, "ocean_time")
  time_units <- nc$dim$ocean_time$units
  nc_close(nc)
  
  origin <- as.POSIXct("1900-01-01 00:00:00", tz = "UTC")
  times <- origin + as.numeric(raw)
  as.Date(times)
}

time_vec_bT <- do.call(c, lapply(bT_files, parse_time_bT))
cat("Time range:", as.character(range(time_vec_bT)), "\n")
# Time range: 1995-01-16 2020-12-16 , correct


#### Find test point
target_lon <- -125
target_lat <- 40

dist_bT <- sqrt((lon2d_bT - target_lon)^2 + (lat2d_bT - target_lat)^2)
# dist_bT[mask_bT == 0] <- NA # excluding only land for rn

idx_bT <- which(dist_bT == min(dist_bT, na.rm = TRUE), arr.ind = TRUE)
ilon_bT <- idx_bT[1]
ilat_bT <- idx_bT[2]

cat("Test point: xi =", ilon_bT, ", eta =", ilat_bT, "\n")
cat("Lon:", lon2d_bT[ilon_bT, ilat_bT],
    "Lat:", lat2d_bT[ilon_bT, ilat_bT], "\n")
# cat("Mask:", mask_bT[ilon_bT, ilat_bT], "\n")


#### Extract time series -------------------------------------------------------

extract_point_ts_bT <- function(file_list, varname, ilon, ilat) {
  ts <- c()
  for (f in file_list) {
    nc  <- nc_open(f)
    val <- ncvar_get(nc, varname,
                     start = c(ilon, ilat, 1, 1),
                     count = c(1, 1, 1, -1))
    ts  <- c(ts, as.numeric(val))
    nc_close(nc)
  }
  return(ts)
}

bT_ts <- extract_point_ts_bT(bT_files,
                             varname = "temp",
                             ilon = ilon_bT,
                             ilat = ilat_bT)

time_vec_bT <- do.call(c, lapply(bT_files, parse_time_bT))

cat("Length:", length(bT_ts), "\n") # 312     
cat("NAs:", sum(is.na(bT_ts)), "\n") # 0             
cat("Temp range:", range(bT_ts, na.rm = TRUE), "\n")  # ~2.137 - 2.6240

plot(time_vec_bT, bT_ts, type = "l",
     main = "Bottom temp at test point",
     ylab = "Temp (deg C)", xlab = NULL)

#### Build DF ------------------------------------------------------------------

df_bT <- data.frame(
  time = time_vec_bT,
  bT = bT_ts,
  month = as.integer(format(time_vec_bT, "%m"))
)

#### Step 1: Compute seasonally-varying 90th percentile threshold -------------------------------------------------
# same as SMHW, P90 3-month window
thresh_90_bT <- compute_threshold(df_bT$bT,
                                  df_bT$month,
                                  percentile = 0.90,
                                  window = 1)

print(data.frame(month = 1:12, bT_threshold = round(thresh_90_bT, 2)))
#    month bT_threshold
# 1      1         2.46
# 2      2         2.40
# 3      3         2.39
# 4      4         2.44
# 5      5         2.50
# 6      6         2.55
# 7      7         2.54
# 8      8         2.46
# 9      9         2.40
# 10    10         2.39
# 11    11         2.44
# 12    12         2.46


df_bT$thresh_90 <- thresh_90_bT[df_bT$month]

#### Step 2: flag exceedances/extremes ------------------------------------------------------

df_bT$exceeds <- ifelse(df_bT$bT > df_bT$thresh_90, 1, 0)

cat("Months exceeding threshold:", sum(df_bT$exceeds),
    "out of", nrow(df_bT), "\n")
# 28 out of 312

#### Step 3: Merge close events --------------------------------------------------------

df_bT$exceeds_merged <- merge_events(df_bT$exceeds, max_gap = 1)

#### Step 4: remove short events ------------------------------------------------------

df_bT$events_bT <-remove_short_events(df_bT$exceeds_merged, min_duration = 2)

cat("Months flagged as bottom MHW:", sum(df_bT$events_bT), "\n")
# 29

#### Step 5: characterize events -------------------------------------------------------
# characterize_events function in the works

df_bT$anomaly <- compute_anomaly(df_bT$bT, df_bT$month)

r_bT <- rle(df_bT$events_bT)
ends_bT <- cumsum(r_bT$lengths)
starts_bT <- ends_bT - r_bT$lengths + 1

bT_events <- data.frame() 
for (i in seq_along(r_bT$values)) {
  if (r_bT$values[i] == 1) {
    idx <- starts_bT[i]:ends_bT[i]
    bT_events <- rbind(bT_events, data.frame(
      start = df_bT$time[starts_bT[i]],
      end = df_bT$time[ends_bT[i]],
      duration_months = length(idx),
      mean_intensity = round(mean(df_bT$anomaly[idx], na.rm = TRUE), 3)
    ))
  }
}

n_years_bT <- as.numeric(diff(range(df_bT$time))) / 365.25
if (nrow(bT_events) > 0) {
  bT_events$freq_per_year <- round(nrow(bT_events) / n_years_bT, 2)
}

print(bT_events)
#        start        end duration_months mean_intensity freq_per_year
# 1 2010-10-16 2010-11-16               2          1.217          0.27
# 2 2014-01-16 2014-06-17               6          1.518          0.27
# 3 2015-07-17 2015-09-16               3          1.480          0.27
# 4 2016-05-17 2016-08-16               4          2.075          0.27
# 5 2017-10-16 2018-01-16               4          1.856          0.27
# 6 2018-12-16 2019-01-16               2          1.563          0.27
# 7 2019-12-16 2020-07-17               8          2.173          0.27


#### Plot bottom MHW ---------------------------------------------------------------------

ggplot(df_bT, aes(x = time)) +
  geom_rect(
    data = df_bT %>% filter(events_bT == 1),
    aes(xmin = time - 15, xmax = time + 15,
        ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.5
  ) +
  geom_line(aes(y = bT), color = "black", linewidth = 0.75) +
  geom_line(aes(y = thresh_90), color = "black",
            linetype = "dashed", linewidth = 0.5) +
  labs(
    title    = paste("Bottom MHW Detection — lon:",
                     round(lon2d_bT[ilon_bT, ilat_bT], 2),
                     "lat:", round(lat2d_bT[ilon_bT, ilat_bT], 2)),
    subtitle = "Dashed = 90th percentile | Red = bottom MHW | depth mask is not available",
    y        = "Bottom Temperature (°C)",
    x        = NULL
  )

# ══════════════════════════════════════════════════════════════════════════════
# PART 3: COMPOUND SMHW + BMHW
# Freeman Section 2.3: compound = 2+ extremes co-occurring for >= 5 days
# Freeman Section 2.4: compound intensity = Z_var1 * Z_var2, then standardized
# ADAPTATION: co-occurrence >= 2 months (daily resolution pending)
# ══════════════════════════════════════════════════════════════════════════════

#### Compare sst and bt events -------------------------------------------------

df_compare <- data.frame(
  time = df$time,
  sst = df$sst,
  bT = df_bT$bT,
  events_SST = df$events,
  events_bT = df_bT$events_bT
)

### Step 1: flag raw co-occurrences ----------------------------------------------
compound_raw <- ifelse(
  df_compare$events_SST == 1 & df_compare$events_bT == 1, 1, 0
)

### Step 2: apply min duration filter to compound events (Freeman 2.3) ----------
df_compare$compound <- remove_short_events(compound_raw, min_duration = 2)

### Step 3: compound intensity = Z_SST * Z_bT, standardized(freeman 2.4) -------------------------------------------------------------------
z_SST <- compute_anomaly(df_compare$sst, as.integer(format(df_compare$time, "%m")))
z_bT  <- compute_anomaly(df_compare$bT,  as.integer(format(df_compare$time, "%m")))

raw_intensity <- z_SST * z_bT
df_compare$compound_intensity <- (raw_intensity - mean(raw_intensity, na.rm = TRUE)) /
  sd(raw_intensity, na.rm = TRUE)

cat("Compound SMHW + BMHW months:", sum(df_compare$compound), "\n")
# 5 compound months


### Step 4: Plot Compound Events ---------------------------------------------------------------
ggplot(df_compare, aes(x = time)) +
  geom_rect(
    data = df_compare %>% filter(events_SST == 1),
    aes(xmin = time - 15, xmax = time + 15,
        ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.3
  ) +
  geom_rect(
    data = df_compare %>% filter(events_bT == 1),
    aes(xmin = time - 15, xmax = time + 15,
        ymin = -Inf, ymax = Inf),
    fill = "blue", alpha = 0.3
  ) +
  # Explicitly outline compound periods
  geom_rect(
    data = df_compare %>% filter(compound == 1),
    aes(xmin = time - 15, xmax = time + 15,
        ymin = -Inf, ymax = Inf),
    fill = NA, color = "magenta", linewidth = 0.9
  ) +
  geom_line(aes(y = scale(sst)[,1],  color = "SST"),         linewidth = 0.75) +
  geom_line(aes(y = scale(bT)[,1],   color = "Bottom Temp"), linewidth = 0.75) +
  scale_color_manual(values = c("SST" = "red", "Bottom Temp" = "blue")) +
  labs(
    title    = "Compound SMHW + BMHW",
    subtitle = paste(
      "Freeman et al. method | Red = SMHW | Blue = BMHW | Magenta outline = compound\n",
      "Compound intensity = Z_SST x Z_bT (standardized) | depth mask + daily resolution pending"
    ),
    y     = "Standardized Anomaly (z-score)",
    x     = NULL,
    color = NULL
  )






















