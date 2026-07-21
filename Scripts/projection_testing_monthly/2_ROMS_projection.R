library(ncdf4)
library(terra)
library(tidync)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(zoo)
library(matrixStats)
library(patchwork)


# Preliminary testing ----------------------------------------------------------

# 1. time units in the SST file — dimension is named "time" here, not "ocean_time"
#    like the oxygen file, so don't assume it matches
nc_sst <- nc_open("D:/LorenzoRFM_ Hollings/Downloads/ROMs/sst_GFDL_1995.nc")
cat("SST time units:", nc_sst$dim$time$units, "\n")

# 2. does "mask" in the SST file share the SST grid shape (176 x 536)?
mask_test <- ncvar_get(nc_sst, "mask")
cat("mask dims:", dim(mask_test), "\n")
nc_close(nc_sst)

# 3. check the bT (bottom temp) file's grid too — is it on the O2 grid,
#    the SST grid, or a third one?
nc_bt <- nc_open("D:/LorenzoRFM_ Hollings/Downloads/ROMs/bT_atl_gfdl_monavg_1995.nc")
cat("bT dims:\n")
print(sapply(nc_bt$dim, function(d) d$len))
lon_bt <- ncvar_get(nc_bt, "lon_rho")
lat_bt <- ncvar_get(nc_bt, "lat_rho")
cat("bT lon range:", range(lon_bt, na.rm=TRUE), " lat range:", range(lat_bt, na.rm=TRUE), "\n")
nc_close(nc_bt)

nc_bt <- nc_open("D:/LorenzoRFM_ Hollings/Downloads/ROMs/bT_atl_gfdl_monavg_1995.nc")
print(names(nc_bt$var))   # does "mask" show up here too?
nc_close(nc_bt)

nc_o2 <- nc_open("D:/LorenzoRFM_ Hollings/Downloads/ROMs/bO2_atl_gfdl_monavg_1995.nc")
print(names(nc_o2$var))
nc_close(nc_o2)
nc_bt$var$temp$missval


# ══════════════════════════════════════════════════════════════════════════════
# MHW/BMHW/BHX PROJECTION ANALYSIS — ROMS VERSION -------------------------------
# Adapted from hindcast script (Freeman et al. methodology)

# CONFIRMED VIA DIRECT INSPECTION (do not re-derive from filenames — verified):
#   - Hindcast:      MOM6, dims ih/jh = 342 x 816, grid in ocean_static.nc,
#                    time = "days since 1993-01-01"
#   - Proj. SST:     ROMS, dims xi_rho/eta_rho = 176 x 536, own longitude/
#                    latitude vars, explicit "mask" variable, var name "sst",
#                    time = "seconds since 1900-01-01"
#   - Proj. bT/bO2:  ROMS, dims xi_rho/eta_rho/s_rho = 176 x 364 x 1 (s_rho
#                    is a singleton — confirmed single-level, NOT a profile;
#                    SHX/HLD is therefore still not reproducible, BHX stands),
#                    own lon_rho/lat_rho vars, NO mask variable — land cells
#                    are the fill value (missval = 1e+37) which ncdf4 auto-
#                    converts to NA on read, var names "temp" (bT) / "oxygen"
#                    (bO2), time = "seconds since 1900-01-01"
#   - bT and bO2 SHARE the same grid (identical shape + lon/lat range)
#   - SST grid is LARGER and does NOT share indices with bT/bO2 grid, despite
#     overlapping domains — regridding (nearest-neighbor) required for any
#     analysis combining SST with bT or bO2 (i.e. SMHW x BMHW, SMHW x BHX).
#   - The MOM6 hindcast grid is a THIRD, independent grid. Fixed-baseline mode
#     (threshold/climatology from hindcast, applied to ROMS projections)
#     requires its own nearest-neighbor mapping from hindcast cells to each
#     ROMS grid.
# ══════════════════════════════════════════════════════════════════════════════

data_dir <- "D:/LorenzoRFM_ Hollings/Downloads/ROMs"

# ══════════════════════════════════════════════════════════════════════════════
# PART 0: Extreme event functions (same as hindcast) ---------------------------
# ══════════════════════════════════════════════════════════════════════════════

compute_threshold <- function(ts, months, percentile = 0.90, window = 1) {
  thresh <- numeric(12)
  for (m in 1:12) {
    window_months <- ((m - 1 - window):(m - 1 + window)) %% 12 + 1
    idx <- which(months %in% window_months)
    thresh[m] <- quantile(ts[idx], probs = percentile, na.rm = TRUE)
  }
  return(thresh)
}

merge_events <- function(exceeds, max_gap = 1) {
  extreme_idx <- which(exceeds == 1)
  if (length(extreme_idx) == 0) return(exceeds)
  gaps <- diff(extreme_idx)
  fill_positions <- which(gaps > 1 & gaps <= max_gap + 1)
  for (pos in fill_positions) {
    start_fill <- extreme_idx[pos] + 1
    end_fill   <- extreme_idx[pos + 1] - 1
    if (start_fill <= end_fill) exceeds[start_fill:end_fill] <- 1
  }
  return(exceeds)
}

remove_short_events <- function(exceeds, min_duration = 2) {
  r <- rle(exceeds)
  r$values[r$values == 1 & r$lengths < min_duration] <- 0
  return(inverse.rle(r))
}

compute_anomaly <- function(ts, months) {
  clim_mean <- tapply(ts, months, mean, na.rm = TRUE)
  clim_sd   <- tapply(ts, months, sd,   na.rm = TRUE)
  (ts - clim_mean[months]) / clim_sd[months]
}

characterize_events <- function(events_flag, anomaly, time_vec) {
  r <- rle(events_flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  out <- data.frame()
  for (i in seq_along(r$values)) {
    if (r$values[i] == 1) {
      idx <- starts[i]:ends[i]
      out <- rbind(out, data.frame(
        start = time_vec[starts[i]], end = time_vec[ends[i]],
        duration_months = length(idx),
        mean_intensity = round(mean(anomaly[idx], na.rm = TRUE), 3)
      ))
    }
  }
  n_years <- as.numeric(diff(range(time_vec))) / 365.25
  cat("Event frequency at this location:", round(nrow(out) / n_years, 2), "events/year\n")
  return(out)
}

# fixed-baseline variants (ref series fits the threshold/climatology, target gets flagged)
compute_threshold_fixed <- function(ref_ts, ref_months, percentile = 0.90, window = 1) {
  compute_threshold(ref_ts, ref_months, percentile, window)
}
compute_anomaly_fixed <- function(ref_ts, ref_months, target_ts, target_months) {
  clim_mean <- tapply(ref_ts, ref_months, mean, na.rm = TRUE)
  clim_sd   <- tapply(ref_ts, ref_months, sd,   na.rm = TRUE)
  (target_ts - clim_mean[target_months]) / clim_sd[target_months]
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 1: Grid utilities (ROMS-specific) ---------------------------------------
# ══════════════════════════════════════════════════════════════════════════════

# Parse ROMS time (seconds since 1900), vs. hindcast's days-since parser
parse_time_roms <- function(file, timevar) {
  nc <- nc_open(file)
  raw <- ncvar_get(nc, timevar)
  time_units <- nc$dim[[timevar]]$units
  nc_close(nc)
  origin <- as.Date(sub("seconds since ", "", time_units))
  origin + as.numeric(raw) / 86400
}

# Find nearest single grid cell to a target real-world lon/lat, on ANY grid
find_nearest_cell <- function(lon2d, lat2d, mask, target_lon_180, target_lat) {
  # grids here use -180:180 convention (longitude/lon_rho), unlike MOM6's 0:360
  dist <- sqrt((lon2d - target_lon_180)^2 + (lat2d - target_lat)^2)
  dist[!mask] <- NA
  idx <- which(dist == min(dist, na.rm = TRUE), arr.ind = TRUE)[1, ]
  list(i = as.integer(idx[1]), j = as.integer(idx[2]))
}

# Build a nearest-neighbor index mapping from a SOURCE grid onto a TARGET grid,
# e.g. mapping every bT/bO2 cell to its nearest SST cell, so event arrays from
# different grids can be combined for compound analysis. Uses coarse lon/lat
# binning instead of brute-force all-pairs distance, since e.g. 176x364 vs
# 176x536 brute force (~23M pairs) is slow but binning keeps it fast.
# Returns a data.frame: source_i, source_j, target_i, target_j (only for
# source cells that are within max_dist_deg of some target cell).
build_nn_mapping <- function(src_lon, src_lat, src_mask,
                             tgt_lon, tgt_lat, tgt_mask,
                             bin_deg = 1, max_dist_deg = 0.5) {
  tgt_idx <- which(tgt_mask, arr.ind = TRUE)
  tgt_lons <- tgt_lon[tgt_mask]
  tgt_lats <- tgt_lat[tgt_mask]
  
  # bin target cells for fast local lookup
  tgt_bin_x <- floor(tgt_lons / bin_deg)
  tgt_bin_y <- floor(tgt_lats / bin_deg)
  tgt_key   <- paste(tgt_bin_x, tgt_bin_y)
  tgt_by_bin <- split(seq_along(tgt_lons), tgt_key)
  
  src_idx <- which(src_mask, arr.ind = TRUE)
  n_src <- nrow(src_idx)
  out <- data.frame(source_i = integer(n_src), source_j = integer(n_src),
                    target_i = integer(n_src), target_j = integer(n_src))
  
  for (k in seq_len(n_src)) {
    si <- src_idx[k, 1]; sj <- src_idx[k, 2]
    slon <- src_lon[si, sj]; slat <- src_lat[si, sj]
    
    bx <- floor(slon / bin_deg); by <- floor(slat / bin_deg)
    candidates <- c()
    for (dx in -1:1) for (dy in -1:1) {
      key <- paste(bx + dx, by + dy)
      if (!is.null(tgt_by_bin[[key]])) candidates <- c(candidates, tgt_by_bin[[key]])
    }
    if (length(candidates) == 0) {
      out[k, ] <- c(si, sj, NA, NA); next
    }
    d <- sqrt((tgt_lons[candidates] - slon)^2 + (tgt_lats[candidates] - slat)^2)
    best <- candidates[which.min(d)]
    if (min(d) > max_dist_deg) {
      out[k, ] <- c(si, sj, NA, NA)
    } else {
      out[k, ] <- c(si, sj, tgt_idx[best, 1], tgt_idx[best, 2])
    }
  }
  out
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 2: Load grid definitions ------------------------------------------------
# ══════════════════════════════════════════════════════════════════════════════

## MOM6 hindcast grid (unchanged from hindcast script) 
grid_file_hcast <- "D:/LorenzoRFM_ Hollings/Downloads/ocean_static.nc"
nc_grid <- nc_open(grid_file_hcast)
hcast_lon2d <- ncvar_get(nc_grid, "geolon")   # 0:360 convention
hcast_lat2d <- ncvar_get(nc_grid, "geolat")
hcast_depth <- ncvar_get(nc_grid, "deptho")
nc_close(nc_grid)
hcast_lon2d_180 <- ifelse(hcast_lon2d > 180, hcast_lon2d - 360, hcast_lon2d)
hcast_mask <- !is.na(hcast_depth) & hcast_depth > 0

## -- ROMS SST grid
nc_sst0 <- nc_open(file.path(data_dir, "sst_GFDL_1995.nc"))
sst_lon <- ncvar_get(nc_sst0, "longitude")     # already -180:180
sst_lat <- ncvar_get(nc_sst0, "latitude")
sst_mask_raw <- ncvar_get(nc_sst0, "mask")     # CONFIRM: 1=ocean/0=land or reverse — check values
nc_close(nc_sst0)
cat("SST mask unique values:", unique(as.vector(sst_mask_raw)), "\n") 
sst_mask <- sst_mask_raw == 1  

## -- ROMS bottom grid (bT/bO2 share this) --
nc_bt0 <- nc_open(file.path(data_dir, "bT_atl_gfdl_monavg_1995.nc"))
bottom_lon <- ncvar_get(nc_bt0, "lon_rho")     # already -180:180
bottom_lat <- ncvar_get(nc_bt0, "lat_rho")
temp0 <- ncvar_get(nc_bt0, "temp", start = c(1, 1, 1, 1), count = c(-1, -1, 1, 1))
nc_close(nc_bt0)
bottom_mask <- !is.na(temp0)  

# ══════════════════════════════════════════════════════════════════════════════
# PART 3: Test point, separately on each grid ----------------------------------
# Same real-world location (-122.02, 36.7) as the hindcast script, found
# independently on all three grids since indices don't correspond across them.
# ══════════════════════════════════════════════════════════════════════════════

target_lon_180 <- -122.02
target_lat     <- 36.7

hcast_pt  <- find_nearest_cell(hcast_lon2d_180, hcast_lat2d, hcast_mask, target_lon_180, target_lat)
sst_pt    <- find_nearest_cell(sst_lon, sst_lat, sst_mask, target_lon_180, target_lat)
bottom_pt <- find_nearest_cell(bottom_lon, bottom_lat, bottom_mask, target_lon_180, target_lat)

cat("Hindcast cell:", hcast_pt$i, hcast_pt$j,
    "| actual lonlat:", hcast_lon2d_180[hcast_pt$i, hcast_pt$j], hcast_lat2d[hcast_pt$i, hcast_pt$j], "\n")
cat("SST cell:", sst_pt$i, sst_pt$j,
    "| actual lonlat:", sst_lon[sst_pt$i, sst_pt$j], sst_lat[sst_pt$i, sst_pt$j], "\n")
cat("Bottom cell:", bottom_pt$i, bottom_pt$j,
    "| actual lonlat:", bottom_lon[bottom_pt$i, bottom_pt$j], bottom_lat[bottom_pt$i, bottom_pt$j], "\n")

# ══════════════════════════════════════════════════════════════════════════════
# PART 4: Point-level loaders (ROMS yearly files) ------------------------------
# ══════════════════════════════════════════════════════════════════════════════

load_yearly_point_sst <- function(data_dir, years, i, j) {
  ts_list <- vector("list", length(years)); time_list <- vector("list", length(years))
  for (k in seq_along(years)) {
    f <- file.path(data_dir, sprintf("sst_GFDL_%d.nc", years[k]))
    if (!file.exists(f)) { warning("Missing: ", f); next }
    nc <- nc_open(f)
    val <- ncvar_get(nc, "sst", start = c(i, j, 1), count = c(1, 1, -1))
    ts_list[[k]] <- as.numeric(val)
    time_list[[k]] <- parse_time_roms(f, "time")
    nc_close(nc)
  }
  ord <- order(do.call(c, time_list))
  list(time = do.call(c, time_list)[ord], ts = do.call(c, ts_list)[ord])
}

load_yearly_point_bottom <- function(data_dir, prefix_fmt, varname, years, i, j) {
  ts_list <- vector("list", length(years)); time_list <- vector("list", length(years))
  for (k in seq_along(years)) {
    f <- file.path(data_dir, sprintf(prefix_fmt, years[k]))
    if (!file.exists(f)) { warning("Missing: ", f); next }
    nc <- nc_open(f)
    # dims are xi_rho, eta_rho, s_rho, ocean_time — s_rho is a singleton, squeeze it
    val <- ncvar_get(nc, varname, start = c(i, j, 1, 1), count = c(1, 1, 1, -1))
    ts_list[[k]] <- as.numeric(val)
    time_list[[k]] <- parse_time_roms(f, "ocean_time")
    nc_close(nc)
  }
  ord <- order(do.call(c, time_list))
  list(time = do.call(c, time_list)[ord], ts = do.call(c, ts_list)[ord])
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 5: Point-level pipeline (same as hind) ----------------------------------
# ══════════════════════════════════════════════════════════════════════════════

run_point_pipeline <- function(time_vec, ts, direction = c("high", "low"),
                               mode = c("moving", "fixed"),
                               fixed_thresh = NULL, ref_ts = NULL, ref_months = NULL,
                               percentile = 0.90, window = 1,
                               max_gap = 1, min_duration = 2, label = "") {
  direction <- match.arg(direction); mode <- match.arg(mode)
  months <- as.integer(format(time_vec, "%m"))
  
  if (mode == "moving") {
    thresh_vec <- compute_threshold(ts, months, percentile, window)
    anomaly    <- compute_anomaly(ts, months)
  } else {
    stopifnot(!is.null(fixed_thresh), !is.null(ref_ts), !is.null(ref_months))
    thresh_vec <- fixed_thresh
    anomaly    <- compute_anomaly_fixed(ref_ts, ref_months, ts, months)
  }
  thresh_ts <- thresh_vec[months]
  exceeds <- if (direction == "high") ifelse(ts > thresh_ts, 1, 0) else ifelse(ts < thresh_ts, 1, 0)
  exceeds[is.na(exceeds)] <- 0
  events <- remove_short_events(merge_events(exceeds, max_gap), min_duration)
  
  cat("---", label, "| mode:", mode, "---\n")
  events_df <- characterize_events(events, anomaly, time_vec)
  list(time = time_vec, ts = ts, months = months, thresh_ts = thresh_ts,
       exceeds = exceeds, events = events, anomaly = anomaly, events_df = events_df)
}

build_compound <- function(pipe_a, pipe_b, label = "") {
  common_time <- as.Date(intersect(as.character(pipe_a$time), as.character(pipe_b$time)))
  ia <- match(common_time, pipe_a$time); ib <- match(common_time, pipe_b$time)
  compound_raw <- ifelse(pipe_a$events[ia] == 1 & pipe_b$events[ib] == 1, 1, 0)
  compound <- remove_short_events(compound_raw, min_duration = 2)
  raw_intensity <- pipe_a$anomaly[ia] * pipe_b$anomaly[ib]
  compound_intensity <- (raw_intensity - mean(raw_intensity, na.rm = TRUE)) / sd(raw_intensity, na.rm = TRUE)
  cat("---", label, "---\n")
  events_df <- characterize_events(compound, compound_intensity, common_time)
  list(time = common_time, compound = compound, intensity = compound_intensity, events_df = events_df)
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 6: Run point-level analysis — GFDL, MOVING BASELINE ONLY FOR NOW ----------
# ══════════════════════════════════════════════════════════════════════════════

years_sst_gfdl <- 1995:2100
years_bt_gfdl  <- 1995:2045   # per earlier file-list flag; confirm actual coverage
years_bo2_gfdl <- 1995:2047

gfdl_sst <- load_yearly_point_sst(data_dir, years_sst_gfdl, sst_pt$i, sst_pt$j)
gfdl_bt  <- load_yearly_point_bottom(data_dir, "bT_atl_gfdl_monavg_%d.nc",  "temp",   years_bt_gfdl,  bottom_pt$i, bottom_pt$j)
gfdl_bo2 <- load_yearly_point_bottom(data_dir, "bO2_atl_gfdl_monavg_%d.nc", "oxygen", years_bo2_gfdl, bottom_pt$i, bottom_pt$j)

gfdl_smhw <- run_point_pipeline(gfdl_sst$time, gfdl_sst$ts, "high", "moving", label = "GFDL SMHW")
gfdl_bmhw <- run_point_pipeline(gfdl_bt$time,  gfdl_bt$ts,  "high", "moving", label = "GFDL BMHW")
gfdl_bhx  <- run_point_pipeline(gfdl_bo2$time, gfdl_bo2$ts, "low",  "moving", percentile = 0.10, label = "GFDL BHX")

# bT x bO2 compound: SAME grid
gfdl_compound_bT_bhx <- build_compound(gfdl_bmhw, gfdl_bhx, "GFDL BMHW x BHX")

# SST x bT or SST x bO2 compound: DIFFERENT grids, but at the POINT level this
# is fine because both pt lookups (sst_pt, bottom_pt) were independently found
# for the SAME real-world lat/lon — no regridding needed for single-point work,
# only for full-grid work
gfdl_compound_sst_bT  <- build_compound(gfdl_smhw, gfdl_bmhw, "GFDL SMHW x BMHW")
gfdl_compound_sst_bhx <- build_compound(gfdl_smhw, gfdl_bhx,  "GFDL SMHW x BHX")



# ══════════════════════════════════════════════════════════════════════════════
# PART 7: Full Grid ------------------------------------------------------------
# ══════════════════════════════════════════════════════════════════════════════








































