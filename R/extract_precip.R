#' Produces daily or hourly precipitation data for a single location ready for use
#' with `microclima::runauto`.
#'
#' @description `extract_precip` takes an nc file containing hourly ERA5
#' climate data, and for a given set of coordinates, produces an (optionally)
#' inverse distance weighted mean of precipitation (at daily or hourly resolution)
#' ready for use with `microclima::runauto`.
#'
#' @param nc character vector containing the path to the nc file. Use the
#' `build_era5_request` and `request_era5` functions to acquire an nc file with
#' the correct set of variables. Data within nc file must span the period
#' defined by start_time and end_time.
#' @param long longitude of the location for which data are required (decimal
#' degrees, -ve west of Greenwich Meridian).
#' @param lat latitude of the location for which data are required (decimal
#' degrees, -ve south of the equator).
#' @param start_time a POSIXlt object indicating the first hour for which data
#' are required.
#' @param end_time a POSIXlt object indicating the last hour for which data
#' are required.
#' @param d_weight logical value indicating whether to apply inverse distance
#' weighting using the 4 closest neighbouring points to the location defined by
#' `long` and `lat`. Default = `TRUE`.
#' @param convert_daily a flag indicating whether the user desires to convert the
#' precipitation vector from hourly to daily averages (TRUE) or remain as hourly
#' values (FALSE). Only daily precipitation will be accepted by `microclima::runauto`.
#'
#' @return a numeric vector of daily or hourly precipitation (mm).
#' @export
#'
#'
extract_precip <- function(nc, long, lat, start_time, end_time,
                                  d_weight = TRUE, convert_daily = TRUE) {

  # Open nc file for error trapping
  nc_dat = ncdf4::nc_open(nc)

  ## Error trapping

  # Check if start_time is after first time observation
  start <- lubridate::ymd_hms("1900:01:01 00:00:00") + (nc_dat$dim$time$vals[1] * 3600)
  if (start_time < start) {
    stop("Requested start time is before the beginning of time series of the ERA5 netCDF.")
  }

  # Check if end_time is before last time observation
  end <- lubridate::ymd_hms("1900:01:01 00:00:00") + (utils::tail(nc_dat$dim$time$vals, n = 1) * 3600)
  if (end_time > end) {
    stop("Requested end time is after the end of time series of the ERA5 netCDF.")
  }

  # Check if requested coordinates are in spatial grid
  if(long < min(nc_dat$dim$longitude$vals) | long > max(nc_dat$dim$longitude$vals)) {
    long_out <- TRUE
  } else {
    long_out <- FALSE
  }

  if(lat < min(nc_dat$dim$latitude$vals) | lat > max(nc_dat$dim$latitude$vals)) {
    lat_out <- TRUE
  } else {
    lat_out <- FALSE
  }

  # close nc file
  ncdf4::nc_close(nc_dat)

  if(long_out & lat_out) {
    stop("Requested coordinates are not represented in the ERA5 netCDF (both longitude and latitude out of range).")
  }
  if(long_out) {
    stop("Requested coordinates are not represented in the ERA5 netCDF (longitude out of range).")
  }
  if(lat_out) {
    stop("Requested coordinates are not represented in the ERA5 netCDF (latitude out of range).")
  }

  if(sum((long %% .25) + (lat %% .25)) == 0 & d_weight == TRUE) {
    message("Input coordinates match ERA5 grid, no distance weighting required.")
    d_weight = FALSE
  }

  if(d_weight == FALSE) {
    long <- plyr::round_any(long, 0.25)
    lat <- plyr::round_any(lat, 0.25)
    dat <- nc_to_df_precip(nc, long, lat, start_time, end_time) %>%
      dplyr::pull(precipitation)
    message("No distance weighting applied, nearest point used.")
  }

  # yes distance weighting - dtr_cor passed to processing function
  if(d_weight == TRUE) {
    focal <- focal_dist(long, lat)
    # collector per weighted neighbour
    focal_collect <- list()
    for(j in 1:nrow(focal)) {
      # applies DTR correction if TRUE
      f_dat <- nc_to_df_precip(nc, focal$x[j], focal$y[j], start_time, end_time) %>%
        dplyr::mutate(inverse_weight = focal$inverse_weight[j])
      focal_collect[[j]] <- f_dat
    }
    # create single weighted vector
    dat <- dplyr::bind_rows(focal_collect, .id = "neighbour") %>%
      dplyr::group_by(obs_time) %>%
      dplyr::summarise(precipitation = weighted.mean(precipitation,
                                                        w = inverse_weight)) %>%
      dplyr::pull(precipitation)
    message("Distance weighting applied.")
  }

  # convert to daily
  if (convert_daily) {
    precip <- matrix(dat,ncol=24,byrow=T) %>%
      rowSums() * 1000
  } else {
    precip <- dat
  }

  return(precip)
}
