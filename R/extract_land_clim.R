#' Produces hourly data of ERA5-Land for a single location
#'
#' @description `extract_clim` takes an nc file containing hourly ERA5-Land climate
#' data, and for a given set of coordinates, produces an (optionally) inverse
#' distance weighted mean of each variable.
#'
#' @param nc character vector containing the path to the nc file. Use the
#' `build_era5_land_request` and `request_era5` functions to acquire an nc file with
#' the correct set of variables. Data within nc file must span the period
#' defined by start_time and end_time.
#' @param long longitude of the location for which data are required (decimal
#' degrees, -ve west of Greenwich Meridian).
#' @param lat latitude of the location for which data are required (decimal
#' degrees, -ve south of the equator).
#' @param start_time a POSIXlt or POSIXct object indicating the first day or hour
#' for which data are required. Encouraged to specify desired timezone as UTC (ERA5
#' data are in UTC by default), but any timezone is accepted.
#' @param end_time a POSIXlt or POSIXct object indicating the last day or hour for
#' which data are required. Encouraged to specify desired timezone as UTC (ERA5
#' data are in UTC by default), but any timezone is accepted.
#' @param d_weight logical value indicating whether to apply inverse distance
#' weighting using the 4 closest neighbouring points to the location defined by
#' `long` and `lat`. Default = `TRUE`.
#' @param cds_version specifies what version of the CDS (Climate Data Store) the
#' ERA5 data were queried from. Either "new" (queries made after Sept 2024) or
#' "legacy" (queries made before Sept 2024). This argument will be deprecated in
#' the future and only queries from the new CDS will be supported.
#'
#' @return a data frame containing hourly values for a suite of climate variables:
#' @return `obs_time` | the date-time (timezone specified in col timezone)
#' @return `temperature` | (degrees celsius)
#' @return `humidity` | specific humidity (kg / kg)
#' @return `pressure` | (Pa)
#' @return `windspeed` | (m / s)
#' @return `winddir` | wind direction, azimuth (degrees from north)
#' @return `szenith` | Solar zenith angle (degrees from a horizontal plane)
#' @return `timezone` | (unitless)
#'
#' @export
#'
#'

extract_land_clim <- function(nc, long, lat, start_time, end_time, d_weight = TRUE,
                              cds_version = "new") {
  # Open nc file for error trapping
  nc_dat <- ncdf4::nc_open(nc)

  ## Error trapping ---------------

  if (!cds_version %in% c("legacy", "new")) {
    stop("Argument `cds_version` must be one of the following values: `new` or`legacy`")
  }

  # Specify the base date-time, which differs between the CDS versions
  if (cds_version == "legacy") {
    base_datetime <- lubridate::ymd_hms("1900:01:01 00:00:00")
  } else {
    base_datetime <- lubridate::ymd_hms("1970:01:01 00:00:00")
  }

  # Confirm that start_time and end_time are date-time objects
  if (any(!class(start_time) %in% c("Date", "POSIXct", "POSIXt", "POSIXlt")) |
    any(!class(end_time) %in% c("Date", "POSIXct", "POSIXt", "POSIXlt"))) {
    stop("`start_time` and `end_time` must be provided as date-time objects.")
  }
  # Confirm that start_time and end_time are same class of date-time objects
  if (any(class(start_time) != class(end_time))) {
    stop("`start_time` and `end_time` must be of the same date-time class.")
  }

  # Check if start_time is after first time observation
  start <- lubridate::ymd_hms("1970:01:01 00:00:00") + (nc_dat$dim$valid_time$vals[1])
  if (start_time < start) {
    stop("Requested start time is before the beginning of time series of the ERA5 netCDF.")
  }

  # Check if end_time is before last time observation
  end <- lubridate::ymd_hms("1970:01:01 00:00:00") + (utils::tail(nc_dat$dim$valid_time$vals, n = 1))
  if (end_time > end) {
    stop("Requested end time is after the end of time series of the ERA5 netCDF.")
  }

  # Check if requested coordinates are in spatial grid
  if (long < min(nc_dat$dim$longitude$vals) | long > max(nc_dat$dim$longitude$vals)) {
    long_out <- TRUE
  } else {
    long_out <- FALSE
  }

  if (lat < min(nc_dat$dim$latitude$vals) | lat > max(nc_dat$dim$latitude$vals)) {
    lat_out <- TRUE
  } else {
    lat_out <- FALSE
  }

  # close nc file
  ncdf4::nc_close(nc_dat)

  if (long_out & lat_out) {
    stop("Requested coordinates are not represented in the ERA5 netCDF (both longitude and latitude out of range).")
  }
  if (long_out) {
    stop("Requested coordinates are not represented in the ERA5 netCDF (longitude out of range).")
  }
  if (lat_out) {
    stop("Requested coordinates are not represented in the ERA5 netCDF (latitude out of range).")
  }

  if (lubridate::tz(start_time) != lubridate::tz(end_time)) {
    stop("start_time and end_time are not in the same timezone.")
  }

  if (lubridate::tz(start_time) != "UTC" | lubridate::tz(end_time) != "UTC") {
    warning("provided times (start_time and end_time) are not in timezone UTC (default timezone of ERA5 data). Output will be provided in timezone UTC however.")
  }

  # Specify hour of end_time as last hour of day, if not specified
  if (lubridate::hour(end_time) == 0) {
    end_time <- as.POSIXlt(paste0(
      lubridate::year(end_time), "-",
      lubridate::month(end_time), "-",
      lubridate::day(end_time),
      " 23:00"
    ), tz = lubridate::tz(end_time))
  }

  if (sum((long %% .1) + (lat %% .1)) == 0 & d_weight == TRUE) {
    message("Input coordinates match ERA5 grid, no distance weighting required.")
    d_weight <- FALSE
  }

  # no distance weighting
  if (d_weight == FALSE) {
    long <- plyr::round_any(long, 0.1)
    lat <- plyr::round_any(lat, 0.1)
    dat <- nc_to_df_land(nc, long, lat, start_time, end_time)
    message("No distance weighting applied, nearest point used.")
  }
  # yes distance weighting
  if (d_weight == TRUE) {
    focal <- focal_dist(long, lat, margin = .1)
    # collector per weighted neighbour
    focal_collect <- list()
    for (j in 1:nrow(focal)) {
      # applies DTR correction if TRUE
      f_dat <- nc_to_df_land(nc, focal$x[j], focal$y[j], start_time, end_time) %>%
        dplyr::mutate(inverse_weight = focal$inverse_weight[j])
      focal_collect[[j]] <- f_dat
    }
    # create single weighted dataframe
    dat <- dplyr::bind_rows(focal_collect, .id = "neighbour") %>%
      dplyr::group_by(obs_time) %>%
      dplyr::summarise_at(
        dplyr::vars(
          temperature, humidity, pressure, windspeed, winddir, szenith
        ),
        weighted.mean,
        w = dplyr::quo(inverse_weight)
      ) %>%
      dplyr::mutate(timezone = lubridate::tz(obs_time))
    message("Distance weighting applied.")
  }
  return(dat)
}
