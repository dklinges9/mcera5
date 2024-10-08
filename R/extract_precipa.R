#' Produces daily or hourly precipitation data across a grid ready for use with
#' several gridded microclimate models.
#'
#' @description `extract_precipa` takes an nc file containing hourly ERA5
#' climate data, and for a given set of coordinates, produces precipitation
#' (at daily or hourly resolution) ready for use with `microclimf::modelina`.
#'
#' @param nc character vector containing the path to the nc file. Use the
#' `build_era5_request` and `request_era5` functions to acquire an nc file with
#' the correct set of variables. Data within nc file must span the period
#' defined by start_time and end_time.
#' @param long_min minimum longitude of the grid for which data are required (decimal
#' @param long_max maximum longitude of the grid for which data are required (decimal
#' degrees, -ve west of Greenwich Meridian).
#' @param lat_min minimum latitude of the grid for which data are required (decimal
#' @param lat_max maximum latitude of the grid for which data are required (decimal
#' degrees, -ve south of the equator).
#' @param start_time a POSIXlt or POSIXct object indicating the first day or hour
#' for which data are required. Encouraged to specify desired timezone as UTC (ERA5
#' data are in UTC by default), but any timezone is accepted.
#' @param end_time a POSIXlt or POSIXct object indicating the last day or hour for
#' which data are required. Encouraged to specify desired timezone as UTC (ERA5
#' data are in UTC by default), but any timezone is accepted.
#' @param convert_daily a flag indicating whether the user desires to convert the
#'  precipitation spatRaster from hourly to daily averages (TRUE) or remain as hourly
#' values (FALSE). Only daily precipitation will be accepted by `microclimf::modelina`.
#' @param cds_version specifies what version of the CDS (Climate Data Store) the
#' ERA5 data were queried from. Either "new" (queries made after Sept 2024) or
#' "legacy" (queries made before Sept 2024). This argument will be deprecated in
#' the future and only queries from the new CDS will be supported.
#'
#' @return a spatRaster of daily or hourly precipitation (mm).
#' @export
#'
#'
extract_precipa <- function(nc, long_min, long_max, lat_min, lat_max,
                            start_time, end_time, convert_daily = TRUE,
                            cds_version = "new") {

  # Open nc file for error trapping
  nc_dat = ncdf4::nc_open(nc)

  ## Error trapping ------------------

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

  if (long_max <= long_min) {
    stop("Maximum longitude must be greater than minimum longitude.")
  }

  if (lat_max <= lat_min) {
    stop("Maximum longitude must be greater than minimum longitude.")
  }

  if (abs(long_min) > 180 | abs(long_max) > 180 |
      abs(lat_min) > 90 | abs(lat_max) > 90) {
    stop("Coordinates must be provided in decimal degrees (longitude between -180 and 180, latitude between -90 and 90).")
  }

  # Check if requested coordinates are in spatial grid
  if (long_min < min(nc_dat$dim$longitude$vals) | long_min > max(nc_dat$dim$longitude$vals) |
      long_max < min(nc_dat$dim$longitude$vals) | long_max > max(nc_dat$dim$longitude$vals)
  ) {
    long_out <- TRUE
  } else {
    long_out <- FALSE
  }

  if (lat_min < min(nc_dat$dim$latitude$vals) | lat_min > max(nc_dat$dim$latitude$vals) |
      lat_max < min(nc_dat$dim$latitude$vals) | lat_max > max(nc_dat$dim$latitude$vals)) {
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

  if (lubridate::tz(start_time) != lubridate::tz(end_time)) {
    stop("start_time and end_time are not in the same timezone.")
  }

  if (lubridate::tz(start_time) != "UTC" | lubridate::tz(end_time) != "UTC") {
    warning("provided times (start_time and end_time) are not in timezone UTC (default timezone of ERA5 data). Output will be provided in timezone UTC however.")
  }

  # Specify hour of end_time as last hour of day, if not specified
  if (lubridate::hour(end_time) == 0) {
    end_time <- as.POSIXlt(paste0(lubridate::year(end_time), "-",
                                  lubridate::month(end_time), "-",
                                  lubridate::day(end_time),
                                  " 23:00"), tz = lubridate::tz(end_time))
  }

  tme <- as.POSIXct(seq(start_time,
                        end_time, by = 3600), tz = lubridate::tz(end_time))

  # Load in netCDF variable
  tp <- terra::rast(nc, subds = "tp")
  # Subset down to desired time period
  tp <- tp[[terra::time(tp) %in% tme]]
  # Subset down to desired spatial extent
  tp <- terra::crop(tp, terra::ext(long_min, long_max, lat_min, lat_max))

  # convert to daily
  if (convert_daily) {
    tp <- terra::tapp(tp, "days", fun = "sum")
  }
  # convert from m to mm of rain
  precip <- tp * 1000

  # Add timesteps back to names
  if (convert_daily) {
    names(precip) <- terra::time(precip)
  } else {
    names(precip) <- paste(terra::time(precip), lubridate::tz(terra::time(precip)))
  }

  return(precip)
}
