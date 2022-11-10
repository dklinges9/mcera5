#' Produces hourly data for a single location ready for use with
#' `microclima::runauto`.
#'
#' @description `extract_clim` takes an nc file containing hourly ERA5 climate
#' data, and for a given set of coordinates, produces an (optionally) inverse
#' distance weighted mean of each variable, ready for use with
#' `microclima::runauto`. Also provides the option to implement a diurnal
#' temperature range correction to air temperatures.
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
#' @param dtr_cor logical value indicating whether to apply a diurnal temperature
#' range correction to air temperature values. Default = `TRUE`.
#' @param dtr_cor_fac numeric value to be used in the diurnal temperature range
#' correction. Default = 1.
#'
#' @return a data frame containing hourly values for a suite of climate variables:
#' @return `obs_time` | the date-time (timezone specified in col timezone)
#' @return `temperature` | (degrees celsius)
#' @return `humidity` | specific humidity (kg / kg)
#' @return `pressure` | (Pa)
#' @return `windspeed` | (m / s)
#' @return `winddir` | wind direction, azimuth (degrees from north)
#' @return `emissivity` | downward long wave radiation flux divided by the sum
#' of net long-wave radiation flux and downward long wave radiation flux (unitless)
#' @return `cloudcover` | (percent)
#' @return `netlong` | Net longwave radiation (MJ / m2 / hr)
#' @return `uplong` | Upward longwave radiation (MJ / m2 / hr)
#' @return `downlong` | Downward longwave radiation (MJ / m2 / hr)
#' @return `rad_dni` | Direct normal irradiance (MJ / m2 / hr)
#' @return `rad_dif` | Diffuse normal irradiance (MJ / m2 / hr)
#' @return `szenith` | Solar zenith angle (degrees from a horizontal plane)
#' @return `timezone` | (unitless)
#'
#' @export
#'
#'

extract_clim <- function(nc, long, lat, start_time, end_time, d_weight = TRUE,
                           dtr_cor = TRUE, dtr_cor_fac = 1) {

  # Open nc file for error trapping
  nc_dat = ncdf4::nc_open(nc)

  ## Error trapping

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

  if(dtr_cor == TRUE & !is.numeric(dtr_cor_fac)) {
    stop("Invalid diurnal temperature range correction value provided.")}

  if(sum((long %% .25) + (lat %% .25)) == 0 & d_weight == TRUE) {
    message("Input coordinates match ERA5 grid, no distance weighting required.")
    d_weight = FALSE
    }

  # no distance weighting - dtr_cor passed to processing function
  if(d_weight == FALSE) {
    long <- plyr::round_any(long, 0.25)
    lat <- plyr::round_any(lat, 0.25)
    dat <- nc_to_df(nc, long, lat, start_time, end_time, dtr_cor = dtr_cor,
                    dtr_cor_fac = dtr_cor_fac)
    message("No distance weighting applied, nearest point used.")
    if(dtr_cor == TRUE) {
      message("Diurnal temperature range correction applied.")
    } else {
      message("No diurnal temperature range correction applied.")
    }
  }

    # yes distance weighting - dtr_cor passed to processing function
    if(d_weight == TRUE) {
      focal <- focal_dist(long, lat)
      # collector per weighted neighbour
      focal_collect <- list()
      for(j in 1:nrow(focal)) {
        # applies DTR correction if TRUE
        f_dat <- nc_to_df(nc, focal$x[j], focal$y[j], start_time, end_time,
                          dtr_cor = dtr_cor, dtr_cor_fac = dtr_cor_fac) %>%
          dplyr::mutate(inverse_weight = focal$inverse_weight[j])
        focal_collect[[j]] <- f_dat
      }
      # create single weighted dataframe
      dat <- dplyr::bind_rows(focal_collect, .id = "neighbour") %>%
        dplyr::group_by(obs_time)%>%
        dplyr::summarise_at(dplyr::vars(temperature, humidity, pressure, windspeed,
                                 winddir, emissivity, cloudcover, netlong,
                                 uplong, downlong, rad_dni, rad_dif, szenith),
                            weighted.mean, w = dplyr::quo(inverse_weight)) %>%
        dplyr::mutate(timezone = "UTC")
      message("Distance weighting applied.")
      if(dtr_cor == TRUE) {
        message("Diurnal temperature range correction applied.")
      } else {
        message("No diurnal temperature range correction applied.")
      }
    }
  return(dat)
}
