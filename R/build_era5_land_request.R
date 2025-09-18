#' Builds a request compatible with the CDS for all microclimate relevant climate
#' variables of ERA5-Land.
#'
#' @description `build_era5_land_request` creates a request or set of requests of
#' ERA5-Land data that can be submitted to the Climate Data Store (CDS) with
#' the `ecmwfr` package.
#' Spatial and temporal extents are defined by the user, and requests are
#' automatically split by month. The following variables are requested:
#' 2m_temperature, 2m_dewpoint_temperature, surface_pressure,
#' 10m_u_component_of_wind, 10m_v_component_of_wind, total_precipitation. Note that
#' these variables are not sufficient to run several microclimate models, e.g.
#' `microclimc` and `microclimf`.
#'
#' @param xmin The minimum longitude to request data for.
#' @param xmax The maximum longitude to request data for.
#' @param ymin The minimum latitude to request data for.
#' @param ymax The maximum latitude to request data for.
#' @param start_time a POSIXlt object indicating the first hour for which data
#' are required.
#' @param end_time a POSIXlt object indicating the last hour for which data
#' are required.
#' @param outfile_name character prefix for .nc files when downloaded (the year,
#' month, and file extension is automatically added). Keep in mind that you may
#' want to submit multiple (or many) requests, and therefore this prefix should
#' adequately describe a unique query (e.g. by its spatial or temporal extent).
#'
#' @export
#'
#'
build_era5_land_request <- function(xmin, xmax, ymin, ymax, start_time, end_time,
                                    outfile_name = "era5_land_out") {
  # input checks
  if (missing(xmin)) {
    stop("xmin is missing")
  }
  if (missing(xmax)) {
    stop("xmax is missing")
  }
  if (missing(ymin)) {
    stop("ymin is missing")
  }
  if (missing(ymax)) {
    stop("ymax is missing")
  }
  if (missing(start_time)) {
    stop("start_time is missing")
  }
  if (missing(end_time)) {
    stop("end_time is missing")
  }

  # round to regular grid
  xmin_r <- plyr::round_any(xmin, .1, f = floor)
  xmax_r <- plyr::round_any(xmax, .1, f = ceiling)
  ymin_r <- plyr::round_any(ymin, .1, f = floor)
  ymax_r <- plyr::round_any(ymax, .1, f = ceiling)
  # area of interest
  ar <- paste0(ymax_r, "/", xmin_r, "/", ymin_r, "/", xmax_r)
  # month/year combos
  ut <- uni_dates(start_time, end_time)

  # iterate over all focal months. ERA5-land must be queried on a monthly basis
  request <- apply(ut, 1, function(time) {


    mon <- time[1]
    yea <- time[2]
    # If first month....
    if (all(time == ut[1,])) {
      # Extract days from start_time to days_in_month()
      days <- format(start_time, "%d") : lubridate::days_in_month(paste0(yea, "-", mon, "-01"))
    }
    # If last month....(importantly, must happen after first month, to overwrite)
    if (all(time == ut[nrow(ut),])) {
      # Extract days from 1 to end_time
      days <- 1 : format(end_time, "%d")
    }

    if (any(time != ut[1,]) & any(time != ut[nrow(ut),])) {
      # Otherwise, get all days in month
      days <- 1 : lubridate::days_in_month(paste0(yea, "-", mon, "-01"))
    }

    # Coerce to character
    days <- as.character(days)
    # add "0" if only one character
    days <- ifelse(nchar(days) < 2, paste0("0", days), days)

      return(list(
        "dataset_short_name" = "reanalysis-era5-land",
        "product_type" = "reanalysis",
        "variable" = c(
          "2m_temperature", "2m_dewpoint_temperature", "surface_pressure",
          "10m_u_component_of_wind", "10m_v_component_of_wind", "total_precipitation",
          "surface_solar_radiation_downwards", "surface_net_thermal_radiation",
          "surface_thermal_radiation_downwards"
        ),
        "year" = as.character(time[2]),
        "month" = as.character(time[1]),
        "day" = days,
        "time" = c(
          "00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00",
          "07:00", "08:00", "09:00", "10:00", "11:00", "12:00", "13:00",
          "14:00", "15:00", "16:00", "17:00", "18:00", "19:00", "20:00",
          "21:00", "22:00", "23:00"
        ),
        "area" = ar,
        "format" = "netcdf",
        "download_format" = "unarchived",
        "target" = paste0(outfile_name, "_", time[2], "_", time[1], ".nc")
      ))
    })

  return(request)
}
