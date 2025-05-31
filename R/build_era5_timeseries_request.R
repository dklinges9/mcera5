#' Builds a request compatible with the CDS for ERA5 time-series climate
#' variables
#'
#' @description `build_era5_timeseries_request` creates a request or set of requests that
#' can be submitted to the Climate Data Store (CDS) with the `ecmwfr` package,
#' to query the (beta) ERA5 hourly time-series data for a single spatial point.
#' NOTE: currently such time-series data does not include all variables necessary
#' for microclimate modeling.
#' Spatial location and temporal duration are defined by the user, and requests are
#' automatically split by month. The following variables are requested:
#' 2m_temperature, 2m_dewpoint_temperature, surface_pressure,
#' 10m_u_component_of_wind, 10m_v_component_of_wind, total_precipitation.
#'
#' @param x The longitude to request data for.
#' @param ymin The latitude to request data for.
#' @param start_time a POSIXlt object indicating the first hour for which data
#' are required.
#' @param end_time a POSIXlt object indicating the last hour for which data
#' are required.
#' @param by_year logical indicating whether to split the request by year, which might be safe for especially long-duration queries which may be lost if connection with the API is unstable (default is TRUE)
#' @param filetype One of "netcdf" or "csv" (ONLY NETCDF FORMAT CURRENTLY OPERATIONAL)
#' @param outfile_name character prefix for .nc files when downloaded (the year, month, and file extension is automatically added). Keep in mind that you may want to submit multiple (or many) requests, and therefore
#' this prefix should adequately describe a unique query (e.g. by its
#' spatial or temporal extent).
#'
#' @export
#'
#'
build_era5_timeseries_request <- function(x, y, start_time, end_time,
                                          by_year = TRUE,
                                          filetype = "netcdf",
                                          outfile_name = "era5_timeseries_out") {
  # input checks
  if (missing(x)) {
    stop("xmin is missing")
  }
  if (missing(y)) {
    stop("ymin is missing")
  }
  if (missing(start_time)) {
    stop("start_time is missing")
  }
  if (missing(end_time)) {
    stop("end_time is missing")
  }

  # round to regular grid
  x_r <- plyr::round_any(x, .25, f = floor)
  y_r <- plyr::round_any(y, .25, f = floor)

  # month/year combos
  ut <- uni_dates(start_time, end_time)

  if (by_year) {
    # iterate over all focal years
    request <- apply(ut %>% dplyr::distinct(yea), 1, function(time) {

      # select months relevant to focal year
      sub_mon <- ut %>%
        dplyr::filter(., yea == time) %>%
        dplyr::select(., mon)

      # Find number of days in last month
      ndays <- days_in_month(max(sub_mon))

      list(
        "dataset_short_name" = "reanalysis-era5-single-levels-timeseries",
        "variable" = c(
          "2m_temperature", "2m_dewpoint_temperature", "surface_pressure",
          "10m_u_component_of_wind", "10m_v_component_of_wind",
          "total_precipitation"
        ),
        "location" = list("latitude" = y_r, "longitude" = x_r),
        "date" = paste0(time, "-", min(sub_mon), "-01", "/",
                        time, "-", max(sub_mon), "-", ndays),
        "format" = filetype,
        "target" = paste0(outfile_name, "_", as.character(time), ".zip")
      )
    })
  } else {
    # request all at once
    request <- list(
        "dataset_short_name" = "reanalysis-era5-single-levels-timeseries",
        "variable" = c(
          "2m_temperature", "2m_dewpoint_temperature", "surface_pressure",
          "10m_u_component_of_wind", "10m_v_component_of_wind",
          "total_precipitation"
        ),
        "location" = list("latitude" = y_r, "longitude" = x_r),
        "date" = paste0(as.Date(start_time), "/",  as.Date(end_time)),
        "format" = filetype,
        "target" = paste0(outfile_name, "_", as.character(time), ".zip")
      )
  }

  return(request)
}

