#' Builds a request compatible with the CDS for all microclimate relevant climate
#' variables.
#'
#' @description `build_era5_request` creates a request or set of requests that
#' can be submitted to the Climate Data Store (CDS) with the `ecmwfr` package.
#' Spatial and temporal extents are defined by the user, and requests are
#' automatically split by month. The following variables are requested:
#' 2m_temperature, 2m_dewpoint_temperature, surface_pressure,
#' 10m_u_component_of_wind, 10m_v_component_of_wind, total_precipitation,
#' total_cloud_cover, mean_surface_net_long_wave_radiation_flux,
#' mean_surface_downward_long_wave_radiation_flux,
#' total_sky_direct_solar_radiation_at_surface,
#' surface_solar_radiation_downwards.
#'
#' @param xmin The minimum longitude to request data for.
#' @param xmax The maximum longitude to request data for.
#' @param ymin The minimum latitude to request data for.
#' @param ymax The maximum latitude to request data for.
#' @param start_time a POSIXlt object indicating the first hour for which data
#' are required.
#' @param end_time a POSIXlt object indicating the last hour for which data
#' are required.
#' @param outfile_name character prefix for .nc files when downloaded (the year, month, and file extension is automatically added). Keep in
#' mind that you may want to submit multiple (or many) requests, and therefore
#' this prefix should adequately describe a unique query (e.g. by its
#' spatial or temporal extent).
#'
#' @export
#'
#'
build_era5_request <- function(xmin, xmax, ymin, ymax, start_time, end_time,
                               outfile_name = "era5_out") {
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
  xmin_r <- plyr::round_any(xmin, .25, f = floor)
  xmax_r <- plyr::round_any(xmax, .25, f = ceiling)
  ymin_r <- plyr::round_any(ymin, .25, f = floor)
  ymax_r <- plyr::round_any(ymax, .25, f = ceiling)
  # area of interest
  ar <- paste0(ymax_r, "/", xmin_r, "/", ymin_r, "/", xmax_r)
  # month/year combos
  ut <- uni_dates(start_time, end_time)

  # iterate over all focal months
  request <- apply(ut, 1, function(time) {
    list(
      "dataset_short_name" = "reanalysis-era5-single-levels",
      "product_type" = "reanalysis",
      "variable" = c(
        "2m_temperature", "2m_dewpoint_temperature", "surface_pressure",
        "10m_u_component_of_wind", "10m_v_component_of_wind",
        "total_precipitation", "total_cloud_cover",
        "mean_surface_net_long_wave_radiation_flux",
        "mean_surface_downward_long_wave_radiation_flux",
        "total_sky_direct_solar_radiation_at_surface",
        "surface_solar_radiation_downwards", "land_sea_mask"
      ),
      "year" = as.character(time[2]),
      "month" = as.character(time[1]),
      "day" = c(
        "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11",
        "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
        "23", "24", "25", "26", "27", "28", "29", "30", "31"
      ),
      "time" = c(
        "00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00",
        "07:00", "08:00", "09:00", "10:00", "11:00", "12:00", "13:00",
        "14:00", "15:00", "16:00", "17:00", "18:00", "19:00", "20:00",
        "21:00", "22:00", "23:00"
      ),
      "area" = ar,
      "format" = "netcdf",
      "target" = paste0(outfile_name, "_", time[2], "_", time[1], ".nc")
    )
  })

  return(request)
}

