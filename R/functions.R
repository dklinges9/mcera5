# calculate humidity from dew point temperature
# tdew = Dewpoint temperature (deg C)
# tc = temperature (deg C)
# p = pressure (pascals)
.humfromdew <- function(tdew, tc, p) {
  pk <- p/1000
  ea <- 0.6108 * exp(17.27 * tdew / (tdew + 237.3))
  e0 <- 0.6108 * exp(17.27 * tc/(tc + 237.3))
  s <- 0.622 * e0/pk
  hr <- (ea/e0) * 100
  hs <- (hr/100) * s
  hs # specific humidity
}

# correct radiation for being an average over the hour rather than on the hour
.rad_calc <- function(rad, tme, long, lat) {

  bound <- function(x, mn = 0, mx = 1) {
    x[x > mx] <- mx
    x[x < mn] <- mn
    x
  }

  tme1h <- as.POSIXlt(tme, tz = "UTC", origin = "1970-01-01 00:00")
  # interpolate times to 10 mins
  tme10min <- as.POSIXlt(seq(from = as.POSIXlt(tme1h[1]-3600),
                             to = as.POSIXlt(tme1h[length(tme1h)]),
                             by = "10 min"), origin = "1970-01-01 00:00",
                         tz = "UTC")
  csr10min <- microclima::clearskyrad(tme10min, lat = lat , long = long,
                                      merid = 0, dst = 0)

  # mean csr per hour (every 6 values)
  csr1h <- tibble::tibble(csr10min) %>%
    dplyr::group_by(group = gl(length(csr10min)/6, 6)) %>%
    dplyr::summarise(csr1h = mean(csr10min)) %>%
    .$csr1h

  # watch out for NAs here - need modding?
  od <- bound(rad / csr1h)
  #if (is.na(od)[1]) od[1] <- mean(od, na.rm = T)
  #if (is.na(od)[length(od)]) od[length(od)] <- mean(od, na.rm = T)
  od[is.na(od)] <- 0

  # hourly time stamps on the half hour
  tme1h_2 <- as.POSIXlt(tme1h - 1800, tz = "UTC", origin = "1970-01-01 00:00")
  # length to interpolate to
  n <- length(tme1h_2) * 2
  # work out od values for every half hour
  od1h <- spline(tme1h_2, od, n = n)$y
  # select just the ones on the hour
  od1h <- od1h[seq(2,length(od1h),2)]
  # csr on the hour
  csr1h_2 <- microclima::clearskyrad(tme1h, lat = lat, long = long, merid = 0, dst = 0)
  # calculate corrected rad on the hour
  rad_out <- csr1h_2 * od1h
  rad_out[is.na(rad_out)] <- 0

  return(rad_out)
}

# calculate position of 4 nearest neighbours to a point, and xy distance to each.
.focal_dist <- function(x, y) {

  # round to nearest 0.25
  x_r <- plyr::round_any(x, .25)
  y_r <- plyr::round_any(y, .25)

  # work out locations of the four neighbour points in the ERA5 dataset
  if(x >= x_r) {
    focal_x <- c(x_r, x_r, x_r + 0.25, x_r + 0.25) } else {
      focal_x <- c(x_r - 0.25, x_r - 0.25, x_r, x_r)
    }
  if(y >= y_r) {
    focal_y <- c(y_r, y_r + 0.25, y_r, y_r + 0.25) } else {
      focal_y <- c(y_r - 0.25, y_r, y_r - 0.25, y_r)
    }

  # work out weighting based on dist between input & 4 surrounding points
  x_dist <- abs(x - focal_x)
  y_dist <- abs(y - focal_y)
  xy_dist <- sqrt(x_dist^2 + y_dist^2)

  focal <- data.frame(x = focal_x, y = focal_y, xy_dist)

  return(focal)
}

# process relevant hourly climate data from an ERA5 nc to data frame
.nc_to_df <- function(nc, x, y, start_time, end_time) {

  dat <- tidync::tidync(nc) %>%
    tidync::hyper_filter(longitude = longitude == x,
                         latitude = latitude == y) %>%
    tidync::hyper_tibble() %>%
    dplyr::mutate(., obs_time = lubridate::ymd_hms("1900:01:01 00:00:00") + (time * 3600),
                  timezone = lubridate::tz(obs_time)) %>% # convert to readable times
    dplyr::filter(., obs_time >= start_time & obs_time < end_time + 1) %>%
    dplyr::rename(., pressure = sp) %>%
    dplyr::mutate(., temperature = t2m - 273.15, # kelvin to celcius
                  humidity = .humfromdew(d2m - 273.15, temperature, pressure),
                  windspeed = sqrt(u10^2 + v10^2),
                  windspeed = microclima::windheight(windspeed, 10, 2),
                  winddir = (atan2(u10, v10) * 180/pi + 180)%%360,
                  cloudcover = tcc * 100,
                  netlong = abs(msnlwrf) * 0.0036,
                  downlong = msdwlwrf * 0.0036,
                  uplong = netlong + downlong,
                  emissivity = downlong/uplong, # converted to MJ m-2 hr-1
                  jd = microclima::julday(lubridate::year(obs_time),
                                          lubridate::month(obs_time),
                                          lubridate::day(obs_time)),
                  si = microclima::siflat(lubridate::hour(obs_time), y, x, jd)) %>%
    dplyr::mutate(., rad_dni = fdir * 0.000001,
                  rad_glbl = ssrd * 0.000001,
                  rad_glbl = .rad_calc(rad_glbl, obs_time, x, y), # fix hourly rad
                  rad_dni = .rad_calc(rad_dni, obs_time, x ,y), # fix hourly rad
                  rad_dif = rad_glbl - rad_dni * si) %>% # converted to MJ m-2 hr-1 from J m-2 hr-1
    dplyr::mutate(., szenith = 90 - microclima::solalt(lubridate::hour(obs_time),
                                                       y, x, jd, merid = 0)) %>%
    dplyr::select(.,obs_time, temperature, humidity, pressure, windspeed,
                  winddir, emissivity, cloudcover, netlong, uplong, downlong,
                  rad_dni, rad_dif, szenith, timezone)

  return(dat)
}

# process relevant precipitation data from an ERA5 nc to data frame
.nc_to_df_precip <- function(nc, x, y, start_time, end_time) {

  dat <- tidync::tidync(nc) %>%
    tidync::hyper_filter(longitude = longitude == x,
                         latitude = latitude == y) %>%
    tidync::hyper_tibble() %>%
    dplyr::mutate(., obs_time = lubridate::ymd_hms("1900:01:01 00:00:00") + (time * 3600),
                  timezone = lubridate::tz(obs_time)) %>% # convert to readable times
    dplyr::filter(., obs_time >= start_time & obs_time < end_time + 1) %>%
    dplyr::rename(., precipitation = tp) %>%
    dplyr::select(., obs_time, precipitation)

  return(dat)
}

# creates a data frame of unique month and year pairs from input start/end times
.uni_dates <- function(start_time, end_time) {

  tme <- as.POSIXlt(seq(start_time,end_time, by = 1))
  mon <- lubridate::month(tme)
  yea <- lubridate::year(tme)
  df <- data.frame(mon,yea) %>%
    dplyr::distinct(.)

  return(df)
}

#' Builds a request compatible with the CDS for all microclimate relevant climate
#' variables.
#'
#' @description `build_era5_request` creates a request or set of requests that
#' can be submitted to the Climate Data Store (CDS) with the `ecmwfr` package.
#' Spatial and temporal extents are defined by the user, and requests are
#' automatically split by year. The following variables requested:
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
#' @param outfile_name character prefix for .nc files when downloaded.
#'
build_era5_request <- function(xmin, xmax, ymin, ymax, start_time, end_time,
                               outfile_name = "era5_out") {

  # input checks
  if(missing(xmin)) { stop("xmin is missing")}
  if(missing(xmax)) { stop("xmax is missing")}
  if(missing(ymin)) { stop("ymin is missing")}
  if(missing(ymax)) { stop("ymax is missing")}
  if(missing(start_time)) { stop("start_time is missing")}
  if(missing(end_time)) { stop("end_time is missing")}

  # round to regular grid
  xmin_r <- plyr::round_any(xmin, .25, f = floor)
  xmax_r <- plyr::round_any(xmax, .25, f = ceiling)
  ymin_r <- plyr::round_any(ymin, .25, f = floor)
  ymax_r <- plyr::round_any(ymax, .25, f = ceiling)
  # area of interest
  ar <- paste0(ymax_r,"/",xmin_r,"/",ymin_r,"/",xmax_r)
  # month/year combos
  ut <- .uni_dates(start_time, end_time)
  request <- list()

  # loop through focal years
  for(i in 1:length(unique(ut$yea))) {

    # focal year
    yr <- unique(ut$yea)[i]
    # select months relevant to focal year
    sub_mon <- ut %>%
      dplyr::filter(., yea == yr) %>%
      dplyr::select(., mon)

    sub_request <- list(
      area = ar,
      product_type = "reanalysis",
      format = "netcdf",
      variable = c("2m_temperature", "2m_dewpoint_temperature", "surface_pressure",
                   "10m_u_component_of_wind", "10m_v_component_of_wind",
                   "total_precipitation", "total_cloud_cover",
                   "mean_surface_net_long_wave_radiation_flux",
                   "mean_surface_downward_long_wave_radiation_flux",
                   "total_sky_direct_solar_radiation_at_surface",
                   "surface_solar_radiation_downwards", "land_sea_mask"),
      year = as.character(yr),
      month = as.character(sub_mon$mon),
      day = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11",
              "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
              "23", "24", "25", "26", "27", "28", "29", "30", "31"),
      time = c("00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00",
               "07:00", "08:00", "09:00", "10:00", "11:00", "12:00", "13:00",
               "14:00", "15:00", "16:00", "17:00", "18:00", "19:00", "20:00",
               "21:00", "22:00", "23:00"),
      dataset = "reanalysis-era5-single-levels",
      target = paste0(outfile_name,"_",yr,".nc"))

    request[[i]] <- sub_request
  }
  return(request)
}

#' Submit request(s) to CDS for ERA5 data download
#'
#' @description Uses the `ecwmfr` package to submit requests to the Climate
#' Data Store (CDS) for ERA5 .nc files. Executes one request at a time.
#'
#' @param request a list of request(s) created with `build_era5_request`.
#' @param uid character vector containing your CDS user ID.
#' @out_path character vector of the location at which to download nc files.
#'
request_era5 <- function(request, uid, out_path) {

  for(req in 1:length(request)) {

    ecmwfr::wf_request(user = as.character(uid),
                       request = request[[req]],
                       transfer = TRUE,
                       path = out_path,
                       verbose = TRUE,
                       time_out = 18000)
  }
}

#' Produces spatially weighted hourly data for a single location ready for use
#' with microclima::runauto
#'
#' @description `point_nc_to_df` takes an nc file or set of nc files containing
#' hourly ERA5 climate data, and for a given set of coordinates, produces an
#' inverse distance weighted mean of each variable, ready for use with
#' microclima::runauto.
#'
#' @param nc character vector containing nc filename(s). Data within nc files
#' must span the period defined by start_time and end_time.
#' @param x longitude of the location for which data are required (decimal
#' degrees, -ve west of Greenwich Meridian).
#' @param y latitude of the location for which data are required (decimal
#' degrees, -ve south of the equator).
#' @param start_time a POSIXlt object indicating the first hour for which data
#' are required.
#' @param end_time a POSIXlt object indicating the last hour for which data
#' are required.
#' @param lsm default TRUE. Use land sea mask filtering?
#' @param lsm_thresh threshold to divide sea (0) and land (1). Default = 0.05.
#'
#' @return a data frame containing hourly values for a suite of climate
#' variables.
#' @export
#'
#' @examples
#'
point_nc_to_df <- function(nc, x, y, start_time, end_time, lsm = TRUE,
                           lsm_thresh = 0.05) {

  # if the desired coords fall right on an ERA5 point - weighting calc not needed
  if(sum((x %% .25) + (y %% .25)) != 0) {

    # if using lsm
    if(lsm == TRUE) {

      # pull lsm_vals for x/y and surrounding points from nc file + apply threshold
      lsm_vals <- tidync::tidync(nc) %>%
        tidync::hyper_filter(longitude = longitude >= x - 0.5 & longitude <= x + 0.5,
                             latitude = latitude >= y -0.5 & latitude <= y + 0.5,
                             lsm > lsm_thresh) %>%
        hyper_tibble() %>%
        dplyr::select(., x = longitude, y = latitude, lsm)

      focal <- .focal_dist(x, y) %>%
        dplyr::mutate(., inverse_weight = 1/sum(1/(1/sum(xy_dist) * xy_dist)) * 1/(1/sum(xy_dist) * xy_dist)) %>%
        dplyr::inner_join(., lsm_vals, by = c("x","y"))

      if(nrow(focal) == 0){
        stop(strwrap(paste0("No land pixels to use in weighting. Provide a different location or reduce the land/sea mask threshold.", collapse = "\n")))
      } else if(nrow(focal) < 4) {
        warning(strwrap(paste0(4-nrow(focal), " neighbours are sea pixels and have been omitted. Optionally provide a different location or reduce the land/sea mask threshold to increase the number of neighbours.")), collapse = "\n")
      }

    } else {

      focal <- .focal_dist(x, y) %>%
        dplyr::mutate(., inverse_weight = 1/sum(1/(1/sum(xy_dist) * xy_dist)) * 1/(1/sum(xy_dist) * xy_dist))
    }

    # collector per nc file
    nc_collect <- list()

    for(i in 1:length(nc)) {

      # collector per weighted neighbour
      focal_collect <- list()

      # gather and manipulate data for each focal point
      for(j in 1:nrow(focal)) {

        dat <- .nc_to_df(nc[i], focal$x[j], focal$y[j], start_time, end_time) %>%
          dplyr::mutate(., inverse_weight = focal$inverse_weight[j])

        # add to collector
        focal_collect[[j]] <- dat
      }

      w_dat <- bind_rows(focal_collect, .id = "neighbour") %>%
        dplyr::group_by(., obs_time)%>%
        dplyr::summarise_at(., vars(temperature, humidity, pressure, windspeed,
                                    winddir, emissivity, cloudcover, netlong,
                                    uplong, downlong, rad_dni, rad_dif, szenith),
                            weighted.mean, w = quo(inverse_weight)) %>%
        dplyr::mutate(., timezone = "UTC")

      nc_collect[[i]] <- w_dat
    }

    # skip weighting calc - use original coords
  } else {

    message("Coordinates match ERA5 grid, no weighting calculation required.")

    # collector per nc file
    nc_collect <- list()

    for(i in 1:length(nc)) {

      dat <- .nc_to_df(nc[i], x, y, start_time, end_time)
      nc_collect[[i]] <- dat
    }
  }

  out <- dplyr::bind_rows(nc_collect)
  return(out)
}

#' Produces spatially weighted hourly data for a single location ready for use
#' with microclima::runauto
#'
#' @description `point_nc_to_df_precip` takes an nc file or set of nc files
#' containing hourly ERA5 climate data, and for a given set of coordinates,
#' produces an inverse distance weighted mean of daily precipitation, ready for
#' use with microclima::runauto.
#'
#' @param nc character vector containing nc filename(s). Data within nc files
#' must span the period defined by start_time and end_time.
#' @param x longitude of the location for which data are required (decimal
#' degrees, -ve west of Greenwich Meridian).
#' @param y latitude of the location for which data are required (decimal
#' degrees, -ve south of the equator).
#' @param start_time a POSIXlt object indicating the first hour for which data
#' are required.
#' @param end_time a POSIXlt object indicating the last hour for which data
#' are required.
#' @param lsm_nc character vector containing the ERA5 land/sea mask nc filename.
#' @param lsm_thresh threshold to divide sea (0) and land (1). Default = 0.05.
#' #'
#' @return a numeric vector of daily precipitation.
#' @export
#'
#' @examples
#'
point_nc_to_df_precip <- function(nc, x, y, start_time, end_time, lsm_nc = NULL,
                                  lsm_thresh = 0.05) {

  if(!is.null(lsm_nc)) {

    # load and correct land sea mask / apply threshold
    lsm <- tidync::tidync(lsm_nc) %>%
      tidync::hyper_tibble() %>%
      dplyr::filter(., lsm > lsm_thresh) %>%
      dplyr::select(x = longitude, y = latitude, lsm) %>%
      dplyr::mutate(., x = dplyr::case_when(
        x > 180 ~ (x - 360),
        TRUE ~ as.numeric(x)))
  }

  # if the desired coords fall right on an ERA5 point - weighting calc not needed
  if(sum((x %% .25) + (y %% .25)) != 0) {

    # if using lsm
    if(!is.null(lsm_nc)) {
      focal <- .focal_dist(x, y) %>%
        dplyr::mutate(., inverse_weight = 1/sum(1/(1/sum(xy_dist) * xy_dist)) * 1/(1/sum(xy_dist) * xy_dist)) %>%
        dplyr::inner_join(., lsm, by = c("x","y"))

      if(nrow(focal) == 0){
        stop(strwrap(paste0("No land pixels to use in weighting. Provide a different location or reduce the land/sea mask threshold.", collapse = "\n")))
      } else if(nrow(focal) < 4) {
        warning(strwrap(paste0(4-nrow(focal), " neighbours are sea pixels and have been omitted. Optionally provide a different location or reduce the land/sea mask threshold to increase the number of neighbours.")), collapse = "\n")
      }

    } else {

      focal <- .focal_dist(x, y) %>%
        dplyr::mutate(., inverse_weight = 1/sum(1/(1/sum(xy_dist) * xy_dist)) * 1/(1/sum(xy_dist) * xy_dist))
    }

    # collector per nc file
    nc_collect <- list()

    for(i in 1:length(nc)) {

      # collector per weighted neighbour
      focal_collect <- list()

      # gather and manipulate data for each focal point
      for(j in 1:nrow(focal)) {

        dat <- .nc_to_df_precip(nc[i], focal$x[j], focal$y[j], start_time, end_time) %>%
          dplyr::mutate(., inverse_weight = focal$inverse_weight[j])

        # add to collector
        focal_collect[[j]] <- dat
      }

      w_dat <- bind_rows(focal_collect, .id = "neighbour") %>%
        dplyr::group_by(., obs_time) %>%
        dplyr::summarise(., precipitation = weighted.mean(precipitation,
                                                          w = inverse_weight))

      nc_collect[[i]] <- w_dat
    }

    # skip weighting calc - use original coords
  } else {

    message("Coordinates match ERA5 grid, no weighting calculation required.")

    # collector per nc file
    nc_collect <- list()

    for(i in 1:length(nc)) {

      dat <- .nc_to_df_precip(nc[i], x, y, start_time, end_time)
      nc_collect[[i]] <- dat
    }
  }

  out <- dplyr::bind_rows(nc_collect) %>%
    dplyr::mutate(., jd = microclima::julday(lubridate::year(obs_time),
                                             lubridate::month(obs_time),
                                             lubridate::day(obs_time))) %>%
    dplyr::group_by(., jd) %>%
    dplyr::summarise(., daily_precip = sum(precipitation)) %>%
    .$daily_precip

  return(out)
}
