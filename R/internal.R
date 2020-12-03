#' function to calculate humidity from dew point temperature
#' @param tdew dewpoint temperature (°C)
#' @param tc air temperature (°C)
#' @param p pressure (Pa)
#' @return specific humidity (Pa)
#' @noRd
humfromdew <- function(tdew, tc, p) {
  pk <- p/1000
  ea <- 0.6108 * exp(17.27 * tdew / (tdew + 237.3))
  e0 <- 0.6108 * exp(17.27 * tc/(tc + 237.3))
  s <- 0.622 * e0/pk
  hr <- (ea/e0) * 100
  hs <- (hr/100) * s
  hs
}

#' function to apply a correction to ERA5 temperatures based on proximity to the
#' coast
#' @param tc air temperature (°C)
#' @param landprop single numeric value indicating proportion of grid cell that
#' is land (0 = all sea, 1 = all land)
#' @param cor_fac correction factor to be applied (default to 1.285 for UK)
#' @return air temperature (°C)
#' @noRd
coastal_correct <- function(tc, landprop, cor_fac = 1.285) {
  td <- matrix(tc, ncol=24, byrow=T)
  tmean <- rep(apply(td, 1, mean), each=24)
  m <- (1 - landprop) * cor_fac + 1
  tdif <- (tc - tmean) * m
  tco <- tmean + tdif
  return(tco)
}


#' function to correct radiation for being an average over the hour rather than
#' on the hour
#' @param rad vector of radiation
#' @param tme vector of times
#' @param long longitude
#' @param lat latitude
#' @return vector of radiation
#' @noRd
rad_calc <- function(rad, tme, long, lat) {

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

#' function to calculate the position of the 4 nearest neighbours to a point,
#' the xy distance and inverse weight of each.
#' @param long longitude
#' @param lat latitude
#' @return data frame of longitude, latitude, xy distance and inverse weight of
#' each neighbouring point
#' @noRd
focal_dist <- function(long, lat) {

  # round to nearest 0.25
  x_r <- plyr::round_any(long, .25)
  y_r <- plyr::round_any(lat, .25)
  # work out locations of the four neighbour points in the ERA5 dataset
  if(long >= x_r) {
    focal_x <- c(x_r, x_r, x_r + 0.25, x_r + 0.25) } else {
      focal_x <- c(x_r - 0.25, x_r - 0.25, x_r, x_r)
    }
  if(lat >= y_r) {
    focal_y <- c(y_r, y_r + 0.25, y_r, y_r + 0.25) } else {
      focal_y <- c(y_r - 0.25, y_r, y_r - 0.25, y_r)
    }
  # work out weighting based on dist between input & 4 surrounding points
  x_dist <- abs(long - focal_x)
  y_dist <- abs(lat - focal_y)
  xy_dist <- sqrt(x_dist^2 + y_dist^2)

  focal <- data.frame(x = focal_x, y = focal_y, xy_dist) %>%
    dplyr::mutate(., inverse_weight = 1/sum(1/(1/sum(xy_dist) * xy_dist)) * 1/(1/sum(xy_dist) * xy_dist))
  return(focal)
}

#' function to  process relevant hourly climate data from an ERA5 nc to
#' a data frame
#' @param nc path to nc file downloaded with `request_era5`
#' @param long longitude
#' @param lat latitude
#' @param start_time start time for data required
#' @param end_time end time for data required
#' @return data frame of hourly climate variables
#' @noRd
nc_to_df <- function(nc, long, lat, start_time, end_time, dtr_cor = FALSE,
                     dtr_cor_fac = 1) {

  dat <- tidync::tidync(nc) %>%
    tidync::hyper_filter(longitude = longitude == long,
                         latitude = latitude == lat) %>%
    tidync::hyper_tibble() %>%
    dplyr::mutate(., obs_time = lubridate::ymd_hms("1900:01:01 00:00:00") + (time * 3600),
                  timezone = lubridate::tz(obs_time)) %>% # convert to readable times
    dplyr::filter(., obs_time >= start_time & obs_time < end_time + 1) %>%
    dplyr::rename(., pressure = sp) %>%
    dplyr::mutate(., temperature = t2m - 273.15, # kelvin to celcius
                  temperature = dplyr::case_when(
                    dtr_cor == TRUE ~ coastal_correct(temperature, lsm, dtr_cor_fac),
                    dtr_cor == FALSE ~ temperature),
                  humidity = humfromdew(d2m - 273.15, temperature, pressure),
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
                  si = microclima::siflat(lubridate::hour(obs_time), y, x, jd, merid = 0)) %>%
    dplyr::mutate(., rad_dni = fdir * 0.000001,
                  rad_glbl = ssrd * 0.000001,
                  rad_glbl = rad_calc(rad_glbl, obs_time, x, y), # fix hourly rad
                  rad_dni = rad_calc(rad_dni, obs_time, x ,y), # fix hourly rad
                  rad_dif = rad_glbl - rad_dni * si) %>% # converted to MJ m-2 hr-1 from J m-2 hr-1
    dplyr::mutate(., szenith = 90 - microclima::solalt(lubridate::hour(obs_time),
                                                       y, x, jd, merid = 0)) %>%
    dplyr::select(.,obs_time, temperature, humidity, pressure, windspeed,
                  winddir, emissivity, cloudcover, netlong, uplong, downlong,
                  rad_dni, rad_dif, szenith, timezone)

  return(dat)
}

# process relevant precipitation data from an ERA5 nc to data frame
#' a function to  process relevant precipitation data from an ERA5 nc to
#' a data frame
#' @param nc path to nc file downloaded with `request_era5`
#' @param long longitude
#' @param lat latitude
#' @param start_time start time for data required
#' @param end_time end time for data required
#' @return vector of daily precipitation values
#' @noRd
nc_to_df_precip <- function(nc, long, lat, start_time, end_time) {

  dat <- tidync::tidync(nc) %>%
    tidync::hyper_filter(longitude = longitude == long,
                         latitude = latitude == lat) %>%
    tidync::hyper_tibble() %>%
    dplyr::mutate(., obs_time = lubridate::ymd_hms("1900:01:01 00:00:00") + (time * 3600),
                  timezone = lubridate::tz(obs_time)) %>% # convert to readable times
    dplyr::filter(., obs_time >= start_time & obs_time < end_time + 1) %>%
    dplyr::rename(., precipitation = tp) %>%
    dplyr::select(., obs_time, precipitation)

  return(dat)
}

#' creates a data frame of unique month and year pairs from input
#' start/end times
#' @param start_time start time
#' @param end_time end time
#' @return data frame of unique months and years
#' @noRd
uni_dates <- function(start_time, end_time) {

  tme <- seq(as.POSIXlt(start_time), as.POSIXlt(end_time), by = 1)
  mon <- lubridate::month(tme)
  yea <- lubridate::year(tme)
  df <- data.frame(mon,yea) %>%
    dplyr::distinct(.)

  return(df)
}
