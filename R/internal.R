#' function to calculate humidity from dew point temperature
#' @param tdew dewpoint temperature (°C)
#' @param tc air temperature (°C)
#' @param p pressure (Pa)
#' @return specific humidity (Kg/Kg)
#' @noRd
humfromdew <- function(tdew, tc, p) {
  pk <- p / 1000
  ea <- 0.6108 * exp(17.27 * tdew / (tdew + 237.3))
  e0 <- 0.6108 * exp(17.27 * tc / (tc + 237.3))
  s <- 0.622 * e0 / pk
  hr <- (ea / e0) * 100
  hs <- (hr / 100) * s
  hs
}

#' function to apply a correction to ERA5 temperatures based on proximity to the
#' coast
#' @param tc air temperature (°C)
#' @param tme a POSIXlt or POSIXct object corresponding to the timeseries of air temperatures
#' @param landprop single numeric value indicating proportion of grid cell that
#' is land (0 = all sea, 1 = all land)
#' @param cor_fac correction factor to be applied (default to 1.285 for UK)
#' @return air temperature (°C)
#' @noRd
coastal_correct <- function(tc, tme, landprop, cor_fac = 1.285) {
  tcdf <- data.frame(
    tc = tc,
    tme = tme
  ) %>%
    dplyr::mutate(
      yday = lubridate::yday(tme),
      year = lubridate::year(tme)
    )
  # Group by yday and year.....
  tcdf <- tcdf %>%
    dplyr::group_by(year, yday) %>%
    dplyr::summarize(tmean = mean(tc, na.rm = T), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::full_join(tcdf, by = dplyr::join_by(year, yday)) %>%
    dplyr::filter(stats::complete.cases(tme, tc))
  # ....to find daily temperature means
  tmean <- dplyr::pull(tcdf, tmean)
  m <- (1 - landprop) * cor_fac + 1
  # Subtract daily means, applying correction factor
  tdif <- (tc - tmean) * m
  # Add back daily means
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
  tme10min <- as.POSIXlt(
    seq(
      from = as.POSIXlt(tme1h[1] - 3600),
      to = as.POSIXlt(tme1h[length(tme1h)]),
      by = "10 min"
    ),
    origin = "1970-01-01 00:00",
    tz = "UTC"
  )
  csr10min <- clearskyrad(tme10min,
    lat = lat, long = long,
    merid = 0, dst = 0
  )
  # mean csr per hour (every 6 values)
  csr1h <- tibble::tibble(csr10min) %>%
    dplyr::group_by(group = gl(length(csr10min) / 6, 6)) %>%
    dplyr::summarize(csr1h = mean(csr10min)) %>%
    .$csr1h
  # watch out for NAs here - need modding?
  od <- bound(rad / csr1h)
  # if (is.na(od)[1]) od[1] <- mean(od, na.rm = T)
  # if (is.na(od)[length(od)]) od[length(od)] <- mean(od, na.rm = T)
  od[is.na(od)] <- 0
  # hourly time stamps on the half hour
  tme1h_2 <- as.POSIXlt(tme1h - 1800, tz = "UTC", origin = "1970-01-01 00:00")
  # length to interpolate to
  n <- length(tme1h_2) * 2
  # work out od values for every half hour
  od1h <- stats::spline(tme1h_2, od, n = n)$y
  # select just the ones on the hour
  od1h <- od1h[seq(2, length(od1h), 2)]
  # csr on the hour
  csr1h_2 <- clearskyrad(tme1h, lat = lat, long = long, merid = 0, dst = 0)
  # calculate corrected rad on the hour
  rad_out <- csr1h_2 * od1h
  rad_out[is.na(rad_out)] <- 0
  return(rad_out)
}

#' function to calculate the position of the 4 nearest neighbours to a point,
#' the xy distance and inverse weight of each.
#' @param long longitude
#' @param lat latitude
#' @param margin distance to the nearest neighbours. Default is 0.25° of ERA5 single levels
#' @return data frame of longitude, latitude, xy distance and inverse weight of
#' each neighbouring point
#' @noRd
focal_dist <- function(long, lat, margin = .25) {
  # round to nearest margin
  x_r <- plyr::round_any(long, margin)
  y_r <- plyr::round_any(lat, margin)
  # work out locations of the four neighbour points in the ERA5 dataset
  if (long >= x_r) {
    focal_x <- c(x_r, x_r, x_r + margin, x_r + margin)
  } else {
    focal_x <- c(x_r - margin, x_r - margin, x_r, x_r)
  }
  if (lat >= y_r) {
    focal_y <- c(y_r, y_r + margin, y_r, y_r + margin)
  } else {
    focal_y <- c(y_r - margin, y_r, y_r - margin, y_r)
  }
  # work out weighting based on dist between input & 4 surrounding points
  x_dist <- abs(long - focal_x)
  y_dist <- abs(lat - focal_y)
  xy_dist <- sqrt(x_dist^2 + y_dist^2)

  focal <- data.frame(x = focal_x, y = focal_y, xy_dist) %>%
    dplyr::mutate(., inverse_weight = 1 / sum(1 / (1 / sum(xy_dist) * xy_dist)) * 1 / (1 / sum(xy_dist) * xy_dist))
  return(focal)
}

#' function to  process relevant hourly climate data from an ERA5 nc to
#' a data frame
#' @param nc nc path to nc file requested with `build_era5_request` and downloaded with `request_era5`
#' @param long longitude
#' @param lat latitude
#' @param start_time start time for data required
#' @param end_time end time for data required
#' @param dtr_cor logical value indicating whether to apply a diurnal temperature
#' range correction to air temperature values. Default = `TRUE`.
#' @param dtr_cor_fac numeric value to be used in the diurnal temperature range
#' correction. Default = 1.
#' @return data frame of hourly climate variables
#' @noRd
nc_to_df <- function(nc, long, lat, start_time, end_time, dtr_cor = TRUE,
                     dtr_cor_fac = 1, reformat= NULL) {
  dat <- tidync::tidync(nc) %>%
    tidync::hyper_filter(
      longitude = longitude == long,
      latitude = latitude == lat
    ) %>%
    tidync::hyper_tibble() %>%
    dplyr::mutate(.,
      obs_time = lubridate::ymd_hms("1900:01:01 00:00:00") + (time * 3600),
      timezone = lubridate::tz(obs_time)
    ) %>% # convert to readable times
    dplyr::filter(., obs_time >= start_time & obs_time < end_time + 1) %>%
    dplyr::rename(., pressure = sp) %>%
    dplyr::mutate(.,
      temperature = t2m - 273.15, # kelvin to celcius
      lsm = dplyr::case_when(
        lsm < 0 ~ 0,
        lsm >= 0 ~ lsm
      ),
      temperature = dplyr::case_when(
        dtr_cor == TRUE ~ coastal_correct(temperature, obs_time, lsm, dtr_cor_fac),
        dtr_cor == FALSE ~ temperature
      ),
      humidity = humfromdew(d2m - 273.15, temperature, pressure),
      windspeed = sqrt(u10^2 + v10^2),
      windspeed = windheight(windspeed, 10, 2),
      winddir = (atan2(u10, v10) * 180 / pi + 180) %% 360,
      cloudcover = tcc * 100,
      netlong = abs(msnlwrf) * 0.0036,
      downlong = msdwlwrf * 0.0036,
      uplong = netlong + downlong,
      emissivity = downlong / uplong, # converted to MJ m-2 hr-1
      jd = julday(
        lubridate::year(obs_time),
        lubridate::month(obs_time),
        lubridate::day(obs_time)
      ),
      si = siflat(lubridate::hour(obs_time), lat, long, jd, merid = 0)
    ) %>%
    dplyr::mutate(.,
      rad_dni = fdir * 0.000001,
      rad_glbl = ssrd * 0.000001,
      rad_glbl = rad_calc(rad_glbl, obs_time, long, lat), # fix hourly rad
      rad_dni = rad_calc(rad_dni, obs_time, long, lat), # fix hourly rad
      rad_dif = rad_glbl - rad_dni * si
    ) %>% # converted to MJ m-2 hr-1 from J m-2 hr-1
    dplyr::mutate(., szenith = 90 - solalt(lubridate::hour(obs_time),
      lat, long, jd,
      merid = 0
    ))

  if(reformat == "micropoint") {
    dat = dat %>%
      dplyr::mutate(raddr = (rad_dni * si)/3600,
                    difrad = rad_dif/0.0036,
                    swdown = raddr + difrad,
                    pres = pressure/1000, # convert to kPa,
                    difrad = rad_dif/0.0036,
                    lwdown = downlong/0.0036) %>%
      rename(temp = temperature)
    dat$relhum = converthumidity(dat$humidity, intype = "specific",
                                 tc = dat$temp, pk = dat$pres)[["relative"]]
    dat$relhum = ifelse(dat$relhum > 100, 100, dat$relhum)
    dat = dat %>% dplyr::select(obs_time, temp, relhum, pres, swdown, difrad, lwdown, windspeed, winddir)
  } else {
    dat = dat %>%
      dplyr::select(
        ., obs_time, temperature, humidity, pressure, windspeed,
        winddir, emissivity, cloudcover, netlong, uplong, downlong,
        rad_dni, rad_dif, szenith, timezone
    )
  }


  return(dat)
}

#' function to  process relevant hourly climate data from an ERA5-Land nc to
#' a data frame
#' @param nc path to nc file requested with `build_era5_land_request` and downloaded with `request_era5`
#' @param long longitude
#' @param lat latitude
#' @param start_time start time for data required
#' @param end_time end time for data required
#' @return data frame of hourly climate variables
#' @noRd
nc_to_df_land <- function(nc, long, lat, start_time, end_time) {
  dat <- tidync::tidync(nc) %>%
    tidync::hyper_filter(
      # ERA5-Land does not return precise coordinate values
      longitude = longitude <= long + 1e-6 & longitude >= long - 1e-6,
      latitude = latitude <= lat + 1e-6 & latitude >= lat - 1e-6
    ) %>%
    tidync::hyper_tibble() %>%
    dplyr::mutate(.,
      obs_time = lubridate::ymd_hms("1900:01:01 00:00:00") + (time * 3600),
      timezone = lubridate::tz(obs_time)
    ) %>% # convert to readable times
    dplyr::filter(., obs_time >= start_time & obs_time < end_time + 1) %>%
    dplyr::rename(., pressure = sp) %>%
    dplyr::mutate(.,
      temperature = t2m - 273.15, # kelvin to celcius
      humidity = humfromdew(d2m - 273.15, temperature, pressure),
      windspeed = sqrt(u10^2 + v10^2),
      windspeed = windheight(windspeed, 10, 2),
      winddir = (atan2(u10, v10) * 180 / pi + 180) %% 360,
      jd = julday(
        lubridate::year(obs_time),
        lubridate::month(obs_time),
        lubridate::day(obs_time)
      ),
      si = siflat(lubridate::hour(obs_time), lat, long, jd, merid = 0)
    ) %>%
    dplyr::mutate(., szenith = 90 - solalt(lubridate::hour(obs_time),
      lat, long, jd,
      merid = 0
    )) %>%
    dplyr::select(
      ., obs_time, temperature, humidity, pressure, windspeed,
      winddir, szenith, timezone
    )

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
    tidync::hyper_filter(
      longitude = longitude == long,
      latitude = latitude == lat
    ) %>%
    tidync::hyper_tibble() %>%
    dplyr::mutate(.,
      obs_time = lubridate::ymd_hms("1900:01:01 00:00:00") + (time * 3600),
      timezone = lubridate::tz(obs_time)
    ) %>% # convert to readable times
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
  date_seq <- seq(lubridate::floor_date(as.Date(start_time), "month"), lubridate::ceiling_date(as.Date(end_time), "month") - lubridate::days(1), by = "month")
  df <- data.frame(
    mon = lubridate::month(date_seq),
    yea = lubridate::year(date_seq)
  )
  return(df)
}

#' Combines a series of netCDFs that all have the same spatial extent and
#' set of variables
#' @param filenames a list of filenames for netCDFs you wish to combine
#' @param combined_name the name of the combined netCDF
#' @noRd
combine_netcdf <- function(filenames, combined_name) {
  files <- lapply(filenames, function(x) {
    ncdf4::nc_open(x)
  })

  # Pull out first file for reference specs
  nc <- files[[1]]
  # Create an empty list to populate
  vars_list <- vector(mode = "list", length = nc$nvars)
  data_list <- vector(mode = "list", length = nc$nvars)
  # One variable at a time
  for (i in 1:length(names(nc$var))) {
    varname <- names(nc$var)[i]
    # Get the variable from each of the netCDFs
    vars_dat <- lapply(files, function(x) {
      ncdf4::ncvar_get(x, varname)
    })

    # Then bind all of the arrays together using abind, flexibly called via do.call
    data_list[[i]] <- do.call(abind::abind, list(
      ... = vars_dat,
      along = 3
    ))

    # To populate the time dimension, need to pull out the time values from each
    # netCDF
    timevals <- lapply(files, function(x) {
      x$dim$time$vals
    })

    # Create a netCDF variable
    vars_list[[i]] <- ncdf4::ncvar_def(
      name = varname,
      units = nc$var[varname][[varname]]$units,
      # Pull dimension names, units, and values from file1
      dim = list(
        # Longitude
        ncdf4::ncdim_def(
          nc$dim$longitude$name, nc$dim$longitude$units,
          nc$dim$longitude$vals
        ),
        # Latitude
        ncdf4::ncdim_def(
          nc$dim$latitude$name, nc$dim$latitude$units,
          nc$dim$latitude$vals
        ),
        # Time
        ncdf4::ncdim_def(
          nc$dim$time$name, nc$dim$time$units,
          # Combination of values of all files
          do.call(c, timevals)
        )
      )
    )
  }

  # Create a new file
  file_combined <- ncdf4::nc_create(
    # Filename from param combined_name
    filename = combined_name,
    # We need to define the variables here
    vars = vars_list
  )


  # And write to it (must write one variable at a time with ncdf4)
  for (i in 1:length(names(nc$var))) {
    ncdf4::ncvar_put(
      nc = file_combined,
      varid = names(nc$var)[i],
      vals = data_list[[i]]
    )
  }

  # Finally, close the file
  ncdf4::nc_close(file_combined)
}
