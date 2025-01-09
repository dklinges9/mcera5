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
                     dtr_cor_fac = 1, format = "microclima") {

  # Extract time dimension
  timedim <- extract_timedim(ncdf4::nc_open(nc))
  # Find basetime from units
  base_datetime <- as.POSIXct(gsub(".*since ", "", timedim$units), tz = "UTC")
  # Extract time values
  nc_datetimes <- timedim$vals
  # If units in hours, multiply by 3600 to convert to seconds
  nc_datetimes <- nc_datetimes * ifelse(
    grepl("hours", timedim$units), 3600, 1
  )

  dat <- tidync::tidync(nc) %>%
    # hyper_filter flexibly, querying the nearest long and lat coords, rather
    # than explicitly filtering to coordinates -- avoids errors in hyper_filter()
    tidync::hyper_filter(
      longitude = longitude == longitude[which.min(abs(longitude - long))],
      latitude = latitude == latitude[which.min(abs(latitude - lat))]
    ) %>%
    tidync::hyper_tibble() %>%
    dplyr::mutate(.,
      obs_time = c(base_datetime + nc_datetimes),
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
      netlong = abs(msnlwrf) * 0.0036,  # converted to MJ m-2 hr-1
      downlong = msdwlwrf * 0.0036,  # converted to MJ m-2 hr-1
      uplong = netlong + downlong,
      emissivity = downlong / uplong,
      precip = tp * 1000,
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

  if (format %in% c("microclima", "NicheMapR")) {
    dat = dat %>%
      dplyr::select(
        ., obs_time, temperature, humidity, pressure, windspeed,
        winddir, emissivity, netlong, uplong, downlong,
        rad_dni, rad_dif, szenith, cloudcover, timezone
      )
  }

  if(format %in% c("microclimc")) {
    dat = dat %>%
      dplyr::mutate(raddr = (rad_dni * si)/0.0036,
                    difrad = rad_dif/0.0036,
                    swrad = raddr + difrad,
                    pres = pressure/1000, # convert to kPa,
                    difrad = rad_dif/0.0036,
                    precip = precip,
                    lwdown = downlong/0.0036) %>%
      dplyr::rename(temp = temperature)
    dat$relhum = converthumidity(dat$humidity, intype = "specific",
                                 tc = dat$temp, pk = dat$pres)[["relative"]]
    dat$relhum = ifelse(dat$relhum > 100, 100, dat$relhum)
    dat = dat %>%
      dplyr::rename(skyem = emissivity) %>%
      dplyr::select(obs_time, temp, relhum, pres, swrad, difrad, skyem,
                                windspeed, winddir, precip)
  }

  if (format %in% c("micropoint", "microclimf")) {
    dat = dat %>%
      dplyr::mutate(raddr = (rad_dni * si)/0.0036,
                    difrad = rad_dif/0.0036,
                    swdown = raddr + difrad,
                    pres = pressure/1000, # convert to kPa,
                    difrad = rad_dif/0.0036,
                    precip = precip,
                    lwdown = downlong/0.0036) %>%
      dplyr::rename(temp = temperature)
    dat$relhum = converthumidity(dat$humidity, intype = "specific",
                                 tc = dat$temp, pk = dat$pres)[["relative"]]
    dat$relhum = ifelse(dat$relhum > 100, 100, dat$relhum)
    dat = dat %>%
      dplyr::select(obs_time, temp, relhum, pres, swdown, difrad, lwdown, windspeed, winddir, precip)
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

  # Extract time dimension
  timedim <- extract_timedim(ncdf4::nc_open(nc))
  # Find basetime from units
  base_datetime <- as.POSIXct(gsub(".*since ", "", timedim$units), tz = "UTC")
  # Extract time values
  nc_datetimes <- timedim$vals
  # If units in hours, multiply by 3600 to convert to seconds
  nc_datetimes <- nc_datetimes * ifelse(
    grepl("hours", timedim$units), 3600, 1
  )

  dat <- tidync::tidync(nc) %>%
    # hyper_filter flexibly, querying the nearest long and lat coords, rather
    # than explicitly filtering to coordinates -- avoids errors in hyper_filter()
    tidync::hyper_filter(
      longitude = longitude == longitude[which.min(abs(longitude - long))],
      latitude = latitude == latitude[which.min(abs(latitude - lat))]
    ) %>%
    tidync::hyper_tibble() %>%
    dplyr::mutate(.,
      obs_time = base_datetime + nc_datetimes,
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

  # Extract time dimension
  timedim <- extract_timedim(ncdf4::nc_open(nc))
  # Find basetime from units
  base_datetime <- as.POSIXct(gsub(".*since ", "", timedim$units), tz = "UTC")
  # Extract time values
  nc_datetimes <- timedim$vals
  # If units in hours, multiply by 3600 to convert to seconds
  nc_datetimes <- nc_datetimes * ifelse(
    grepl("hours", timedim$units), 3600, 1
  )

  dat <- tidync::tidync(nc) %>%
    # hyper_filter flexibly, querying the nearest long and lat coords, rather
    # than explicitly filtering to coordinates -- avoids errors in hyper_filter()
    tidync::hyper_filter(
      longitude = longitude == longitude[which.min(abs(longitude - long))],
      latitude = latitude == latitude[which.min(abs(latitude - lat))]
    ) %>%
    tidync::hyper_tibble() %>%
    dplyr::mutate(.,
      obs_time = base_datetime + nc_datetimes,
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
  date_seq <- seq(
    lubridate::floor_date(as.Date(start_time), "month"),
    lubridate::ceiling_date(as.Date(end_time), "month") - lubridate::days(1), by = "month"
    )
  df <- data.frame(
    mon = lubridate::month(date_seq),
    yea = lubridate::year(date_seq)
  )
  return(df)
}

#' Function to find the longest shared substring among a vector of strings, used
#'to identify file_prefix for a list of requests
#' @param strings a list or vector of strings to search for a common substring
#' @noRd
shared_substring <- function(strings) {

  # Helper function to get all substrings of a string
  get_substrings <- function(string) {
    substrings <- character(0)
    n <- nchar(string)

    for (i in 1:n) {
      for (j in i:n) {
        substrings <- c(substrings, substr(string, i, j))
      }
    }

    return(substrings)
  }

  # Start with substrings of the first string in the vector
  common_substrings <- get_substrings(strings[1])

  # Intersect common substrings with each subsequent string in the vector
  for (i in 2:length(strings)) {
    substrings_i <- get_substrings(strings[i])
    common_substrings <- intersect(common_substrings, substrings_i)

    # If no common substring is found at any point, return an empty string
    if (length(common_substrings) == 0) {
      return("")
    }
  }

  # If there are common substrings, return the longest one
  longest_substring <- common_substrings[which.max(nchar(common_substrings))]
  return(longest_substring)
}

#' Function to extract the time dimension from a `ncdf4` object. Stays flexible
#' with exact name of time dimension to support both old and new CDS
#' @param nc a `ncdf4` object
#' @noRd
extract_timedim <- function(nc) {
  # Extract time dimension
  # Specifically: pull out the first dimension that has 'time' in its name
  return(nc$dim[grepl("time", names(nc$dim))][[1]])
}

#' Function to bind a series of netCDFs, stored inside a .zip file, which all
#' have the same spatial extent and time dimension, but different sets of variables
#' @param nc_zip character vector of the file path and name of zip file containing the netCDF files downloaded from the CDS. Must end with ".zip"
#' @param combined_name character vector of the file path and desired name of combined netCDF file. Must end with ".nc"
#' @noRd
bind_zipped_netcdf <- function(nc_zip, combined_name) {
  # Confirm that nc_zip exists
  if (!file.exists(nc_zip)) {
    stop("File does not exist at the path provided to `nc_zip`")
  }
  # Check that nc_zip includes ".zip"
  if (substr(nc_zip, nchar(nc_zip)-4+1, nchar(nc_zip)) != ".zip") {
    stop("Value provided to argument `nc_zip` must end with .zip")
  }
  # Check that combined_name includes ".nc"
  if (substr(combined_name, nchar(combined_name)-3+1, nchar(combined_name)) != ".nc") {
    stop("Value provided to argument `combined_name` must end with .nc")
  }
  unzip(nc_zip, exdir = gsub(".zip", "", nc_zip))
  filenames <- list.files(gsub(".zip", "", nc_zip), full.names = TRUE)
  files <- lapply(filenames, function(x) {
    ncdf4::nc_open(x)
  })
  # Pull out first file for reference specs
  nc_example <- ncdf4::nc_open(filenames[1])
  lat <- ncdf4::ncvar_get(nc_example, "latitude")
  lon <- ncdf4::ncvar_get(nc_example, "longitude")
  time_dim <- extract_timedim(nc_example)
  ncdf4::nc_close(nc_example)

  # Define the dimensions
  lat_dim <- ncdf4::ncdim_def("latitude", "degrees_north", lat)
  lon_dim <- ncdf4::ncdim_def("longitude", "degrees_east", lon)
  time_dim <- ncdf4::ncdim_def(time_dim$name,
                               time_dim$units,
                               time_dim$vals)

  # Create an empty list to populate
  vars_list <- list()
  data_list <- list()

  # Loop through input files to extract variables and their metadata
  for (file in filenames) {
    nc <- ncdf4::nc_open(file)
    varnames <- names(nc$var)

    # Remove unnecessary vars
    varnames <- varnames[!grepl("number|expver", varnames)]

    for (var_name in varnames) {
      # Extract variable data and attributes
      data_list[[var_name]] <- ncdf4::ncvar_get(nc, var_name)
      var_unit <- nc$var[var_name][[var_name]]$units
      var_longname <- nc$var[[var_name]]$longname

      # Define the variable
      vars_list[[var_name]] <- ncdf4::ncvar_def(
        name = var_name,
        units = var_unit,
        dim = list(lon_dim, lat_dim, time_dim),
        longname = var_longname,
        missval = nc$var[[var_name]]$missval
      )
    }

    ncdf4::nc_close(nc)
  }

  # Create the new netCDF file and write the variables
  file_combined <- ncdf4::nc_create(combined_name, vars = vars_list)

  # Write the data into the new file
  for (var_name in names(data_list)) {
    ncdf4::ncvar_put(
      nc = file_combined,
      varid = var_name,
      vals = data_list[[var_name]])
  }

  # Close the output file
  ncdf4::nc_close(file_combined)

}
