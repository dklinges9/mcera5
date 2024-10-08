#' Produces hourly data across a grid ready for use with
#' several gridded microclimate models.
#'
#' @description `extract_clima` takes an nc file containing hourly ERA5 climate
#' data, and for a given set of coordinates, reconverts climate variables to make
#' ready for use with `microclimf::modelina` by default. Also provides the option
#' to implement a diurnal temperature range correction to air temperatures.
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
#' @param dtr_cor logical value indicating whether to apply a diurnal temperature
#' range correction to air temperature values. Default = `TRUE`.
#' @param dtr_cor_fac numeric value to be used in the diurnal temperature range
#' correction. Default = 1.285, based on calibration against UK Met Office
#' observations.
#' @param format specifies what microclimate package extracted climate data will
#' be used for. Data will be formatted accordingly. Default is "microclimf".
#' Options: "microclima", "NicheMapR", "microclimc", "microclimf", "micropoint".
#' Note: of these options, only "microclimf" accepts as input an array of climate
#' variables. For all other models you will need to iterate through each spatial
#' point to run the model
#' @param cds_version specifies what version of the CDS (Climate Data Store) the
#' ERA5 data were queried from. Either "new" (queries made after Sept 2024) or
#' "legacy" (queries made before Sept 2024). This argument will be deprecated in
#' the future and only queries from the new CDS will be supported.
#'
#' @return Returns a list of wrapped spatRasters containing hourly values for a
#' suite of climate variables. The returned climate variables depends on the
#' value of argument `format`.
#'
#' If format is "microclima" or "NicheMapR":
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
#' If format is "microclimc":
#' @return `obs_time` | POSIXlt object of dates and times
#' @return `temp` | temperature (degrees C)
#' @return `relhum` | relative humidity (percentage)
#' @return `pres` | atmospheric press (kPa)
#' @return `swrad` | Total incoming shortwave radiation (W / m^2)
#' @return `difrad` | Diffuse radiation (W / m^2)
#' @return `skyem` | Sky emissivity (0-1)
#' @return `lwdown` | Total downward longwave radiation (W / m^2)
#' @return `windspeed` | Wind speed (m/s)
#' @return `winddir` | Wind direction (decimal degrees)
#' @return `precip` | precipitation (mm)
#'
#' If format is "micropoint" or "microclimf":
#' @return `obs_time` | POSIXlt object of dates and times
#' @return `temp` | temperature (degrees C)
#' @return `relhum` | relative humidity (percentage)
#' @return `pres` | atmospheric press (kPa)
#' @return `swdown` | Total incoming shortwave radiation (W / m^2)
#' @return `difrad` | Diffuse radiation (W / m^2)
#' @return `skyem` | Sky emissivity (0-1)
#' @return `lwdown` | Total downward longwave radiation (W / m^2)
#' @return `windspeed` | Wind speed (m/s)
#' @return `winddir` | Wind direction (decimal degrees)
#' @return `precip` | precipitation (mm)
#'
#' @export
extract_clima <- function(
    nc, long_min, long_max, lat_min, lat_max, start_time, end_time,
    dtr_cor = TRUE, dtr_cor_fac = 1.285,
    format = "microclimf", cds_version = "new") {

  # Open nc file for error trapping
  nc_dat = ncdf4::nc_open(nc)

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
  start <- base_datetime + (nc_dat$dim$valid_time$vals[1])
  if (start_time < start) {
    stop("Requested start time is before the beginning of time series of the ERA5 netCDF.")
  }

  # Check if end_time is before last time observation
  end <- base_datetime + (utils::tail(nc_dat$dim$valid_time$vals, n = 1))
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
  if (long_min < (min(nc_dat$dim$longitude$vals) - 0.125) |
      long_max > (max(nc_dat$dim$longitude$vals) + 0.125)
  ) {
    long_out <- TRUE
  } else {
    long_out <- FALSE
  }

  if (lat_min < (min(nc_dat$dim$latitude$vals) - 0.125) |
      lat_min > (max(nc_dat$dim$latitude$vals) + 0.125)
  ) {
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

  # Check that `format` is an accepted value
  if (!format %in% c("NicheMapR", "microclima", "microclimc", "micropoint", "microclimf")) {
    stop("Argument `format` must be one of the following values: `NicheMapR`, `microclima`, `microclimc`, `micropoint`, `microclimf`")
  }

  if (dtr_cor == TRUE & !is.numeric(dtr_cor_fac)) {
    stop("Invalid diurnal temperature range correction value provided.")
  }

  tme <- as.POSIXct(seq(start_time,
                        end_time, by = 3600), tz = lubridate::tz(end_time))

  ## Load in and subset netCDF variables --------------

  varname_list <- c("t2m", "d2m", "sp", "u10" , "v10",  "tp", "tcc", "msnlwrf",
                    "msdwlwrf", "fdir", "ssrd", "lsm")

  var_list <- lapply(varname_list, function(v) {

    if (v == "lsm") {
      # only need one timestep for land-sea mask
      r <- terra::rast(nc, subds = v)[[1]]
    } else {
      # For all others, subset down to desired time period
      r <- terra::rast(nc, subds = v)
      r <- r[[terra::time(r) %in% tme]]
    }

    # Name layers as timesteps
    names(r) <- paste(terra::time(r), lubridate::tz(terra::time(r)))
    # Subset down to desired spatial extent
    r <- terra::crop(r, terra::ext(long_min, long_max, lat_min, lat_max))
    return(r)
  })

  names(var_list) <- varname_list

  t2m <- var_list$t2m
  d2m <- var_list$d2m
  sp <- var_list$sp
  u10 <- var_list$u10
  v10 <- var_list$v10
  tcc <- var_list$tcc
  msnlwrf <- var_list$msnlwrf
  msdwlwrf <- var_list$msdwlwrf
  fdir <- var_list$fdir
  ssrd <- var_list$ssrd
  prec <- var_list$tp * 1000 # convert from mm to metres
  lsm <- var_list$lsm
  temperature <- t2m - 273.15 # kelvin to celcius
  ## Coastal correction ----------

  # see land-sea mask here:
  # https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

  # Only conduct if there are cells with proximity to water
  if (any(terra::values(lsm) < 1)) {
    # Calculate daily average
    # Indices to associate each layer with its yday
    ind <- rep(1:(dim(temperature)[3]/24), each = 24)
    # Average across days
    tmean <- terra::tapp(temperature, ind, fun = mean, na.rm = T)
    # Repeat the stack 24 times to expand back out to original timeseries
    tmean <- rep(tmean, 24)
    # Sort according to names so that the stack is now in correct order: each
    # daily mean, repeated 24 times
    # your pasted command properly sorts the names, X1 to X365 (or X366)
    tmean <- tmean[[paste0("X", sort(rep(seq(1:(dim(temperature)[3]/24)), 24)))]]
    m <- (1 - lsm) * dtr_cor_fac + 1
    tdif <- (temperature - tmean) * m
    temperature <- tmean + tdif
  }

  humidity <- humfromdew(d2m - 273.15,
                         temperature,
                         sp)
  windspeed = sqrt(u10^2 + v10^2)
  windspeed = windheight(windspeed, 10, 2)
  winddir = (terra::atan2(u10, v10) * 180/pi + 180)%%360
  cloudcover = tcc * 100
  netlong = abs(msnlwrf) * 0.0036 # Convert to MJ/m^2/hr
  downlong = msdwlwrf * 0.0036 # Convert to MJ/m^2/hr
  uplong = netlong + downlong
  emissivity = downlong/uplong
  jd = julday(lubridate::year(tme),
              lubridate::month(tme),
              lubridate::day(tme))
  rad_dni = fdir * 0.000001 # Convert form J/m^2 to MJ/m^2
  rad_glbl = ssrd * 0.000001 # Convert form J/m^2 to MJ/m^2
  ## si processing -----------------
  # use t2m as template of dimensions for iterating through
  si <- t2m
  foo <- t2m[[1]]
  coords <- as.data.frame(terra::crds(t2m[[1]]))
  # # Get a vect (points) of the scene coords
  # foo <- coordinates(raster(foo))
  # points <- terra::vect(foo, crs = terra::crs(t2m, proj = T), type = "points")
  # # And coerce points to dataframe
  # points <- terra::geom(points, df = TRUE)
  # Create a template with dimensions (x * y) and length(tme)
  out <- array(NA, dim = c(nrow(coords), length(tme)))
  for (i in 1:nrow(coords)) {
    out[i,] <- siflat(lubridate::hour(tme),
                      lat = coords$y[i],
                      long = coords$x[i],
                      jd)
  }
  si <- terra::setValues(si, out)

  # Calc rad_dif
  rad_dif = rad_glbl - rad_dni * si
  ## szenith processing -----------------
  # use t2m as template of dimensions for iterating through
  szenith <- t2m
  # Create a template with dimensions (x * y) and length(tme)
  out <- array(NA, dim = c(nrow(coords), length(tme)))
  for (i in 1:nrow(coords)) {
    out[i,] <- solalt(lubridate::hour(tme),
                      lat = coords$y[i],
                      long = coords$x[i],
                      julian = jd)
  }
  szenith <- terra::setValues(szenith, out)

  ## Add timesteps back to names ---------------
  # Only necessary for temperature at the moment, all other variables retain info
  names(temperature) <- names(t2m)
  terra::time(temperature) <- terra::time(t2m)

  ## Format of output use ----------
  ## Equivalent of hourlyncep_convert
  if (format %in% c("microclimc", "micropoint", "microclimf")) {
    pres <- sp / 1000
    ## Convert humidity from specific to relative
    relhum <- humidity
    terra::values(relhum) <- converthumidity(h = terra::as.array(humidity),
                                             intype = "specific",
                                             tc  = terra::as.array(temperature),
                                             pk = terra::as.array(pres))$relative
    relhum[relhum > 100] <- 100
    raddr <- (rad_dni * si)/0.0036 # convert back to W/m^2
    difrad <- rad_dif/0.0036 # convert from MJ/hr to W/m^2
    swrad <- raddr + difrad
    downlong <- downlong / 0.0036 # convert from MJ/hr to W/m^2
  }
  # Return list - SpatRasters now wrapped as won't store as list if saved otherwise'
  if (format %in% c("microclimc", "micropoint", "microclimf")) {
    return(list(temp = terra::wrap(temperature),
                relhum = terra::wrap(relhum),
                pres = terra::wrap(pres),
                swdown = terra::wrap(swrad),
                difrad = terra::wrap(difrad),
                lwdown = terra::wrap(downlong),
                windspeed = terra::wrap(windspeed),
                winddir = terra::wrap(winddir),
                precip = terra::wrap(prec)))
  } else {
    return(list(obs_time = tme,
                temperature = terra::wrap(temperature),
                humidity = terra::wrap(humidity),
                pressure = terra::wrap(sp),
                windspeed = terra::wrap(windspeed),
                winddir = terra::wrap(winddir),
                emissivity = terra::wrap(emissivity),
                cloudcover = terra::wrap(tcc),
                netlong = terra::wrap(netlong),
                uplong = terra::wrap(uplong),
                downlong = terra::wrap(downlong),
                rad_dni = terra::wrap(rad_dni),
                rad_dif = terra::wrap(rad_dif),
                szenith = terra::wrap(szenith)))
  }
}
