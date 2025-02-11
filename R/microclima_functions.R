#' calculates clearsky radiation
#'
#' @description
#' `clearskyrad` is used to calculate clear-sky shortwave irradiance using the
#'  Crawford & Duchon (1999) method
#' @param tme a single value or vector of POSIXlt objects indicating the time(s)
#' for which clearksy radiation is required.
#' @param lat latitude in decimal degrees
#' @param long longitude in decimal degrees
#' @param h an optional single value or vector of specific humidities (\ifelse{html}{\out{kg kg<sup>{-1}</sup> }}{\eqn{kg kg^{-1}}}).
#' @param tc an optional single value or vector of temperatures (ºC).
#' @param p an optional single value or vector of pressures (Pa).
#' @param G an optional single value or vector describing he moisture profile
#' in the atmosphere (per Smith 1966).
#' @param Ie an optional single value for extra-terrestrail radiation to permit adjustment for
#' sun-earth distances (see details).
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @seealso The function [cloudfromrad()] uses this function to return 1 - ratio of measured to clearsky
#' radiation for input when computing longwave radiation when using [longwavetopo()] or [longwaveveg()]
#' @export
#'
#' @details
#' The units returned are the same as for `Ie`, with the default option in W / M^2.
#' If no values for `p` are provided, a default value of 101300 Pa, typical of
#' sea-level pressure, is assumed. The method used is that detailed in
#' Crawford & Duchon (1999) Quarterly Journal of the Royal Meteorological
#' Society 122: 1127-1151. The method is not greatly sensitive to humidity,
#' temperature and pressure so approximate values can be provided, or the defaults
#' chosen, if these data are unavailable.
#'
#' @return a single value or vector of clearksy radiation.
#'
#' @examples
#' tme <- as.POSIXlt(c(0:23) * 3600, origin = "2010-05-23 00:00", tz = "GMT")
#' Io <- clearskyrad(tme, 50, -5, 0.007953766 , 11)
#' plot(Io ~ as.POSIXct(tme), type = "l")
clearskyrad <- function(tme, lat, long, h = 0.00697, tc = 15, p = 101300, G = 2.78, Ie = 1352.778,
                        merid = round(long/15, 0) * 15, dst = 0) {
  jd <- julday(tme$year + 1900, tme$mon + 1, tme$mday)
  lt <- tme$hour + tme$min / 60 + tme$sec / 3600
  sa <- solalt(lt, lat, long, jd, merid, dst)
  sa[sa < 0] <- NA
  z <- (90 - sa) * (pi / 180)
  m <- 35 * cos(z) * ((1224 * cos(z)^2 + 1)^(-0.5))
  TrTpg <- 1.021 - 0.084 * (m * 0.000949 * 0.01 * p + 0.051)^0.5
  pk <- p / 1000
  e0 <- 0.6108 * exp(17.27 * tc/(tc + 237.3))
  ws <- 0.622 * e0/pk
  rh <- (h/ws) * 100
  rh <- ifelse(rh > 100, 100, rh)
  xx <- log(rh / 100) + ((17.27 * tc) / (237.3 + tc))
  Td <- (237.3 * xx) / (17.27 - xx)
  u <- exp(0.1133 - log(G + 1) + 0.0393 * Td)
  Tw <- 1 - 0.077 * (u * m) ^ 0.3
  Ta <- 0.935 * m
  od <- TrTpg * Tw * Ta
  Ic <- Ie * (cos(z)) * od
  Ic[Ic > Ie] <- NA
  Ic[Ic < 0] <- NA
  Ic
}

#' Applies height correction to wind speed measurements
#'
#' @description `windheight` is used to to apply a height correction to wind speed measured at a specified height above ground level to obtain estimates of wind speed at a desired height above the ground.
#'
#' @param ui numeric value(s) of measured wind speed (\ifelse{html}{\out{m s<sup>-1</sup> }}{\eqn{ m s^{-1}}}) at height `zi` (m).
#' @param zi a numeric value idicating the height (m) above the ground at which `ui` was measured.
#' @param zo a numeric value indicating the height (m) above ground level at which wind speed is desired.
#'
#' @details Thus function assumes a logarithmic height profile to convert
#' wind speeds. It performs innacurately when `uo` is lower than 0.2
#' and a warning is given. If `uo` is below ~0.08 then the logairthmic height
#' profile cannot be used, and `uo` is converted to 0.1 and a warning given.
#'
#' @seealso The function [windcoef()] calculates a topographic or vegetation sheltering effect.
#'
#' @return numeric values(s) of wind speed at height specified by `zo` (\ifelse{html}{\out{m s<sup>-1</sup> }}{\eqn{m s^{-1}}}).
#' @export
#'
#' @examples
#' windheight(3, 10, 1) # good
#' windheight(3, 10, 0.15) # performs poorly. Warning given
#' windheight(3, 10, 0.05) # cannot calculate. ui converted and warning given
#'
windheight <- function(ui, zi, zo) {
  if (zo < 0.2 & zo > (5.42 / 67.8)) {
    warning("Wind-height profile function performs poorly when wind
            height is below 20 cm")
  }
  if (zo <= (5.42 / 67.8)) {
    warning(paste("wind-height profile cannot be calculated for height ",
                  zo * 100, " cm"))
    print("Height converted to 10 cm")
    zo <- 0.1
  }
  uo <- ui * log(67.8 * zo - 5.42) / log(67.8 * zi - 5.42)
  uo
}

#' Calculates the astronomical Julian day
#'
#' @description `julian` is used to calculate the astronomical Julian day (days since since January 1, 4713 BCE at noon UTC) from a given year, month and day.
#'
#' @param year year (AD).
#' @param month month in numeric form (1-12).
#' @param day days of the month (1-31).
#' @param hour hours (decimal, 0-23).
#' @param min minutes (decimal, 0-59).
#' @param sec seconds (decimal, 0-59).
#' @param dst an optional numeric value specifying the time zone expressed as hours different from GMT (-ve to west).
#'
#' @return Julian Day. I.e. the number of days since January 1, 4713 BCE at noon UTC.
#' @export
#'
#' @examples
#' jd1 <- julday(2010, 1, 31)
#' jd2 <- julday(2010, 1, 31, 11, 0, 0)
#' jd1 - jd2
julday <- function(year, month, day, hour = 12, min = 0, sec = 0, dst = 0) {
  day_decimal <- day + (hour - dst + (min + sec / 60) / 60) / 24
  monthadj <- month + (month < 3) * 12
  yearadj <- year + (month < 3) * -1
  julian_day <- trunc(365.25 * (yearadj + 4716)) + trunc(30.6001 *
                                                           (monthadj + 1)) + day_decimal - 1524.5
  B <- (2 - trunc(yearadj / 100) + trunc(trunc(yearadj / 100) / 4))
  julian_day <- julian_day + (julian_day > 2299160) * B
  julian_day
}

#' Calculates the solar index for a flat surface
#'
#' @param localtime local time (decimal hour, 24 hour clock).
#' @param lat latitude of the location for which the solar altitude is required (decimal degrees, -ve south of the equator).
#' @param long longitude of the location for which the solar altitude is required (decimal degrees, -ve west of Greenwich meridian).
#' @param julian Julian day expressed as an integer as returned by [julday()].
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#'
#' @return a numeric value representing the solar altitude (º).
#' @export
#'
#' @examples
#' # solar index at noon GMT on 21 June 2010, Porthleven, Cornwall
#' jd <- julday (2010, 6, 21) # Julian day
#' siflat(12, 50.08, -5.31, jd)
#'
siflat <- function(localtime, lat, long, julian, merid = round(long / 15, 0) * 15, dst = 0){
  saltitude <- solalt(localtime, lat, long, julian, merid, dst)
  alt <- saltitude * (pi/180)
  index <- cos(pi/2 - alt)
  index[index < 0] <- 0
  index
}

#' Calculates the solar altitude
#'
#' @description `solalt` is used to calculate the solar altitude at any given location from the local time.
#'
#' @param localtime local time (decimal hour, 24 hour clock).
#' @param lat latitude of the location for which the solar altitude is required (decimal degrees, -ve south of the equator).
#' @param long longitude of the location for which the solar altitude is required (decimal degrees, -ve west of Greenwich meridian).
#' @param julian Julian day expressed as an integer as returned by [julday()].
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#'
#' @return a numeric value representing the solar altitude (º).
#' @export
#'
#' @examples
#' # solar altitude at noon on 21 June 2010, Porthleven, Cornwall
#' jd <- julday (2010, 6, 21) # Julian day
#' solalt(12, 50.08, -5.31, jd)
solalt <- function(localtime, lat, long, julian, merid = round(long / 15, 0) * 15, dst = 0) {
  stime <- solartime(localtime, long, julian, merid, dst)
  tt <- 0.261799 * (stime - 12)
  declin <- (pi * 23.5 / 180) * cos(2 * pi * ((julian - 171) / 365.25))
  sinh <- sin(declin) * sin(lat * pi / 180) + cos(declin) *
    cos(lat * pi / 180) * cos(tt)
  sa <- (180 * atan(sinh / sqrt(1 - sinh * sinh))) / pi
  sa
}

#' Calculates the solar time
#'
#' @description `solartime` is used to calculate the solar time. I.e. the time that would be measured by a sundial.
#'
#' @param localtime local time (decimal hour, 24 hour clock).
#' @param long longitude of the location for which the solar time is required (decimal degrees, -ve west of Greenwich meridian).
#' @param julian Julian day expressed as an integer as returned by [julday()].
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT).
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#'
#' @return the solar time. I.e. the times that would be measured by a sundial (hours).
#' @export
#'
#' @details
#' ‘solartime’ accounts for two factors: firstly, east or west component of the analemma, namely
#' the angular offset of the Sun from its mean position on the celestial sphere as viewed from
#' Earth due the eccentricity of the Earth's orbit and the obliquity due to tilt of the Earth's
#' rotational axis. These two factors have different wavelengths, amplitudes and phases,
#' that vary over geological timescales. The equations used here are those derived by Milne. BY default,
#' local are used, with the meridian set to round(long / 15, 0) * 15.
#'
#'
#' @examples
#' jd <- julday (2010, 6, 21) # Julian day
#' solartime(12, -5, jd) # solartime at noon on 21 June 2010, 5ºW
solartime <- function(localtime, long, julian, merid = round(long / 15, 0) * 15, dst = 0) {
  m <- 6.24004077 + 0.01720197 * (julian -  2451545)
  eot <- -7.659 * sin(m) + 9.863 * sin (2 * m + 3.5932)
  st <- localtime + (4 * (long - merid) + eot) / 60 - dst
  st
}


#' Convert between different measures of humidity
#'
#' @description converts between different humidity measures, namely
#' vapour pressure or relative, absolute or specific humidity.
#'
#' @param h humidity value(s). Units as follows: specific humidity (kg/kg),
#' absolute humidity (kg/m^3),
#' relative humidity (Percentage), vapour pressure (kPa).
#' @param intype a character string description of the humidity type of `h`. One of "relative", "absolute" or "specific".
#' @param tc A numeric value specifying the temperature (deg C).
#' @param pk An optional numeric value specifying the atmospheric pressure (kPa).
#'
#' @details This function converts between vapour pressure and specific,
#' relative and absolute humidity, based on sea-level pressure and
#' temperature. It returns a list of relative, absolute and specific
#' humidity and vapour pressure. If one or more of the
#' relative humidity values exceeds 100\% a warning is given.
#'
#' @return a list of numeric humidity values with the following components:
#' @return `relative` relative humidity (Percentage).
#' @return `absolute`  absolute humidity (kg/m^3).
#' @return `specific` specific humidity (kg/kg).
#' @return `vapour_pressure` vapour pressure (kPa).
#' @export
#'
#' @examples
#' converthumidity(90, 'relative', 20)
#' converthumidity(0.01555486, 'absolute', 20)
#' converthumidity(0.01292172, 'specific', 20)
converthumidity <- function(h, intype = "relative", tc = 11, pk = 101.3) {
  tk <- tc + 273.15
  if (intype != "specific" & intype != "relative" & intype !=
      "absolute") {
    warning("No valid input humidity specified. Humidity assumed to be\n relative")
    intype <- "relative"
  }
  e0 <- 0.6108 * exp(17.27 * tc/(tc + 237.3))
  ws <- (18.02 / 28.97) * (e0 / pk)
  if (intype == "specific") {
    hr <- (h/ws) * 100
  }
  if (intype == "absolute") {
    ea <- (tk * h) / 2.16679
    hr <- (ea/e0) * 100
  }
  if (intype == "relative")
    hr <- h
  if (max(hr, na.rm = T) > 100)
    warning(paste("Some relative humidity values > 100%",
                  max(hr, na.rm = T)))
  if (intype == "vapour pressure") {
    hr <- (ea / e0) * 100
  }
  ea <- e0 * (hr / 100)
  hs <- (hr / 100) * ws
  ha <- 2.16679 * (ea/tk)
  return(list(relative = hr, absolute = ha, specific = hs,
              vapour_pressure = ea))
}

#' Calculates the diffuse fraction from incoming shortwave radiation
#'
#' @description `difprop` calculates proportion of incoming shortwave radiation that is diffuse radiation using the method of Skartveit et al. (1998) Solar Energy, 63: 173-183.
#'
#' @param rad a vector of incoming shortwave radiation values (either \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} or \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}})
#' @param julian the Julian day as returned by [julday()]
#' @param localtime a single numeric value representing local time (decimal hour, 24 hour clock)
#' @param lat a single numeric value representing the latitude of the location for which partitioned radiation is required (decimal degrees, -ve south of equator).
#' @param long a single numeric value representing the longitude of the location for which partitioned radiation is required (decimal degrees, -ve west of Greenwich meridian).
#' @param hourly specifies whether values of `rad` are hourly (see details).
#' @param watts a logical value indicating  whether the units of `rad` are \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}} (TRUE) or \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} (FALSE).
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param corr an optional numeric value representing a correction to account for over- or under-estimated diffuse proportions. Values > 1 will apportion a greater ammount of total radiation as diffuse than originally calculated by the formula.
#'
#' @return a vector of diffuse fractions (either \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} or \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}}).
#' @export
#'
#' @details
#' The method assumes the environment is snow free. Both overall cloud cover and heterogeneity in
#' cloud cover affect the diffuse fraction. Breaks in an extensive cloud deck may primarily
#' enhance the beam irradiance, whereas scattered clouds may enhance the diffuse irradiance and
#' leave the beam irradiance unaffected.  In consequence, if hourly data are available, an index
#' is applied to detect the presence of such variable/inhomogeneous clouds, based on variability
#' in radiation for each hour in question and values in the preceding and deciding hour.  If
#' hourly data are unavailable, an average variability is determined from radiation intensity.
#'
#' @examples
#' rad <- c(5:42) / 0.036 # typical values of radiation in W/m^2
#' jd <- julday(2017, 6, 21) # julian day
#' dfr <- difprop(rad, jd, 12, 50, -5)
#' plot(dfr ~ rad, type = "l", lwd = 2, xlab = "Incoming shortwave radiation",
#'      ylab = "Diffuse fraction")
difprop <- function(rad, julian, localtime, lat, long, hourly = TRUE,
                    watts = TRUE, merid = round(long / 15, 0) * 15, dst = 0,
                    corr = 1) {
  if (watts) rad <- rad * 0.0036
  sa <- solalt(localtime, lat, long, julian, merid, dst)
  alt <- sa * (pi / 180)
  k1 <- 0.83 - 0.56 * exp(- 0.06 * sa)
  si <- cos(pi / 2 - alt)
  si[si < 0] <- 0
  k <- rad / (4.87 * si)
  k[is.na(k)] <- 0
  k <- ifelse(k > k1, k1, k)
  k[k < 0] <- 0
  rho <- k / k1
  if (hourly) {
    rho <- c(rho[1], rho, rho[length(rho)])
    sigma3  <- 0
    for (i in 1:length(rad)) {
      sigma3[i] <- (((rho[i + 1] - rho[i]) ^ 2 + (rho[i + 1] - rho[i + 2]) ^ 2)
                    / 2) ^ 0.5
    }
  } else {
    sigma3a <- 0.021 + 0.397 * rho - 0.231 * rho ^ 2 - 0.13 *
      exp(-1 * (((rho - 0.931) / 0.134) ^ 2) ^ 0.834)
    sigma3b <- 0.12 + 0.65 * (rho - 1.04)
    sigma3 <- ifelse(rho <= 1.04, sigma3a, sigma3b)
  }
  k2 <- 0.95 * k1
  d1 <- ifelse(sa > 1.4, 0.07 + 0.046 * (90 - sa) / (sa + 3), 1)
  K <- 0.5 * (1 + sin(pi * (k - 0.22) / (k1 - 0.22) - pi / 2))
  d2 <- 1 - ((1 - d1) * (0.11 * sqrt(K) + 0.15 * K + 0.74 * K ^ 2))
  d3 <- (d2 * k2) * (1 - k) / (k * (1 - k2))
  alpha <- (1 / sin(alt)) ^ 0.6
  kbmax <- 0.81 ^ alpha
  kmax <- (kbmax + d2 * k2 / (1 - k2)) / (1 + d2 * k2 / (1 - k2))
  dmax <- (d2 * k2) * (1 - kmax) / (kmax * (1 - k2))
  d4 <- 1 - kmax * (1 - dmax) / k
  d <- ifelse(k <= kmax, d3, d4)
  d <- ifelse(k <= k2, d2, d)
  d <- ifelse(k <= 0.22, 1, d)
  kX <- 0.56 - 0.32 * exp(-0.06 * sa)
  kL <- (k - 0.14) / (kX - 0.14)
  kR <- (k - kX) / 0.71
  delta <- ifelse(k >= 0.14 & k < kX, -3 * kL ^ 2 *(1 - kL) * sigma3 ^ 1.3, 0)
  delta <- ifelse(k >= kX & k < (kX + 0.71), 3 * kR * (1 - kR) ^ 2 * sigma3 ^
                    0.6, delta)
  d[sigma3 > 0.01] <- d[sigma3 > 0.01] + delta[sigma3 > 0.01]
  d[rad == 0] <- 0.5
  d[sa < 0] <- 1
  # apply correction
  dif_val <- rad * d
  dif_val_adj <- dif_val * corr
  d <- dif_val_adj /rad
  d[d > 1] <- 1
  d[d < 0] <- 1
  d[is.na(d)] <- 0.5
  d
}

#' calculates cloud cover form shortwave radiation
#'
#' @description
#' `cloudfromrad` is used to derive a cloud cover index as 1 - the ratio of measured to clearksy radiation
#' missing values (e.g. at night) are interpolated
#' @param rad a single numeric value or vector of global solar irradiance(s). Units must be the same
#' as those for `Ie`. Default is Watts /m^2.
#' @param tme a single value or vector of POSIXlt objects indicating the time(s)
#' for which clearksy radiation is required.
#' @param lat latitude in decimal degrees
#' @param long longitude in decimal degrees
#' @param h a single value or vector of specific humidities (\ifelse{html}{\out{kg kg<sup>{-1}</sup> }}{\eqn{kg kg^{-1}}}).
#' @param tc a single value or vector of temperatures (ºC).
#' @param p an optional single value or vector of pressures (Pa).
#' @param G an optional single value or vector describing he moisture profile
#' in the atmosphere (per Smith 1966).
#' @param Ie an optional single value for extra-terrestrial radiation to permit adjustment for
#' sun-earth distances (see details).
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @import zoo
#' @export
#'
#' @seealso The function [clearskyrad()] is uses to derive clear-sky irradiance. Can be used
#' to derive cloud covers for computing longwave radiation when using [longwavetopo()] or [longwaveveg()]
#'
#' @details
#' If no values for `p` are provided, a default value of 101300 Pa, typical of
#' sea-level pressure, is assumed. The method used is that detailed in
#' Crawford & Duchon (1999) Quarterly Journal of the Royal Meteorological
#' Society 122: 1127-1151. The method is not greatly sensitive to humidity,
#' temperature and pressure so approximate values can be provided, or the defcloudfromradaults
#' chosen, if these data are unavailable.
#'
#' @return a single value or vector of clearksy radiation.
#'
#' @examples
#' tme <- as.POSIXlt(c(0:23) * 3600, origin = "2010-05-23 00:00", tz = "GMT")
#' rad <- clearskyrad(tme, 50, -5, 0.007953766 , 11) * 0.75
#' cfc <- cloudfromrad(rad, tme, 50, -5, 0.007953766 , 11)
#' plot(cfc ~ as.POSIXct(tme), type = "l") # should be 0.25
cloudfromrad <- function(rad, tme, lat, long, h = 0.00697, tc = 15, p = 101300, G = 2.78,
                         Ie = 1352.778, merid = round(long/15, 0) * 15, dst = 0) {
  Ic <- clearskyrad(tme, lat, long, h, tc, p, G, Ie, merid, dst)
  s <- rad / Ic
  s[s > 1] <- 1
  s[s < 0] <- 0
  if (is.na(s[1])) s[1] <- mean(s, na.rm = T)
  if (is.na(s[length(s)])) s[length(s)] <- mean(s, na.rm = T)
  # Replaced method `s <- na.approx(s)` used in microclima to avoid a dependency on
  # package `zoo`
  s <- approx(seq_along(s), s, xout = seq_along(s), method = "linear")$y
  cfc <- 1 - s
  cfc
}
