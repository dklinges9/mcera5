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
