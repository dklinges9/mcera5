# This snippet resolves the notes about "no visible binding" for variable names
# when running R CMD check or devtools::check(). See this link in Hadley and
# Jenny's book for info:
# https://r-pkgs.org/package-within.html#echo-a-working-package

utils::globalVariables(c(
  ".", "yea", "mon", "abind", "obs_time", "weighted.mean", "inverse_weight",
  "precipitation", "longitude", "latitude", "time", "sp", "t2m", "d2m",
  "temperature", "pressure", "u10", "v10", "windspeed", "tcc", "msnlwrf",
  "msdwlwrf", "netlong", "downlong", "uplong", "jd", "fdir", "ssrd", "rad_glbl",
  "rad_dni", "si", "humidity", "winddir", "emissivity", "cloudcover", "rad_dif",
  "szenith", "timezone", "tp"))
