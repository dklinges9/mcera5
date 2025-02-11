
library(terra)
library(mcera5)
my_nc <- "user_files/era5land_20241224_2000.nc"
# my_nc <- "era5_20250109_2000.nc"
nc <- rast(my_nc)
extr <- ext(nc)
x <- mean(c(extr[1], extr[2]))
y <- mean(c(extr[3], extr[4]))
source("R/internal.R")
# Open nc file to get datetimes
nc_dat = ncdf4::nc_open(my_nc)
timedim <- extract_timedim(nc_dat)
# Find basetime from units
base_datetime <- as.POSIXct(gsub(".*since ", "", timedim$units), tz = "UTC")
# Extract time values
nc_datetimes <- c(timedim$vals)
# If units in hours, multiply by 3600 to convert to seconds
nc_datetimes <- nc_datetimes * ifelse(
  grepl("hours", timedim$units), 3600, 1
)
# Find first timestep
first_timestep <- nc_datetimes[1]
last_timestep <- tail(nc_datetimes[6])

st_time <- base_datetime + first_timestep
en_time <- base_datetime + last_timestep

x <- -75
y <- 42
out <- extract_clim(nc = my_nc, long = x, lat = y,
                         start_time = st_time, end_time = en_time,
                         d_weight = TRUE)

out <- extract_land_clim(nc = my_nc, long = x, lat = y,
                         start_time = st_time, end_time = en_time,
                         d_weight = TRUE)

out <- extract_land_clim(nc = my_nc, long = x, lat = y,
                         start_time = st_time, end_time = en_time,
                         d_weight = FALSE)

head(out)
