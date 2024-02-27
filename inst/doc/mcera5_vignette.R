## ----setup, include=FALSE-------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      collapse = TRUE,
                      comment = "#>")
# This option allows the displayed title ("Downloading ERA5 data for use...") to
# be different from the \VignetteIndexEntry{} title:
options(rmarkdown.html_vignette.check_title = FALSE)

## ----packages, warning = FALSE, message = FALSE---------------------------------------------------------------
library(mcera5)
library(dplyr)
library(ecmwfr)
library(ncdf4)
library(curl)
library(keyring)
library(abind)
library(lubridate)
library(tidync)
library(microclima) # remotes::install_github("ilyamaclean/microclima")
library(NicheMapR) # remotes::install_github("mrke/NicheMapR")
library(microctools) # remotes::install_github("ilyamaclean/microctools")

## ----funs, include = FALSE------------------------------------------------------------------------------------
#############
##FUNCTIONS##
#############

files_source <- list.files(here::here("R/"), full.names = T)
sapply(files_source,source)

## ----creds, eval = FALSE--------------------------------------------------------------------------------------
#  # assign your credentials
#  uid <- "*****"
#  cds_api_key <- "********-****-****-****-************"
#  
#  ecmwfr::wf_set_key(user = uid,
#                     key = cds_api_key,
#                     service = "cds")

## ----build request--------------------------------------------------------------------------------------------

# bounding coordinates (in WGS84 / EPSG:4326)
xmn <- -4
xmx <- -2
ymn <- 49
ymx <- 51

# temporal extent
st_time <- as.POSIXlt("2010-02-26 00:00", tz = "UTC")
en_time <- as.POSIXlt("2010-03-01 23:00", tz = "UTC")


# Set a unique prefix for the filename (here based on spatial
# coordinates), and the file path for downloaded .nc files
file_prefix <- "era5_-4_-2_49_51"
file_path <- getwd()


# build a request (covering multiple years)
req <- build_era5_request(xmin = xmn, xmax = xmx, 
                          ymin = ymn, ymax = ymx,
                          start_time = st_time,
                          end_time = en_time,
                          outfile_name = file_prefix)

## ----list_view------------------------------------------------------------------------------------------------
str(req)

## ----send_request, eval = FALSE-------------------------------------------------------------------------------
#  request_era5(request = req, uid = uid, out_path = file_path)

## ----process_clim, eval = FALSE-------------------------------------------------------------------------------
#  # list the path of the .nc file for a given year
#  my_nc <- paste0(getwd(), "/era5_-4_-2_49_51_2010.nc")
#  
#  # for a single point (make sure it's within the bounds of your .nc file)
#  x <- -3.654600
#  y <- 50.640369
#  
#  # gather all hourly variables
#  clim_point <- extract_clim(nc = my_nc, long = x, lat = y,
#                              start_time = st_time, end_time = en_time)
#  head(clim_point)

## ----process_precip, eval = FALSE-----------------------------------------------------------------------------
#  # gather daily precipitation (we specify to convert precipitation from hourly
#  # to daily, which is already the default behavior)
#  precip_point <- extract_precip(nc = my_nc, long = x, lat = y,
#                                     start_time = st_time,
#                                     end_time = en_time,
#                                     convert_daily = TRUE)

## ----runauto_example, eval = FALSE----------------------------------------------------------------------------
#  # create a 200 x 200 30 m spatial resoltuion DEM for location
#  r <- microclima::get_dem(lat = y, long = x, resolution = 30)
#  
#  # change date format to fit runauto requirements
#  temps <- microclima::runauto(r = r, dstart = "26/02/2010",dfinish = "01/03/2010",
#                               hgt = 0.1, l = NA, x = NA,
#                               habitat = "Barren or sparsely vegetated",
#                               hourlydata = clim_point,
#                               dailyprecip = precip_point,
#                               plot.progress= FALSE)

## ----hourlyncep_convert_example, eval = FALSE-----------------------------------------------------------------
#  climdata <- hourlyncep_convert(climdata = clim_point, lat = y, long = x)

## ----process_clima, eval = FALSE------------------------------------------------------------------------------
#  # list the path of the .nc file for a given year
#  my_nc <- paste0(getwd(), "/era5_-4_-2_49_51_2010.nc")
#  
#  # 4 corners of a spatial grid (make sure it's within the bounds of your .nc file)
#  long_min <- -3.7
#  long_max <- -2.9
#  lat_min <- 50.1
#  lat_max <- 50.8
#  
#  # gather all hourly variables
#  clim_grid <- extract_clima(nc, long_min, long_max, lat_min, lat_max,
#                            start_time = st_time,
#                            end_time = en_time,
#                            dtr_cor = TRUE, dtr_cor_fac = 1.285,
#                            reformat = TRUE)
#  
#  str(clim_grid)
#  
#  
#  precip_grid <- extract_precipa(nc, long_min, long_max, lat_min, lat_max, start_time,
#                               end_time, convert_daily = TRUE)
#  
#  str(precip_grid)
#  

