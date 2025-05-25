#' Combines a series of netCDFs that all have the same spatial extent,
#' set of variables, and units of time
#' @param filenames a list of filenames for netCDFs you wish to combine
#' @param combined_name the name of the combined netCDF, must end with ".nc"
#'
#' @return No object is returned, the combined netCDF file is written to disk
#' with the name specified by `combined_name`
#' @export
combine_netcdf <- function(filenames, combined_name) {
  # Check that combined_name includes ".nc"
  if (substr(combined_name, nchar(combined_name)-3+1, nchar(combined_name)) != ".nc") {
    stop("Value provided to argument `combined_name` must end with .nc")
  }

  files <- lapply(filenames, function(x) {
    ncdf4::nc_open(x)
  })

  # Pull out first file for reference specs
  nc <- files[[1]]
  # Remove file metadata from vars
  varnames <- names(nc$var)
  varnames <- varnames[!grepl("number|expver", varnames)]
  # Create an empty list to populate
  vars_list <- vector(mode = "list", length = length(varnames))
  data_list <- vector(mode = "list", length = length(varnames))
  # One variable at a time
  for (i in 1:length(varnames)) {
    varname <- varnames[i]
    # Get the variable from each of the netCDFs
    vars_dat <- lapply(files, function(x) {
      ncdf4::ncvar_get(x, varname)
    })

    # Then bind all of the arrays together using abind, flexibly called via do.call
    data_list[[i]] <- do.call(abind::abind, list(
      ... = vars_dat,
      # If the ERA5 netCDF is a slice that has just 1 value of latitude, or 1 value
      # of longitude, then vars_dat will return a 2-dimensional array rather
      # than an expected 3-dimensional array. But the time dimension should always
      # best last, and it is the time dimension we specify for argument 'along'
      # So set 'along' as the length of dimensions, to pull the last dimension
      along = length(dim(vars_dat[[1]]))
    ))

    # To populate the time dimension, need to pull out the time values from each
    # netCDF
    timevals <- lapply(files, function(x) {
      extract_timedim(x)$vals
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
          extract_timedim(nc)$name, extract_timedim(nc)$units,
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
  for (i in 1:length(varnames)) {
    ncdf4::ncvar_put(
      nc = file_combined,
      varid = varnames[i],
      vals = data_list[[i]]
    )
  }

  # Finally, close the file
  ncdf4::nc_close(file_combined)
}
