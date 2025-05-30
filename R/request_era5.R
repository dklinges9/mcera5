#' Submit request(s) to CDS for ERA5 data download
#'
#' @description Uses the `ecwmfr` package to submit requests to the Climate
#' Data Store (CDS) for ERA5 .nc files. Executes one request at a time.
#'
#' @param request a list of request(s) created with `build_era5_request`.
#' @param uid character vector containing your CDS user ID.
#' @param out_path character vector of the location at which to download nc files.
#' @param overwrite TRUE for overwriting a file of the same path. Default FALSE.
#' @param combine TRUE for combining downloaded files into one file. Default is
#' TRUE as current behavior is to download multiple monthly files, to stay below
#' CDS API limits.
#'
#' @export
#'
#'
request_era5 <- function(request, uid, out_path, overwrite = FALSE,
                         combine = TRUE) {

  if (length(request) == 1 & combine) {
    cat("Your request will all be queried at once and does not need to be combined.\n")
  }

  for(req in 1:length(request)) {

    # Check whether file already exists for requested out_path
    if (file.exists(paste0(out_path, "/", request[[req]]$target)) & !overwrite) {
      if (length(request) > 1) {
        stop("Filename already exists within requested out_path in request ", req, " of request series. Use overwrite = TRUE if you wish to overwrite this file.")
      } else {
        stop("Filename already exists within requested out_path. Use overwrite = TRUE if you wish to overwrite this file.")
      }
    }

    ecmwfr::wf_request(user = as.character(uid),
                       request = request[[req]],
                       transfer = TRUE,
                       path = out_path,
                       verbose = TRUE,
                       time_out = 18000)

    if (file.exists(paste0(out_path, "/", request[[req]]$target))) {
      if (length(request) > 1) {
        cat("ERA5 netCDF file", req, "successfully downloaded.\n")
      } else {
        cat("ERA5 netCDF file successfully downloaded.\n")
      }

      # If target file is a .zip, but request isn't timeseries, then unzip and bind together netCDF files
      if (grepl(".zip", request[[req]]$target) & !grepl("timeseries", request[[req]]$dataset_short_name)) {
        cat("Downloaded file is a .zip, so unzipping and binding together contents...\n")
        # Designate input path (nc_zip) and output path (combined_name) that are
        # handed to bind_zippped_netcdf()
        nc_zip <- paste0(out_path, "/", request[[req]]$target)
        combined_name <- paste0(out_path, "/", gsub(".zip", ".nc", request[[req]]$target))
        bind_zipped_netcdf(nc_zip, combined_name)
      }

      # if (grepl(".zip", request[[req]]$target) & grepl("timeseries", request[[req]]$dataset_short_name)) {
      #   cat("Downloaded file is a .zip, so unzipping contents...\n")
      #   unzip(nc_zip, exdir = out_path)
      # }
    }
  }
  if (length(request) > 1 & combine) {
    cat("Now combining netCDF files...\n")
    # Get list of filenames
    fnames <- lapply(request, function(x) {
      # Replace .zip with .nc, in case targets were .zip files
      x$target <- gsub(".zip", ".nc", x$target)
      return(x$target)
    })
    # Recover file prefix from filenames
    file_prefix <- shared_substring(fnames)
    # If last character is "_", remove
    file_prefix <- ifelse(substr(file_prefix, nchar(file_prefix), nchar(file_prefix)) == "_",
                          substr(file_prefix, 1, nchar(file_prefix) - 1),
                          file_prefix)
    # Combine netCDF files
    combine_netcdf(filenames = paste0(out_path, "/", fnames),
                   combined_name = paste0(out_path, "/", file_prefix, ".nc"))
    cat("Finished.\n")
  }
}
