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
    }
  }
  if (length(request) > 1 & combine) {
    cat("Now combining netCDF files...\n")
    # Get list of filenames
    fnames <- lapply(request, function(x) {x$target})
    # Recover file prefix from filenames
    file_prefix <- shared_substring(fnames)
    # If last character is "_", remove
    file_prefix <- ifelse(substr(file_prefix, nchar(file_prefix), nchar(file_prefix)) == "_",
                          substr(file_prefix, 1, nchar(file_prefix) - 1),
                          file_prefix)
    # Combine netCDF files
    combine_netcdf(filenames = fnames,
                   combined_name = paste0(out_path, "/", file_prefix, ".nc"))
    cat("Finished.\n")
  }
}
