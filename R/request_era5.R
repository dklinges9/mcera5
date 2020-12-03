#' Submit request(s) to CDS for ERA5 data download
#'
#' @description Uses the `ecwmfr` package to submit requests to the Climate
#' Data Store (CDS) for ERA5 .nc files. Executes one request at a time.
#'
#' @param request a list of request(s) created with `build_era5_request`.
#' @param uid character vector containing your CDS user ID.
#' @param out_path character vector of the location at which to download nc files.
#'
#' @export
#'
#'
request_era5 <- function(request, uid, out_path) {

  for(req in 1:length(request)) {

    ecmwfr::wf_request(user = as.character(uid),
                       request = request[[req]],
                       transfer = TRUE,
                       path = out_path,
                       verbose = TRUE,
                       time_out = 18000)
  }
}
