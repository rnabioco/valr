#' @param min_overlap minimum overlap in base pairs required for the operation.
#'   Set to `1` to exclude book-ended intervals (matching bedtools behavior), or
#'   `0` to include them (legacy valr behavior). The default will change from
#'   `0` to `1` in a future version.
