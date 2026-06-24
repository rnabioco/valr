#' @param min_overlap minimum overlap in base pairs required for the operation.
#'   Defaults to `1`, which excludes book-ended intervals (those that touch but
#'   do not overlap), matching bedtools behavior. Set to `0` to include
#'   book-ended intervals (the legacy valr behavior).
