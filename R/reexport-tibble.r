#' Deprecated tibble re-exports
#'
#' `r lifecycle::badge("deprecated")`
#'
#' `data_frame()` and `as_data_frame()` were re-exported from tibble, where they
#' have since been deprecated. Use [tibble::tibble()] and [tibble::as_tibble()]
#' instead.
#'
#' @param ... arguments passed on to [tibble::tibble()] or [tibble::as_tibble()]
#'
#' @keywords internal
#' @name reexports-deprecated
NULL

#' @rdname reexports-deprecated
#' @export
data_frame <- function(...) {
  lifecycle::deprecate_warn("0.9.2", "data_frame()", "tibble::tibble()")
  tibble::tibble(...)
}

#' @rdname reexports-deprecated
#' @export
as_data_frame <- function(...) {
  lifecycle::deprecate_warn("0.9.2", "as_data_frame()", "tibble::as_tibble()")
  tibble::as_tibble(...)
}

#' @importFrom tibble tribble
#' @export
tibble::tribble

#' @importFrom tibble tibble
#' @export
tibble::tibble

#' @importFrom tibble as_tibble
#' @export
tibble::as_tibble
