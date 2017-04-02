#' Tibble for intervals.
#'
#' Inherits from \code{\link[tibble]{tibble}} and enforces \code{chrom},
#' \code{start} and \code{end} columns.
#'
#' @param ... params for \code{\link[tibble]{tibble}}
#' @param validate check valid column names
#'
#' @export
tbl_interval <- function(..., validate = TRUE) {
  out <- tibble::tibble(...)
  if (validate) {
    out <- check_interval(out)
  }
  class(out) <- union('tbl_ivl', class(out))
  out
}

#' @rdname tbl_interval
#' @export
tbl_ivl <- tbl_interval

#' @export
as_tbl_ivl <- function (x, ...) {
  UseMethod("as_tbl_ivl")
}

#' @export
as_tbl_ivl.tbl_df <- function(x, ...) {
  class(x) <- union('tbl_ivl', class(x))
  x
}

#' Tibble for reference sizes.
#'
#' Inherits from \code{\link[tibble]{tibble}} and enforces \code{chrom} and
#' \code{size} columns. Identical to UCSC "chromSizes" files.
#'
#' @param ... params for \code{\link[tibble]{tibble}}
#' @param validate check valid column names
#'
#' @export
tbl_sizes <- function(..., validate = TRUE) {
  out <- tibble::tibble(...)
  if (validate) {
    out <- check_sizes(out)
  }
  class(out) <- union('tbl_szs', class(out))
  out
}

#' @rdname tbl_sizes
#' @export
tbl_szs <- tbl_sizes

#' @export
as_tbl_szs <- function (x, ...) {
  UseMethod("as_tbl_szs")
}

#' @export
as_tbl_szs.tbl_df <- function(x, ...) {
  class(x) <- union('tbl_szs', class(x))
  x
}

# Validity checks ---------------------------------------------------

check_interval <- function(x) {
  expect_names <- c('chrom', 'start', 'end')
  x
}

check_sizes <- function(x) {
  expect_names <- c('chrom', 'size')
  x
}
