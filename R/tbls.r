#' create an interval tibble
#'
#' Inherits from \code{[tibble](tibble)} and enforces \code{chrom}, \code{start}
#' and \code{end} columns.
#'
#' @inheritParams tibble::tibble
#'
#' @export
tbl_interval <- function(...) {
  out <- tibble::tibble(...)
  class(out) <- union('tbl_ivl', class(out))
  out
}

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

#' create a sizes tibble
#'
#' Inherits from \code{[tibble](tibble)} and enforces \code{chrom} and
#' \code{size} columns. Identical to UCSC "chromSizes" files.
#'
#' @inheritParams tibble::tibble
#'
#' @export
tbl_sizes <- function(...) {
  out <- tibble::tibble(...)
  class(out) <- union('tbl_szs', class(out))
  out
}

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
