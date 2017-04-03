#' Tibble for intervals.
#'
#' @details Required column names:
#'
#' \itemize{
#'   \item{\code{chrom}}
#'   \item{\code{start}}
#'   \item{\code{end}}}
#'
#' @param x A \code{data_frame}
#' @param ... params for \code{\link[tibble]{tibble}}
#' @param .validate check valid column names
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1',  1,      50,
#'   'chr1',  10,     75,
#'   'chr1',  100,    120
#' )
#'
#' is.tbl_interval(x)
#' x <- tbl_interval(x)
#' is.tbl_interval(x)
#'
#' @export
tbl_interval <- function(x, ..., .validate = TRUE) {
  out <- tibble::as_tibble(x, ...)
  if (.validate) {
    out <- check_interval(out)
  }
  class(out) <- union('tbl_ivl', class(out))
  out
}

#' Test if the object is a tbl_interval.
#'
#' @param x An object
#' @return \code{TRUE} if the object inherits from the \code{tbl_interval} class.
#' @export
is.tbl_interval <- function(x) {
  "tbl_ivl" %in% class(x)
}

#' Tibble for reference sizes.
#'
#' Equivalent to information in UCSC "chromSizes" files.
#'
#' @details Required column names:
#'
#' \itemize{
#'   \item{\code{chrom}}
#'   \item{\code{size}}}
#'
#' @param x A \code{data_frame}
#' @param ... params for \code{\link[tibble]{tibble}}
#' @param .validate check valid column names
#'
#' @examples
#' genome <- tibble::tribble(
#'   ~chrom, ~size,
#'   'chr1', 1e6,
#'   'chr2', 1e7
#' )
#'
#' is.tbl_sizes(genome)
#' genome <- tbl_sizes(genome)
#' is.tbl_sizes(genome)
#'
#' @export
tbl_sizes <- function(x, ..., .validate = TRUE) {
  out <- tibble::as_tibble(x, ...)
  if (.validate) {
    out <- check_sizes(out)
  }
  class(out) <- union('tbl_szs', class(out))
  out
}

#' Test if the object is a tbl_sizes.
#'
#' @param x An object
#' @return \code{TRUE} if the object inherits from the \code{tbl_sizes} class.
#' @export
is.tbl_sizes <- function(x) {
  "tbl_szs" %in% class(x)
}

# Validity checks ---------------------------------------------------

check_interval <- function(x) {
  expect_names <- c('chrom', 'start', 'end')
  check_names(x, expect_names)
  x
}

check_sizes <- function(x) {
  expect_names <- c('chrom', 'size')
  check_names(x, expect_names)

  # check for unique refs
  chroms <- x[['chrom']]
  dups <- duplicated(chroms)

  if (any(dups)) {
    chrom_dups <- chroms[dups]
    if (length(chrom_dups) > 1) {
      err <- paste0(chrom_dups, collapse = ', ')
    } else {
      err <- chrom_dups[1]
    }
    stop(sprintf("duplicate chroms in genome: %s", err))
  }

  x
}

check_names <- function(x, expected) {

  names_missing <- setdiff(expected, names(x))

  if (length(names_missing != 0)) {
    if (length(names_missing) > 1) {
      err <- paste0(names_missing, collapse = ', ')
    } else {
      err <- names_missing[1]
    }

    stop(sprintf("expected %d required names, missing: %s",
                 length(expected), err))
  }

}
