#' Tibble for intervals.
#'
#' Required column names are \code{chrom}, \code{start} and \code{end}
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

#' Construct an tbl_interval using tribble formatting.
#'
#' @param ... data for \code{\link[tibble]{tribble}}
#'
#' @return \code{\link{tbl_interval}}
#'
#' @examples
#' trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1',  1,      50,
#'   'chr1',  10,     75
#' )
#'
#' @export
trbl_interval <- function(...) {

  out <- tibble::tribble(...)
  out <- tbl_interval(out)
  out
}

#' Test if the object is a tbl_interval.
#'
#' @param x An object
#' @return \code{TRUE} if the object inherits from the \code{\link{tbl_interval}} class.
#' @export
is.tbl_interval <- function(x) {
  "tbl_ivl" %in% class(x)
}

#' Tibble for reference sizes.
#'
#' Equivalent to information in UCSC "chromSizes" files. Required column names are:
#' \code{chrom} and \code{size}
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
#' is.tbl_genome(genome)
#' genome <- tbl_genome(genome)
#' is.tbl_genome(genome)
#'
#' @export
tbl_genome <- function(x, ..., .validate = TRUE) {
  out <- tibble::as_tibble(x, ...)
  if (.validate) {
    out <- check_genome(out)
  }
  class(out) <- union('tbl_gnm', class(out))
  out
}

#' Construct a tbl_genome using tribble formatting.
#'
#' @param ... for \code{\link[tibble]{tribble}}
#'
#' @return \code{\link{tbl_genome}}
#'
#' @examples
#' trbl_genome(
#'   ~chrom, ~size,
#'   'chr1', 1e6
#' )
#'
#' @export
trbl_genome <- function(...) {

  out <- tibble::tribble(...)
  out <- tbl_genome(out)
  out
}

#' Test if the object is a tbl_genome.
#'
#' @param x An object
#' @return \code{TRUE} if the object inherits from the \code{\link{tbl_genome}} class.
#' @export
is.tbl_genome <- function(x) {
  "tbl_gnm" %in% class(x)
}

# Validity checks ---------------------------------------------------

check_interval <- function(x) {
  expect_names <- c('chrom', 'start', 'end')
  check_names(x, expect_names)
  x
}

check_genome <- function(x) {
  expect_names <- c('chrom', 'size')
  check_names(x, expect_names)

  # check for unique refs
  chroms <- x[['chrom']]
  dups <- duplicated(chroms)

  if (any(dups)) {
    stop(sprintf("duplicate chroms in genome: %s",
                 paste0(chroms[dups], collapse = ', ')))
  }
  x
}

check_names <- function(x, expected) {
  missing <- setdiff(expected, names(x))
  if (length(missing) != 0) {
    stop(sprintf("expected %d required names, missing: %s",
                 length(expected),
                 paste0(missing, collapse = ', ')))
  }
}
