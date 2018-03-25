#' Tibble for intervals.
#'
#' Required column names are `chrom`, `start` and `end`.
#'
#' @param x A `data_frame`
#' @param ... params for [tibble::tibble()]
#' @param .validate check valid column names
#'
#' @rdname tbl_interval
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1',  1,     50,
#'   'chr1',  10,    75,
#'   'chr1',  100,   120
#' )
#'
#' is.tbl_interval(x)
#'
#' x <- tbl_interval(x)
#' is.tbl_interval(x)
#'
#' @export
tbl_interval <- function(x, ..., .validate = TRUE) {
  out <- tibble::as_tibble(x, ...)
  if (.validate) {
    out <- check_interval(out)
  }
  class(out) <- union("tbl_ivl", class(out))
  out
}

#' Coerce objects to tbl_intervals.
#'
#' This is an S3 generic. valr includes methods to coerce tbl_df and GRanges
#' objects.
#'
#' @param x object to convert to tbl_interval.
#'
#' @return [tbl_interval()]
#'
#' @examples
#' \dontrun{
#' gr <- GenomicRanges::GRanges(
#'         seqnames = S4Vectors::Rle(
#'                      c("chr1", "chr2", "chr1", "chr3"),
#'                      c(1, 1, 1, 1)),
#'         ranges   = IRanges::IRanges(
#'                      start = c(1, 10, 50, 100),
#'                      end = c(100, 500, 1000, 2000),
#'                      names = head(letters, 4)),
#'         strand   = S4Vectors::Rle(
#'                      c("-", "+"), c(2, 2))
#'       )
#'
#' as.tbl_interval(gr)
#'
#' # There are two ways to convert a tbl_interval to GRanges:
#'
#' gr <- GenomicRanges::GRanges(
#'         seqnames = S4Vectors::Rle(x$chrom),
#'         ranges   = IRanges::IRanges(
#'                      start = x$start + 1,
#'                      end = x$end,
#'                      names = x$name),
#'         strand   = S4Vectors::Rle(x$strand)
#'         )
#' # or:
#'
#' gr <- GenomicRanges::makeGRangesFromDataFrame(dplyr::mutate(x, start = start +1))
#'
#' }
#'
#' @export
as.tbl_interval <- function(x) {
  UseMethod("as.tbl_interval")
}

#' @export
#' @rdname as.tbl_interval
as.tbl_interval.tbl_df <- function(x) {
  tbl_interval(x)
}

#' @export
#' @rdname as.tbl_interval
as.tbl_interval.data.frame <- function(x) {
  tbl_interval(x)
}

#' @export
#' @rdname as.tbl_interval
as.tbl_interval.GRanges <- function(x) {
  # https://www.biostars.org/p/89341/
  res <- tibble(
    chrom = as.character(x@seqnames),
    start = x@ranges@start - 1,
    end = x@ranges@start - 1 + x@ranges@width,
    name = rep(".", length(x)),
    score = rep(".", length(x)),
    strand = as.character(x@strand)
  )

  res <- mutate(res, strand = ifelse(strand == "*", ".", strand))
  tbl_interval(res)
}

#' Construct a tbl_interval using tribble formatting.
#'
#' @rdname tbl_interval
#'
#' @return [tbl_interval()]
#
#' @export
trbl_interval <- function(...) {
  out <- tibble::tribble(...)
  out <- as.tbl_interval(out)
  out
}

#' Test if the object is a tbl_interval.
#'
#' @param x An object
#' @return `TRUE` if the object inherits from the [tbl_interval()] class.
#' @export
is.tbl_interval <- function(x) {
  "tbl_ivl" %in% class(x)
}

#' Tibble for reference sizes.
#'
#' Equivalent to information in UCSC "chromSizes" files. Required column names are:
#' `chrom` and `size`
#'
#' @param x A `data_frame`
#' @param ... params for [tibble::tibble()]
#' @param .validate check valid column names
#'
#' @rdname tbl_genome
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
  class(out) <- union("tbl_gnm", class(out))
  out
}

#' Coerce objects to tbl_genome.
#'
#' This is an S3 generic. valr includes methods to coerce tbl_df and data.frame
#' objects.
#'
#' @param x object to convert to tbl_genome.
#'
#' @return [tbl_genome()]
#'
#' @export
as.tbl_genome <- function(x) {
  UseMethod("as.tbl_genome")
}

#' @export
#' @rdname as.tbl_genome
as.tbl_genome.tbl_df <- function(x) {
  tbl_genome(x)
}

#' @export
#' @rdname as.tbl_genome
as.tbl_genome.data.frame <- function(x) {
  tbl_genome(x)
}

#' Construct a tbl_genome using tribble formatting.
#'
#' @return [tbl_genome()]
#'
#' @rdname tbl_genome
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
#' @return `TRUE` if the object inherits from the [tbl_genome()] class.
#' @export
is.tbl_genome <- function(x) {
  "tbl_gnm" %in% class(x)
}

# Validity checks ---------------------------------------------------

check_interval <- function(x) {
  expect_names <- c("chrom", "start", "end")
  check_names(x, expect_names)
  x
}

check_genome <- function(x) {
  expect_names <- c("chrom", "size")
  check_names(x, expect_names)

  # check for unique refs
  chroms <- x[["chrom"]]
  dups <- duplicated(chroms)

  if (any(dups)) {
    stop(sprintf(
      "duplicate chroms in genome: %s",
      paste0(chroms[dups], collapse = ", ")
    ))
  }
  x
}

check_names <- function(x, expected) {
  missing <- setdiff(expected, names(x))
  if (length(missing) != 0) {
    stop(sprintf(
      "expected %d required names, missing: %s",
      length(expected),
      paste0(missing, collapse = ", ")
    ))
  }
}
