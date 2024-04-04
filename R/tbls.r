# Validity checks ---------------------------------------------------
#' Bed-like data.frame requirements for valr functions
#'
#' @name ivl_df
NULL

#' Bed-like data.frame requirements for valr functions
#' @rdname ivl_df
#' @name genome_df
NULL

#' Check bed-like data.frame requirements for valr compatibility
#'
#' Required column names for interval dataframes are
#' `chrom`, `start` and `end`. Internally interval dataframes are
#' validated using `check_interval()`
#'
#' @param x A `data.frame` or `tibble::tibble`
#' @rdname  ivl_df
#'
#' @examples
#' # using tibble
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 1, 50,
#'   "chr1", 10, 75,
#'   "chr1", 100, 120
#' )
#'
#' check_interval(x)
#'
#' # using base R data.frame
#' x <- data.frame(
#'   chrom = "chr1",
#'   start = 0,
#'   end = 100,
#'   stringsAsFactors = FALSE
#' )
#'
#' check_interval(x)
#'
#' @export
check_interval <- function(x) {
  expect_names <- c("chrom", "start", "end")
  check_names(x, expect_names)

  if (!tibble::is_tibble(x)) {
    x <- tibble::as_tibble(x)
  }
  x
}


#' Check genome file data.frame requirements for valr compatibility
#'
#' Required column names for genome dataframes are
#' `chrom` and `size`. Internally genome dataframes are
#' validated using `check_genome()`.
#'
#' @param x A `data.frame` or `tibble::tibble`
#' @rdname ivl_df
#'
#' @examples
#' # example genome input
#'
#' x <- tibble::tribble(
#'   ~chrom, ~size,
#'   "chr1", 1e6
#' )
#'
#' check_genome(x)
#'
#' @export
check_genome <- function(x) {
  expect_names <- c("chrom", "size")
  check_names(x, expect_names)

  # check for unique refs
  chroms <- x[["chrom"]]
  dups <- unique(chroms[duplicated(chroms)])

  if (length(dups) > 0) {
    cli::cli_abort(
      "duplicate chroms in genome: {dups}"
    )
  }

  if (!tibble::is_tibble(x)) {
    x <- tibble::as_tibble(x)
  }

  x
}

check_names <- function(x, expected) {
  missing <- setdiff(expected, names(x))
  if (length(missing) != 0) {
    n <- length(expected)
    cli::cli_abort(
      "expected {n} required names, missing: {missing}",
    )
  }
}


#' Convert Granges to bed tibble
#'
#' @param x GRanges object to convert to bed tibble.
#'
#' @return [tibble::tibble()]
#'
#' @examples
#' \dontrun{
#' gr <- GenomicRanges::GRanges(
#'   seqnames = S4Vectors::Rle(
#'     c("chr1", "chr2", "chr1", "chr3"),
#'     c(1, 1, 1, 1)
#'   ),
#'   ranges = IRanges::IRanges(
#'     start = c(1, 10, 50, 100),
#'     end = c(100, 500, 1000, 2000),
#'     names = head(letters, 4)
#'   ),
#'   strand = S4Vectors::Rle(
#'     c("-", "+"), c(2, 2)
#'   )
#' )
#'
#' gr_to_bed(gr)
#'
#' # There are two ways to convert a bed-like data.frame to GRanges:
#'
#' gr <- GenomicRanges::GRanges(
#'   seqnames = S4Vectors::Rle(x$chrom),
#'   ranges = IRanges::IRanges(
#'     start = x$start + 1,
#'     end = x$end,
#'     names = x$name
#'   ),
#'   strand = S4Vectors::Rle(x$strand)
#' )
#' # or:
#'
#' gr <- GenomicRanges::makeGRangesFromDataFrame(dplyr::mutate(x, start = start + 1))
#' }
#'
#' @export
gr_to_bed <- function(x) {
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
  res
}
