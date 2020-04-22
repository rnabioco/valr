#' Convert Granges to bed tibble
#'
#' @param x GRanges object to convert to bed tibble.
#'
#' @return [tibble::tibble()]
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
#' gr_to_bed(gr)
#'
#' # There are two ways to convert a bed-like data.frame to GRanges:
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
