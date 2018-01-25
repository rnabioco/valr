#' Read genome files.
#'
#' Genome files (UCSC "chromSize" files) contain chromosome name and size
#' information. These sizes are used by downstream functions to identify
#' computed intervals that have coordinates outside of the genome bounds.
#'
#' @param path containing chrom/contig names and sizes, one-pair-per-line,
#'   tab-delimited
#'
#' @return [tbl_genome()], sorted by `size`
#'
#' @note URLs to genome files can also be used.
#'
#' @family read functions
#'
#' @examples
#' read_genome(valr_example('hg19.chrom.sizes.gz'))
#'
#' \dontrun{
#' # `read_genome` accepts a URL
#' read_genome('https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes')
#' }
#'
#' @export
read_genome <- function(path) {
  colnames <- c("chrom", "size")
  genome <- suppressMessages(readr::read_tsv(path, col_names = colnames))
  genome <- arrange(genome, desc(size))
  genome <- as.tbl_genome(genome)
  genome
}

#' Select intervals bounded by a genome.
#'
#' Used to remove out-of-bounds intervals, or trim interval coordinates using a
#' `genome`.
#'
#' @param x [tbl_interval()]
#' @param genome [tbl_genome()]
#' @param trim adjust coordinates for out-of-bounds intervals
#'
#' @return [tbl_interval()]
#'
#' @family utilities
#'
#' @examples
#' x <- trbl_interval(
#'  ~chrom, ~start, ~end,
#'  "chr1", -100,   500,
#'  "chr1", 100,    1e9,
#'  "chr1", 500,    1000
#' )
#'
#' genome <- read_genome(valr_example('hg19.chrom.sizes.gz'))
#'
#' # out-of-bounds are removed by default ...
#' bound_intervals(x, genome)
#'
#' # ... or can be trimmed within the bounds of a genome
#' bound_intervals(x, genome, trim = TRUE)
#'
#' @export
bound_intervals <- function(x, genome, trim = FALSE) {
  if (!is.tbl_interval(x)) x <- as.tbl_interval(x)

  res <- left_join(x, genome, by = "chrom")
  if (trim) {
    res <- mutate(
      res,
      start = ifelse(start < 0,
        0,
        pmin(start, size - 1)
      ),
      end = ifelse(end > size,
        size,
        pmax(1, end)
      )
    )
    res <- select(res, -size)
  } else {
    res <- filter(res, start >= 0 & start < size & end <= size & end > 0)
    res <- select(res, -size)
  }

  if (any(res$start == res$end)) {
    warning(paste0(
      sum(res$start == res$end),
      " intervals generated with same start and end were discarded"
    ))
  }

  res <- res[res$start != res$end, ]

  res
}
