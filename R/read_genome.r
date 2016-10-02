#' read genome files (i.e., UCSC "chromSize" files)
#' 
#' Genome files contain chromosome name and size information. These sizes are 
#' used by downstream functions to identify computed intervals that have ends
#' outside the defined size.
#' 
#' @param filename file containing chrom/contig names and sizes,
#'   one-pair-per-line, tab-delimited
#'   
#' @return \code{data_frame} with colnames \code{chrom} and \code{size}, sorted
#'   by \code{size}
#'   
#' @family read functions
#' 
#' @examples
#' genome <- valr_example('hg19.chrom.sizes.gz')
#' read_genome(genome)
#' 
#' \dontrun{
#' # can also read from a URL
#' read_genome('https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes')
#' }
#' 
#' @export
read_genome <- function(filename) {
  colnames <- c('chrom', 'size')
  genome <- suppressMessages(readr::read_tsv(filename, col_names = colnames))
  genome <- arrange(genome, desc(size))
  genome
}

#' Select intervals bounded by a genome.
#' 
#' @param x a tbl of intervals
#' @param genome a tbl of chrom sizes
#' @param trim adjust coordinates for out-of-bounds intervals
#'
#' @return \code{data_frame}
#' 
#' @rdname read_genome
#'  
#' @examples
#' x <- tibble::tribble(
#'  ~chrom, ~start, ~end,
#'  "chr1", -100,   500,
#'  "chr1", 100,    1e9,
#'  "chr1", 500,    1000
#' )
#' 
#' genome <- read_genome(valr_example('hg19.chrom.sizes.gz'))
#' 
#' bound_intervals(x, genome)
#' bound_intervals(x, genome, trim = TRUE)
#' 
#' @export
bound_intervals <- function(x, genome, trim = FALSE) {
  res <- left_join(x, genome, by = "chrom") 
  if (trim) {
    res <- mutate(res,
                  start = ifelse(start < 1, 1, start),
                  end = ifelse(end > size, size, end))
    res <- select(res, -size)
  } else {
    res <- filter(res, start >= 1 & end <= size)
    res <- select(res, -size)
  }
   
  res
}
