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
#' @examples
#' genome <- system.file('extdata', 'hg19.chrom.sizes.gz', package = 'Rbedtools')
#' read_genome(genome)
#' read_genome('https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes')
#' 
#' @export
read_genome <- function(filename) {
  colnames <- c('chrom', 'size')
  genome <-
    read_tsv(filename, col_names = colnames) %>%
    arrange(desc(size))
  genome
}

#' Select intervals bounded by a genome.
#' 
#' @param bed_tbl a tbl of intervals
#' @param genome a tbl of chrom sizes
#' @param trim adjust coordinates for out-of-bounds intervals
#'
#' @return \code{data_frame}
#'  
#' @examples
#' bed_tbl <- tibble::frame_data(
#'  ~chrom, ~start, ~end,
#'  "chr1", -100,   500,
#'  "chr1", 100,    1e9,
#'  "chr1", 500,    1000
#' )
#' 
#' genome <- read_genome('https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes')
#' 
#' bound_intervals(bed_tbl, genome)
#' bound_intervals(bed_tbl, genome, trim = TRUE)
#' 
#' @export
bound_intervals <- function(bed_tbl, genome, trim = FALSE) {
  if (trim) {
   res <-bed_tbl %>% 
    left_join(genome, by = "chrom") %>%
      mutate(start = ifelse(start < 1, 1, start),
             end = ifelse(end > size, size, end)) %>%
      select(-size)
  } else {
   res <-bed_tbl %>% 
     left_join(genome, by = 'chrom') %>%
     filter(start >= 1 & end <= size) %>%
     select(-size)
  }
   
   res
}
