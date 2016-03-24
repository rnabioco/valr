#' read genome files (i.e., UCSC "chromSize" files)
#' 
#' Genome files contain chromosome name and size information. These sizes are 
#' used by downstream functions to identify computed intervals that have ends
#' outside the defined size.
#' 
#' @param filename file containing chrom/contig names and sizes,
#'   one-pair-per-line, tab-delimited
#'   
#' @return \code{data.frame} with colnames \code{chrom} and \code{size}, sorted
#'   by \code{chrom}
#'   
#' @examples
#' read_genome('hg19.chrom.sizes.gz')  
#' @export
read_genome <- function(filename) {
  colnames <- c('chrom', 'size')
  genome <- read_tsv(filename, col_names = colnames)
  genome <- tbl_df(genome) %>% arrange(chrom)
  genome
}

#' select bounded intervals
#' 
#' @param bed_tbl a tbl of intervals
#' @param genome a tbl of chrom sizes
#' 
#' @export
bound_intervals <- function(bed_tbl, genome) {
   res <-
     bed_tbl %>% 
     left_join(genome, by = 'chrom') %>%
     filter(start >= 1 & end <= size) %>%
     select(-size)
   
   res
}
