#' read genome files (i.e., UCSC "chromSize" files)
#' 
#' Genome files contain chromosome name and size information. These sizes are 
#' used by downstream functions to identify intervals that have ends outside the
#' defined size.
#' 
#' @param filename file containing chrom/contig names and sizes, one-pair-per-line, 
#'   tab-delimited
#'   
#' @return \code{dplyr::tbl_df} with colnames \code{chrom} and \code{size},
#'   sorted by \code{chrom}
#'   
#' @examples
#' genome <- read_genome('hg19.chrom.sizes.gz')  
#' 
#' @export
read_genome <- function(filename) {
  colnames <- c('chrom', 'size')
  genome <- readr::read_tsv(filename, col_names = colnames)
  genome <- dplyr::tbl_df(genome) %>% dplyr::arrange(chrom)
  genome
}