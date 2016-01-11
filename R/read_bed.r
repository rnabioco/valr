#' Methods to read BED and related files
#' 
#' @param filename file or URL
#' @param colnames override standard names with vector of names
#' @param ... options to pass to readr::read_tsv
#' 
read_bed <- function(filename, colnames = bed12_colnames, ...) {
  bed_raw <- readr::read_tsv(filename, col_names = colnames)
  bed_tbl <- dplyr::tbl_df(bed_raw)
  bed_tbl
}

read_bedgraph <- function(filename, colnames = bedgraph_colnames, ...) {
  bedgraph_raw <- readr::read_tsv(filename, col_names = colnames)
  bedgraph_tbl <- dplyr::tbl_df(bedgraph_raw)
  bedgraph_tbl
}

bed12_colnames <- c('chrom', 'start', 'stop', 'name',
                    'score', 'strand', 'cds_start', 'cds_end',
                    'item_rgb', 'exon_count', 'exon_sizes', 'exon_starts')

bedgraph_colnames <- c('chrom', 'start', 'stop', 'count')
