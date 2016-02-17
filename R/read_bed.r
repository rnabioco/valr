#' @title Read BED and related files
#'
#' @description \code{read_bed} reads BED files and \code{read_bedgraph} reads bedGraph files.
#' Return values are \code{dplyr::tbl_df}s that are sorted by chrom and start unless otherwise specified.
#'  
#' @name read_bed
#'   
#' @param filename file or URL
#' @param n_fields number fields in the BED file
#' @param col_names column names to use for \code{dplyr::tbl_df}
#' @param sort sort the tbl by chrom and start
#' @param ... options to pass to \code{readr::read_tsv}
<<<<<<< Updated upstream
#' 
#' @details https://genome.ucsc.edu/FAQ/FAQformat.html#format1
=======
#'   
#' @details \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}
>>>>>>> Stashed changes
#'   
#' @examples
#' 
#' # \code{read_bed} assumes 3 field BED format. 
#' bed3_tbl <- read_bed('3fields.bed.gz')
#' bed6_tbl <- read_bed('6fields.bed.gz', n_fields = 6)
#' bed12_tbl <- read_bed('12fields.bed.gz', n_fields = 12)
#' 
#' # \code{readr} will generate warnings if extra columns are identified in the ouput.
#' Empty columns are dropped (TODO).
#' 
#' # Returned \code{tbl_df} is sorted by chrom and start unless specified
#' sorted_bed_tbl <- read_bed('unsorted.bed.gz')
#' unsorted_bed_tbl <- read_bed('unsorted.bed.gz', sort = FALSE)
#' 
#' # `chrom` and `strand` are converted to factors unless specified
#' bed_tbl <- read_bed('3fields.bed.gz')
#' bed_tbl <- read_bed('3fields.bed.gz', factor_cols = FALSE)
#' 
#' @return \code{dplyr::tbl_df}
#'   
#' @export
read_bed <- function(filename, n_fields = 3, col_names = bed12_colnames,
                     sort = TRUE, factor_cols = TRUE, ...) {
  
  colnames <- col_names[1:n_fields]
  
  bed_raw <- readr::read_tsv(filename, col_names = colnames, ...)
  bed_tbl <- dplyr::tbl_df(bed_raw) 
  
  if (sort) {
    bed_tbl <- bed_tbl %>% dplyr::arrange(chrom, start)
  }

  # factorize chrom and strand
  if (factor_cols) {
    
    bed_tbl <- bed_tbl %>% dplyr::mutate(chrom = as.factor(chrom))
    
    if ('strand' %in% colnames(bed_tbl)) {
      bed_tbl <- bed_tbl %>% dplyr::mutate(strand = as.factor(strand))
    }
  } 
  
  bed_tbl
}

#' @rdname read_bed
#' 
#' @param filename file name or URL
#' @param ... params to pass to \code{read_bed}
#' 
#' @details https://genome.ucsc.edu/goldenPath/help/bedgraph.html
#' 
#' @examples
#' 
#' bedgraph_tbl <- read_bedgraph('test.bg.gz')
#' 
#' @export
read_bedgraph <- function(filename, ...) {
  bedgraph_tbl <- read_bed(filename, n_fields = 4,
                           col_names = bedgraph_colnames, ...)
  bedgraph_tbl
}

# thickStart renamed to cds_start
# thickEnd renamed to cds_end
bed12_colnames <- c('chrom', 'start', 'end', 'name',
                    'score', 'strand', 'cds_start', 'cds_end',
                    'item_rgb', 'exon_count', 'exon_sizes', 'exon_starts')

bedgraph_colnames <- c('chrom', 'start', 'end', 'value')

