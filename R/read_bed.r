#' @title Read BED and related files
#'   
#' @description \code{read_bed} reads BED files and \code{read_bedgraph} reads
#'   bedGraph files. Return value are \code{data.frame} that is sorted by chrom
#'   and start unless otherwise specified.
#'   
#' @name read_bed
#'   
#' @param filename file or URL
#' @param n_fields number fields in the BED file
#' @param col_types column type spec for \code{readr::read_tsv}
#' @param sort sort the tbl by chrom and start
#' @param factor_cols factor the \code{chrom} and \code{strand} columns
#' @param ... options to pass to \code{readr::read_tsv}
#'   
#' @return \code{data.frame}

#' @details \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}
#'   
#' @examples
#' 
#' # \code{read_bed} assumes 3 field BED format. 
#' bed3_tbl <- read_bed('3fields.bed.gz')
#' 
#' bed6_tbl <- read_bed('6fields.bed.gz', n_fields = 6)
#' 
#' # read BED12 format
#' bed12_tbl <- read_bed12('12fields.bed.gz')
#' 
#' # Returned \code{data.frame} is sorted by \code{chrom}, \code{start} and \code{end} unless specified
#' sorted_bed_tbl <- read_bed('unsorted.bed.gz')
#' unsorted_bed_tbl <- read_bed('unsorted.bed.gz', sort = FALSE)
#' 
#' # \code{chrom} and \code{strand} are converted to factors unless specified
#' bed_tbl <- read_bed('3fields.bed.gz', factor_cols = FALSE)
#'   
#' @export
read_bed <- function(filename, n_fields = 3, col_types = bed12_coltypes,
                     sort = TRUE, factor_cols = TRUE, ...) {
  
  coltypes <- col_types[1:n_fields]
  colnames <- names(coltypes)
  
  bed_raw <- read_tsv(filename, col_names = colnames, col_types = coltypes, ...)
  bed_tbl <- tbl_df(bed_raw) 

  # factorize chrom and strand
  if (factor_cols) {
    bed_tbl <- bed_tbl %>% mutate(chrom = as.factor(chrom))
    
    if ('strand' %in% colnames(bed_tbl)) {
      bed_tbl <- bed_tbl %>% mutate(strand = as.factor(strand))
    }
  } 
  
  if (sort) {
    bed_tbl <- bed_sort(bed_tbl)
  }
 
  bed_tbl
}

#' @rdname read_bed
#' 
#' @param filename file name or URL
#' @param ... params to pass to \code{read_bed}
#' 
#' @details https://genome.ucsc.edu/FAQ/FAQformat.html#format1
#' 
#' @examples
#' 
#' bed12_tbl <- read_bed12('mm9.bed12.gz')
#' 
#' @export
read_bed12 <- function(filename, ...) {
  bed12_tbl <- read_bed(filename, n_fields = 12)
  bed12_tbl
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
  # load as bed4, rename `value` column and covert to double
  bedgraph_tbl <- read_bed(filename, n_fields = 4) %>%
    rename(value = name) %>%
    mutate(value = as.double(value))
  bedgraph_tbl
}

# thickStart renamed to cds_start
# thickEnd renamed to cds_end
bed12_coltypes <- list(chrom = col_character(),
                       start = col_double(),
                       end = col_double(),
                       name = col_character(),
                       score = col_character(),
                       strand = col_character(),
                       cds_start = col_double(),
                       cds_end = col_double(),
                       item_rgb = col_character(),
                       exon_count = col_double(),
                       exon_sizes = col_character(),
                       exon_starts = col_character())
