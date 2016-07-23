#' @title Read BED and related files.
#'   
#' @description \code{read_bed} reads BED files and \code{read_bedgraph} reads 
#'   bedGraph files. Return value are \code{data.frame} that is sorted by chrom 
#'   and start unless otherwise specified.
#'   
#' @param filename file or URL
#' @param n_fields number fields in the BED file
#' @param col_types column type spec for \code{readr::read_tsv}
#' @param sort sort the tbl by chrom and start
#' @param factor_cols factor the \code{chrom} and \code{strand} columns
#' @param ... options to pass to \code{readr::read_tsv}
#'   
#' @return \code{data_frame}
#' 
#' @family read data
#'   
#' @details \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}
#'   
#' @examples
#' # read_bed assumes 3 field BED format. 
#' bed3_path <- system.file('extdata', '3fields.bed.gz', package = 'valr')
#' bed3_tbl <- read_bed(bed3_path)
#' 
#' bed6_path <- system.file('extdata', '6fields.bed.gz', package = 'valr')
#' bed6_tbl <- read_bed(bed6_path, n_fields = 6)
#' 
#' # Result is sorted by chrom and start unless `sort = FALSE`
#' unsorted_bed_tbl <- read_bed(bed3_path, sort = FALSE)
#' 
#' # use `is_sorted` 
#' is_sorted(unsorted_bed_tbl)
#' 
#' # chrom and strand are converted to factors unless specified
#' bed_tbl <- read_bed(bed3_path, factor_cols = FALSE)
#' is.factor(bed_tbl$chrom)
#' 
#' @export
read_bed <- function(filename, n_fields = 3, col_types = bed12_coltypes,
                     sort = TRUE, factor_cols = FALSE, ...) {
  
  coltypes <- col_types[1:n_fields]
  colnames <- names(coltypes)
  
  bed_raw <- read_tsv(filename, col_names = colnames, col_types = coltypes, ...)
  bed_tbl <- as_data_frame(bed_raw) 

  # factorize chrom and strand
  if (factor_cols) {
    bed_tbl <- mutate(bed_tbl, chrom = as.factor(chrom))
    
    if ('strand' %in% colnames(bed_tbl)) {
      bed_tbl <- mutate(bed_tbl, strand = as.factor(strand))
    }
  } 
  
  if (sort) {
    bed_tbl <- bed_sort(bed_tbl)
  }
 
  bed_tbl
}

#' @rdname read_bed
#' 
#' @details \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}
#' 
#' @family read data
#' 
#' @examples
#' bed12_path <- system.file('extdata', 'mm9.bed12.gz', package = 'valr')
#' bed12_tbl <- read_bed12(bed12_path)
#' 
#' @export
read_bed12 <- function(filename, ...) {
  bed12_tbl <- read_bed(filename, n_fields = 12)
  bed12_tbl
}

#' @rdname read_bed
#' 
#' @details \url{https://genome.ucsc.edu/goldenPath/help/bedgraph.html}
#' 
#' @family read data
#' 
#' @examples
#' bedgraph_path <- system.file('extdata', 'test.bg.gz', package = 'valr')
#' bedgraph_tbl <- read_bedgraph(bedgraph_path)
#' 
#' @export
read_bedgraph <- function(filename, ...) {
  # load as bed4, rename `value` column and covert to double
  bedgraph_tbl <- read_bed(filename, n_fields = 4, sort = FALSE) %>%
    rename(value = name) %>%
    mutate(value = as.double(value))
  bedgraph_tbl <- bed_sort(bedgraph_tbl)
  bedgraph_tbl
}

#' @rdname read_bed
bed12_coltypes <- list(chrom = col_character(),
                       start = col_integer(),
                       end = col_integer(),
                       name = col_character(),
                       score = col_character(),
                       strand = col_character(),
                       cds_start = col_integer(),
                       cds_end = col_integer(),
                       item_rgb = col_character(),
                       exon_count = col_integer(),
                       exon_sizes = col_character(),
                       exon_starts = col_character())
