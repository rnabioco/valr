#' Convert BED12 to individual exons in BED6.
#' 
#' After conversion to BED6 format, the \code{score} column contains the exon 
#' number, with respect to strand (i.e., the first exon for \code{-} strand
#' genes will have larger start and end coordinates).
#' 
#' @param x tbl in BED12 format
#' @family utils
#'   
#' @examples
#' bed12_path <- valr_example('mm9.bed12.gz')
#' x <- read_bed12(bed12_path)
#' bed12_to_exons(x) 
#' 
#' @export
bed12_to_exons <- function(x) {
 
  if (! ncol(x) == 12) 
    stop('expected 12 column input', call. = FALSE)
  
  res <- tidyr::unnest(x, .exon_size = stringr::str_split(stringr::str_replace(exon_sizes, ',$', ''), ','),
                          .exon_start = stringr::str_split(stringr::str_replace(exon_starts, ',$', ''), ',')) 
  res <- mutate(res,
                .exon_size = as.double(.exon_size),
                .exon_start = as.double(.exon_start))
  res <- group_by(res, name)
  res <- mutate(res,
                .exon_num = ifelse(strand == '+',
                                          row_number(),
                                          rev(row_number())))
  res <- mutate(res,
                .start = start + .exon_start,
                .end = .start + .exon_size,
                .score = .exon_num)
  res <- select(res, chrom, .start, .end, name, .exon_num, strand)
  res <- rename(res, start = .start, end = .end, score = .exon_num)
  
  res <- ungroup(res)
  res <- bed_sort(res)
  
  res
}
