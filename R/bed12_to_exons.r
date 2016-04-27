#' Convert BED12 to individual exons in BED6.
#' 
#' After conversion to BED6 format, the \code{score} column contains the exon 
#' number, with respect to strand (i.e., exon 1 for `-` strand genes will have 
#' larger start and end coordinates).
#' 
#' @param bed12_tbl tbl in BED12 format
#'   
#' @examples
#' bed12_path <- system.file('extdata', 'mm9.bed12.gz', package = 'Rbedtools')
#' bed12_tbl <- read_bed12(bed12_path)
#' bed6_exon_tbl <- bed12_to_exons(bed12_tbl) 
#' 
#' # first exons
#' subset(bed6_exon_tbl, score == 1)
#' 
#' @export
bed12_to_exons <- function(bed12_tbl) {
  
  assert_that(ncol(bed12_tbl) == 12)
  
  res <- bed12_tbl %>%
    tidyr::unnest(.exon_size = str_split(str_replace(exon_sizes, ',$', ''), ','),
                  .exon_start = str_split(str_replace(exon_starts, ',$', ''), ',')) %>%
    mutate(.exon_size = as.double(.exon_size),
           .exon_start = as.double(.exon_start)) %>%
    group_by(name) %>%
    mutate(.exon_num = ifelse(strand == '+',
                              row_number(),
                              rev(row_number()))) %>%
    mutate(.start = start + .exon_start,
           .end = .start + .exon_size,
           .score = .exon_num) %>%
    select(chrom, .start, .end, name, .exon_num, strand) %>%
    rename(start = .start, end = .end, score = .exon_num) %>%
    ungroup() %>%
    bed_sort()
  
  res
}
