#' Convert BED12 to individual exons in BED6.
#'
#' After conversion to BED6 format, the `score` column contains the exon
#' number, with respect to strand (i.e., the first exon for `-` strand
#' genes will have larger start and end coordinates).
#'
#' @param x [tbl_interval()]
#' @family utilities
#'
#' @examples
#' x <- read_bed12(valr_example('mm9.refGene.bed.gz'))
#' bed12_to_exons(x)
#'
#' @export
bed12_to_exons <- function(x) {

  if (!is.tbl_interval(x)) x <- tbl_interval(x)

  if (! ncol(x) == 12)
    stop("expected 12 column input", call. = FALSE)

  res <- bed12toexons_impl(x)
  res <- tibble::as_tibble(res)
  res <- arrange(res, chrom, start)

  res
}
