#' Convert BED12 to individual exons in BED6.
#'
#' After conversion to BED6 format, the `score` column contains the exon
#' number, with respect to strand (i.e., the first exon for `-` strand
#' genes will have larger start and end coordinates).
#'
#' @param x [ivl_df]
#'
#' @family utilities
#'
#' @examples
#' x <- read_bed12(valr_example('mm9.refGene.bed.gz'))
#'
#' bed12_to_exons(x)
#'
#' @export
bed12_to_exons <- function(x) {
  x <- check_interval(x)

  if (!ncol(x) == 12) {
    stop("expected 12 column input", call. = FALSE)
  }

  res <- bed12toexons_impl(x)
  res <- res[, c("chrom", "start", "end", "name", "score", "strand")]
  res <- bed_sort(res)

  res
}
