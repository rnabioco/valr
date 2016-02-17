#' Merge overlapping intervals
#' 
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/merge.html}
#' 
#' @examples 
#' interval_tbl <- dplyr::tibble(
#'  ~chrom, ~start, ~end,
#'  "chr1", 100,    200,
#'  "chr1", 150,    250,
#'  "chr2", 200,    400
#' )
#' 
#' bed_merge(interval_tbl)
#' 
#' @export
bed_merge <- function(interval_tbl) {
  merge_result <- merge_cpp(interval_tbl)
  merge_result
}