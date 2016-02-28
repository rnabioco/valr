#' Merge overlapping intervals
#' 
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/merge.html}
#' 
#' @examples 
#' bed_tbl <- dplyr::tibble(
#'  ~chrom, ~start, ~end,
#'  "chr1", 100,    200,
#'  "chr1", 150,    250,
#'  "chr2", 200,    400
#' )
#' 
#' bed_tbl <- bed_sort(bed_tbl) 
#' bed_merge(bed_tbl)
#' 
#' @export
bed_merge <- function(bed_tbl) {
  
  if (!attr(bed_tbl, "sorted")) {
    stop("bed_merge expects sorted tbl, use bed_sort")
  }
  
  res <- merge_cpp(bed_tbl)
  res <- bed_sort(res)

  res  
}