#' instersect BED intervals
#' 
#' @param bed_tbl_a BED intervals 
#' @param bed_tbl_b BED intervals 
#' 
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html}
#'  
#' @export
bed_intersect <- function(bed_tbl_a, bed_tbl_b) {

  if ( ! is_sorted(df_a) ) {
    df_a <- bed_sort(df_a)
  }
  if ( ! is_sorted(df_b) ) {
    df_b <- bed_sort(df_b)
  }
 
  res <- intersect_impl(df_a, df_b)
  
  res
}
