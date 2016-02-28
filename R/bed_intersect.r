#' instersect BED intervals
#' 
#' @param df_a BED intervals 
#' @param df_b BED intervals 
#' 
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html}
#'  
#' @export
bed_intersect <- function(df_a, df_b, by_chrom = TRUE) {
 
  if (!attr(df_a, "sorted") || !attr(df_b, "sorted")) {
    stop("bed_intersection requires sorted intervals, use bed_sort") 
  }
  
  res <- df_a %>%
    group_by(chrom) %>%
    intersect_cpp(., df_b)
  
  res <- intersect_cpp(df_a, df_b)
  
  res
}
