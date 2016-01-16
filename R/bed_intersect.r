#' instersect BED intervals
#' 
#' @param df_a BED intervals 
#' @param df_b BED intervals 
#' 
#' @export
bedtools_intersect <- function(df_a, df_b, by_chrom = TRUE) {
  
  res <- df_a %>%
    group_by(chrom) %>%
    intersect_cpp(., df_b)
  res <- intersect_cpp(df_a, df_b)
  res
}