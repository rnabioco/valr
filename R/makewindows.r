#'
#' makewindows.R
#' 
#' @description dplyr implementation of BEDtools makewindows
#' 
#' @param bed_df BED data in \code{dplyr::tbl_df} format
#' @param reverse reverse window numbers
#' 
#' @return \code{dplyr::tbl_df} with `win.num` column
#' 
bedtools_makewindows <- function(bed_df, reverse = FALSE) {
  
  makewindows_result <- bed_df %>%
    group_by(chrom, start, end) %>%
    mutate(win.num = cut(end - start, size))
    
  makewindows_result 
}