#'
#' map.r
#' 
#' @description impelementation of BEDtools `map` function with dplyr
#' 
#' @param a intervals in tbl_df format
#' @param b signal in tbl_df format
#' @param col_name bare column name
#' @param operation operation to perform on intersected intervals. One of mean, median,
#'        sum, min, max, absmin, absmax 
#'        
#' @return \code{dplyr::tbl_df}
#' 
#' @export
bedtools_map <- function(a, b) {
  
  intersection <- intersect_(a, b)
  
  map_result <- intersection %>%
    group_by(chrom, .start_a, .start_b) %>%
    summarize(result = operation(col_name)) %>%
    select(chrom, col_spec, result)
  
  map_result 
}

#' @rdname bedtools_map
#' @export
bt_map <- bedtools_map