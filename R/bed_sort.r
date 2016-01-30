#' Sort intervals
#'
#' @param intervals tbl of intervals
#' @param size sort by interval size
#' @param by_chrom sort within chromosome
#' @param reverse reverse sort order
#' 
#' @examples 
#' bed_df <- dplyr::tibble(
#'    ~chrom, ~start, ~end,
#'    "chr8", 500, 1000,
#'    "chr8", 1000, 5000,
#'    "chr8", 100, 200,
#'    "chr1", 100, 300,
#'    "chr1", 100, 200
#' )
#' 
#' # sort by chrom and coordinates
#' bed_df %>% bed_sort()
#' 
#' # reverse sort order
#' bed_df %>% bed_sort(reverse = TRUE)
#' 
#' # sort by interval size
#' bed_df %>% bed_sort(size = TRUE)
#' 
#' # sort by interval size within chrom
#' bed_df %>% bed_sort(size = TRUE, by_chrom = TRUE)
#' 
#' @export
bed_sort <- function(intervals, size = FALSE,
                     by_chrom = FALSE, reverse = FALSE) {

  if (size) {
    
    res <- intervals %>% mutate(.size = end - start) 
    
    if (by_chrom) {
       res <- res %>% group_by(chrom) 
    }
    
    if (reverse) {
      res <- res %>% arrange(desc(.size))
    } else {       
      res <- res %>% arrange(.size)
    }
    
    # remove .size column and groups in result
    res <- res %>% select(-.size) %>% ungroup()
    
  } else {
  
    # sort by coordinate 
    if (reverse) {
      res <- intervals %>%
        arrange(desc(chrom), desc(start), desc(end))
    } else {
      res <- intervals %>%
        arrange(chrom, start, end)
    } 
  } 
  
  res
}
