#' Sort intervals
#'
#' @param intervals tbl of intervals
#' @param by_size sort by interval size
#' @param by_chrom sort within chromosome
#' @param reverse reverse sort order
#' 
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/sort.html}
#'
#' @examples
#' bed_tbl <- tibble::frame_data(
#'    ~chrom, ~start, ~end,
#'    "chr8", 500, 1000,
#'    "chr8", 1000, 5000,
#'    "chr8", 100, 200,
#'    "chr1", 100, 300,
#'    "chr1", 100, 200
#' )
#' 
#' # sort by chrom and start
#' bed_sort(bed_tbl)
#' 
#' # reverse sort order
#' bed_sort(bed_tbl, reverse = TRUE)
#' 
#' # sort by interval size
#' bed_sort(bed_tbl, by_size = TRUE)
#' 
#' # sort by decreasing interval size
#' bed_sort(bed_tbl, by_size = TRUE, reverse = TRUE)
#' 
#' # sort by interval size within chrom
#' bed_sort(bed_tbl, by_size = TRUE, by_chrom = TRUE)
#' 
#' @export
bed_sort <- function(intervals, by_size = FALSE,
                     by_chrom = FALSE, reverse = FALSE) {

  if (by_size) {
    
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
    res <- res %>% select(-.size)
    
  } else {
  
    if (by_chrom) {
       res <- res %>% group_by(chrom) 
    }
    
    # sort by coordinate 
    if (reverse) {
      res <- intervals %>%
        arrange(chrom, desc(start))
    } else {
      res <- intervals %>%
        arrange(chrom, start)
    } 
  } 
 
  # remove groups in result 
  res <- res %>% ungroup() %>% tbl_df
 
  # add `sorted` attribute 
  attr(res, "sorted") <- TRUE
  
  res
}

#' determine whether tbl has been previously sorted
#' 
#' @export
is_sorted <- function(bed_tbl) {
  
  sorted_attr <- attr(bed_tbl, "sorted")
  
  if (is.null(sorted_attr) || ! sorted_attr) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}
