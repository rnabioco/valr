#' Sort a tbl of intervals.
#' 
#' Multiple sorting parameters can be combined. note that \code{by_chrom} sorts
#' within a chrom, not by chrom.
#' 
#' Sorting strips groups from the input.
#' 
#' @param x tbl of intervals
#' @param by_size sort by interval size
#' @param by_chrom sort within chromosome
#' @param reverse reverse sort order
#'   
#' @seealso
#'   \url{http://bedtools.readthedocs.org/en/latest/content/tools/sort.html}
#'   
#' @examples
#' x <- tibble::tribble(
#'    ~chrom, ~start, ~end,
#'    "chr8", 500, 1000,
#'    "chr8", 1000, 5000,
#'    "chr8", 100, 200,
#'    "chr1", 100, 300,
#'    "chr1", 100, 200
#' )
#' 
#' # sort by chrom and start
#' bed_sort(x)
#' 
#' # reverse sort order
#' bed_sort(x, reverse = TRUE)
#' 
#' # sort by interval size
#' bed_sort(x, by_size = TRUE)
#' 
#' # sort by decreasing interval size
#' bed_sort(x, by_size = TRUE, reverse = TRUE)
#' 
#' # sort by interval size within chrom
#' bed_sort(x, by_size = TRUE, by_chrom = TRUE)
#' 
#' @export
bed_sort <- function(x, by_size = FALSE,
                     by_chrom = FALSE, reverse = FALSE) {

  if (by_size) {
    
    res <- dplyr::mutate(x, .size = end - start) 
    
    if (by_chrom) {
       res <- dplyr::group_by(res, chrom) 
    }
    
    if (reverse) {
      res <- dplyr::arrange(res, desc(.size))
    } else {       
      res <- dplyr::arrange(res, .size)
    }
    
    # remove .size column and groups in result
    res <- dplyr::select(res, -.size)
    
  } else {
  
    if (by_chrom) {
       res <- dplyr::group_by(x, chrom) 
    }
    
    # sort by coordinate 
    if (reverse) {
      res <- dplyr::arrange(x, chrom, desc(start))
    } else {
      res <- dplyr::arrange(x, chrom, start)
    } 
  } 
 
  # remove groups in result 
  res <- dplyr::ungroup(res) 
  res <- tibble::as_tibble(res)
 
  # add `sorted` attribute 
  attr(res, "sorted") <- TRUE
  
  res
}

#' Ask whether tbl is sorted.
#' 
#' @param x tbl of intervals
#' 
#' @export
is_sorted <- function(x) {
  
  sorted_attr <- attr(x, "sorted")
  
  if (is.null(sorted_attr) || ! sorted_attr) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}
