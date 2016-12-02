#' Sorting tbls of intervals
#' 
#' Interval tbls can be sorted using \code{\link[dplyr]{arrange}}. See examples 
#' for sorting by chrom and start, and by size using
#' \code{\link[dplyr]{mutate}}.
#' 
#' @name sorting
#'   
#' @examples
#' 
#' # unsorted tbl
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 150,    500,
#'   'chr1', 1,      100,
#'   'chr2', 500,    1000,
#'   'chr2', 100,    300
#' ) 
#' 
#' # sort by start
#' dplyr::arrange(x, start)
#' 
#' # sort by descending start
#' dplyr::arrange(x, desc(start))
#' 
#' # sort by chrom and start
#' dplyr::arrange(x, chrom, start)
#' 
#' # sort by size
#' x <- dplyr::mutate(x, .size = end - start)
#' dplyr::arrange(x, .size)
#' 
NULL

#' Sort a tbl of intervals.
#' 
#' Multiple sorting parameters can be combined. note that \code{by_chrom} sorts 
#' within a chrom, not by chrom.
#' 
#' Sorting strips groups from the input.
#' 
#' @section Deprecated function: \code{bed_sort} has been deprecated. Instead
#'   use \code{\link[dplyr]{arrange}}. See \link{sorting} for details.
#'   
#' @param x tbl of intervals
#' @param by_size sort by interval size
#' @param by_chrom sort within chromosome
#' @param reverse reverse sort order
#'   
#' @seealso 
#' \url{http://bedtools.readthedocs.org/en/latest/content/tools/sort.html}
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
bed_sort <- function(x, by_size = FALSE, by_chrom = FALSE, reverse = FALSE) {

  .Deprecated(msg = 'bed_sort() is deprecated. see examples in `?sorting`.')
  
  if (by_size) {
    
    res <- mutate(x, .size = end - start) 
    
    if (by_chrom) {
       res <- group_by(res, chrom) 
    }
    
    if (reverse) {
      res <- arrange(res, desc(.size))
    } else {       
      res <- arrange(res, .size)
    }
    
    # remove .size column and groups in result
    res <- select(res, -.size)
    
  } else {
  
    if (by_chrom) {
       res <- group_by(x, chrom) 
    }
    
    # sort by coordinate 
    if (reverse) {
      res <- arrange(x, chrom, desc(start))
    } else {
      res <- arrange(x, chrom, start)
    } 
  } 
 
  # add `sorted` attribute 
  attr(res, "sorted") <- TRUE
  
  res
}

