#' Compute the relative distance between query intervals and reference intervals
#' 
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param detail report relative distances for each interval. Default = FALSE
#' 
#' @return \code{data_frame}
#' 
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/reldist.html}
#' 
#' @examples
#' x <- tibble::frame_data(
#' ~chrom,   ~start,    ~end,
#' "chr1",    75,       125
#'   )
#' 
#' y <- tibble::frame_data(
#'   ~chrom,   ~start,    ~end,
#'   "chr1",    50,       100,
#'   "chr1",    100,       150
#'   )
#' 
#' bed_reldist(x, y)
#' 
#' 
#' @export

bed_reldist <- function(x, y, detail = FALSE) {
  
  if ( ! is_sorted(x) )
    x <- bed_sort(x)
  if ( ! is_sorted(y) )
    y <- bed_sort(y)
 
  x <- dplyr::group_by(x, chrom, add = TRUE)
  y <- dplyr::group_by(y, chrom, add = TRUE)
  
  res <- reldist_impl(x, y)
  
  if (!detail){
    res$reldist <- floor(res$reldist * 100) / 100
    total_ivls <- nrow(res)
    res <- dplyr::group_by(res, reldist)
    res <- dplyr::summarize(res, 
                     counts = n(),
                     total = total_ivls,
                     freq = counts / total)
  } 
  
  res
  
}
  
  
