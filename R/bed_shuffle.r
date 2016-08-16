#' shuffle input intervals
#' 
#' @param x tbl of intervals
#' @param genome chrom sizes
#' @param incl tbl of included intervals
#' @param excl tbl of excluded intervals
#' @param max_tries maximum tries to identify a bounded region
#' 
#' @return \code{data_frame}
#' 
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/shuffle.html}
#' 
#' @examples 
#' x <- tibble::frame_data(
#' )
#' 
#' genome <- tibble::frame_data(
#' )
#' 
#' @export
bed_shuffle <- function(x, genome, incl = NULL, excl = NULL, max_tries = 100) {
 
  # create empty data frames if incl / excl are null 
  if (is.null(incl))
    incl <- data_frame()
  if (is.null(excl)) 
    excl <- data_frame()
  
  res <- shuffle_impl(x, genome, incl, excl, max_tries)
  res
}  
