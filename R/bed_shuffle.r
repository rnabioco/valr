#' shuffle input intervals
#' 
#' @param x tbl of intervals
#' @param genome chrom sizes
#' @param incl tbl of included intervals
#' @param excl tbl of excluded intervals
#' @param max_tries maximum tries to identify a bounded region
#' @param within shuffle within chromosomes
#' @param seed seed for reproducible intervals
#' 
#' @return \code{data_frame}
#' 
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/shuffle.html}
#' 
#' @examples 
#' x <- tibble::tribble(
#' )
#' 
#' genome <- tibble::tribble(
#' )
#' 
#' @export
bed_shuffle <- function(x, genome, incl = NULL, excl = NULL, max_tries = 100, seed = 0) {

  # flatten incl and excl
  if (!is.null(incl))
      incl <- bed_merge(incl)
  if (!is.null(excl))
      excl <- bed_merge(excl)
      
  # find the intervals to be sampled from. case where only incl intervals are
  # defined is not eval explicitly
  if (is.null(incl) && is.null(excl)) {
    incl <- genome
  } else if (is.null(incl) && !is.null(incl)) {
    incl <- bed_subtract(genome, excl)
  } else if (!is.null(incl) && !is.null(excl)) {
    incl <- bed_subtract(incl, excl)  
  }
 
  res <- shuffle_impl(x, incl, max_tries, seed)
 
  res 
}  
