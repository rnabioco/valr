#' Identify intervals in a genome not covered by a query.
#' 
#' @param x tbl of intervals
#' @param genome chrom sizes
#' 
#' @return \code{data_frame}
#' 
#' @examples 
#' genome <- tibble::frame_data(
#'    ~chrom,  ~size,
#'    "chr1", 500,
#'    "chr2", 600,
#'    "chr3", 800
#' ) 
#' 
#' x <- tibble::frame_data(
#'    ~chrom, ~start, ~end,
#'    "chr1", 100,    300,
#'    "chr1", 200,    400,
#'    "chr2", 1,      100,
#'    "chr2", 200,    400,
#'    "chr3", 500,    600
#' )
#' 
#' # intervals not covered by x
#' bed_complement(x, genome)
#' 
#' @export
bed_complement <- function(x, genome) {

  if ( ! is_merged(x) ) {
    res <- bed_merge(x)
  } 
  
  res <- dplyr::group_by(res, chrom)
  
  res <- complement_impl(res, genome)

  res <- bed_sort(res)
 
  res 
}
