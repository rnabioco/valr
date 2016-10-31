#' Cluster neighboring intervals.
#' 
#' Returned \code{data_frame} contains a new \code{.id} column that can be used for
#' grouping along with \code{chrom}. Default \code{max_dist} of 0 means that
#' overlapping and book-ended intervals will be clustered.
#' 
#' @param x tbl of intervals
#' @param max_dist maximum distance between clustered intervals.
#'   
#' @template groups
#' 
#' @return \code{data_frame}
#' @family single-set-ops
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/cluster.html} 
#'
#' @examples
#' x <- tibble::tribble(
#'  ~chrom, ~start, ~end,
#'  "chr1", 100,  200,
#'  "chr1", 180,  250,
#'  "chr1", 250,  500,
#'  "chr1", 501,  1000
#' )
#' 
#' bed_cluster(x)
#' 
#' # glyph illustrating clustering of overlapping and book-ended intervals
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 1,      10,
#'   'chr1', 5,      20,
#'   'chr1', 30,     40,
#'   'chr1', 40,     50,
#'   'chr1', 80,     90
#' )
#' 
#' bed_glyph(bed_cluster(x), label = '.id')
#'
#' @export
bed_cluster <- function(x, max_dist = 0) {

  res <- group_by(x, chrom, add = TRUE)
  
  res <- merge_impl(res, max_dist)
    
  res <- group_by(res, chrom)
  res <- mutate(res, .id = dense_rank(.id_merge))
  res <- select(res, -.id_merge, -.overlap_merge)
  res <- ungroup(res)
  
  res
}
