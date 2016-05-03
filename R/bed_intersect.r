#' Identify intersecting intervals.
#' 
#' @param x tbl of intervals 
#' @param y tbl of intervals 
#' @param strand intersect intervals on same strand
#' @param strand_opp intersect intervals on opposite strands
#' @param suffix_x suffix for intersected intervals from x (except chrom)
#' @param suffix_y suffix for intersected intervals from y (except chrom)
#'
#' @note Book-ended intervals have \code{.overlap} values of 0 in the output.
#'  
#' @examples 
#' x <- tibble::frame_data(
#' ~chrom, ~start, ~end,
#' "chr1", 100,    500,
#' "chr2", 200,    400,
#' "chr2", 300,    500,
#' "chr2", 800,    900
#' )
#' 
#' y <- tibble::frame_data(
#' ~chrom, ~start, ~end, ~value,
#' "chr1", 150,    400,  100,
#' "chr1", 500,    550,  100,
#' "chr2", 230,    430,  200,
#' "chr2", 350,    430,  300
#' )
#'
#' bed_intersect(x, y)
#'  
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html}
#'  
#' @export
bed_intersect <- function(x, y, strand = FALSE, strand_opp = FALSE,
                          suffix_x = '.x', suffix_y = '.y') {
 
  if ( ! is_sorted(x) )
    x <- bed_sort(x)
  if ( ! is_sorted(y) )
    y <- bed_sort(y)
 
  if (is.null(groups(x)) || groups(x) != "chrom")
    x <- group_by(x, chrom)
  if (is.null(groups(y)) || groups(y) != "chrom")
    y <- group_by(y, chrom)

  res <- intersect_impl(x, y, suffix_x, suffix_y)
  
  if (strand) {
    if (! 'strand' %in% colnames(res)){
      stop("strand arg specified on unstranded data_frame", .Call = FALSE)
    }
     res <- filter(res, strand.x == strand.y) 
  } else if (strand_opp) {
     res <- filter(res, strand.x != strand.y) 
  }
  
  res
}


