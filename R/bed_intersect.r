#' Identify intersecting intervals.
#' 
#' @param x tbl of intervals 
#' @param y tbl of intervals 
#' @param strand intersect intervals on same strand
#' @param strand_opp intersect intervals on opposite strands
#' @param suffix colname suffixes in output
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
bed_intersect <- function(x, y, strand = FALSE, strand_opp = FALSE, suffix = c('.x', '.y')) {
 
  if ( ! is_sorted(x) )
    x <- bed_sort(x)
  if ( ! is_sorted(y) )
    y <- bed_sort(y)
 
  if (is.null(groups(x)) || groups(x) != "chrom")
    x <- group_by(x, chrom)
  if (is.null(groups(y)) || groups(y) != "chrom")
    y <- group_by(y, chrom)

  # dplyr::check_suffix
  if (!is.character(suffix) || length(suffix) != 2) {
    stop("`suffix` must be a character vector of length 2.", call. = FALSE)
  }
  
  suffix <- list(x = suffix[1], y = suffix[2])

  res <- intersect_impl(x, y, suffix$x, suffix$y)
  
  if (strand) {
    if (!('strand' %in% colnames(x) && 'strand' %in% colnames(y))){
      stop("`strand` specified on unstranded data_frame", call. = FALSE)
    }
     res <- filter(res, strand.x == strand.y) 
  } else if (strand_opp) {
     res <- filter(res, strand.x != strand.y) 
  }
  
  res
}


