#' Compute coverage of intervals.
#' 
#' @param x tbl of intervals 
#' @param y tbl of intervals 
#' @param ... extra arguments (not used)
#'
#' @note Book-ended intervals are counted as overlapping.
#' 
#' @template groups
#' 
#' @family multi-set-ops
#' 
#' @return original \code{x} \code{data_frame} with the following additional
#'   columns:
#' \itemize{ 
#'   \item{\code{.ints}}{ number of x intersections}
#'   \item{\code{.cov}}{ per-base coverage of x intervals}
#'   \item{\code{.len}}{ total length of y intervals covered by x intervals} 
#'   \item{\code{.frac}}{ .len scaled by total of y intervals}
#'   }
#
#' @examples 
#' x <- tibble::tribble(
#' ~chrom, ~start, ~end, ~strand,
#' "chr1", 100,    500,  '+',
#' "chr2", 200,    400,  '+',
#' "chr2", 300,    500,  '-',
#' "chr2", 800,    900,  '-'
#' )
#' 
#' y <- tibble::tribble(
#' ~chrom, ~start, ~end, ~value, ~strand,
#' "chr1", 150,    400,  100,    '+',
#' "chr1", 500,    550,  100,    '+',
#' "chr2", 230,    430,  200,    '-',
#' "chr2", 350,    430,  300,    '-'
#' )
#'
#' bed_coverage(x, y)
#'  
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/coverage.html}
#'  
#' @export
bed_coverage <- function(x, y, ...) {
  
  if ( ! is_sorted(x) )
    x <- bed_sort(x)
  if ( ! is_sorted(y) )
    y <- bed_sort(y)
  
  x <- group_by(x, chrom, add = TRUE)
  y <- group_by(y, chrom, add = TRUE)
  
  res <- coverage_impl(x, y)
  
  res
}
