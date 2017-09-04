#' Compute coverage of intervals.
#'
#' @param x [tbl_interval()]
#' @param y [tbl_interval()]
#' @param ... extra arguments (not used)
#'
#' @note Book-ended intervals are included in coverage calculations.
#'
#' @template groups
#'
#' @family multiple set operations
#'
#' @return
#' [tbl_interval()] with the following additional columns:
#'
#'   - `.ints` number of `x` intersections
#'   - `.cov` per-base coverage of `x` intervals
#'   - `.len` total length of `y` intervals covered by `x` intervals
#'   - `.frac` `.len` scaled by the number of `y` intervals
#
#' @examples
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end, ~strand,
#'   "chr1", 100,    500,  '+',
#'   "chr2", 200,    400,  '+',
#'   "chr2", 300,    500,  '-',
#'   "chr2", 800,    900,  '-'
#' )
#'
#' y <- trbl_interval(
#'   ~chrom, ~start, ~end, ~value, ~strand,
#'   "chr1", 150,    400,  100,    '+',
#'   "chr1", 500,    550,  100,    '+',
#'   "chr2", 230,    430,  200,    '-',
#'   "chr2", 350,    430,  300,    '-'
#' )
#'
#' bed_coverage(x, y)
#'
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/coverage.html}
#'
#' @export
bed_coverage <- function(x, y, ...) {
  if (!is.tbl_interval(x)) x <- as.tbl_interval(x)
  if (!is.tbl_interval(y)) y <- as.tbl_interval(y)

  x <- bed_sort(x)
  x <- group_by(x, chrom, add = TRUE)

  y <- bed_sort(y)
  y <- group_by(y, chrom, add = TRUE)

  res <- coverage_impl(x, y)

  res
}
