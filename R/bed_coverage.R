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
  y <- bed_sort(y)

  # establish grouping with shared groups (and chrom)
  groups_xy <- shared_groups(x, y)
  groups_xy <- unique(as.character(c("chrom", groups_xy)))
  groups_vars <- rlang::syms(groups_xy)

  # type convert grouping factors to characters if necessary and ungroup
  x <- convert_factors(x, groups_xy)
  y <- convert_factors(y, groups_xy)

  x <- group_by(x, !!! groups_vars)
  y <- group_by(y, !!! groups_vars)

  if (utils::packageVersion("dplyr") < "0.7.99.9000"){
    x <- update_groups(x)
    y <- update_groups(y)
  }

  grp_indexes <- shared_group_indexes(x, y)
  res <- coverage_impl(x, y,
                       grp_indexes$x,
                       grp_indexes$y)

  res
}
