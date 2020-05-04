#' Identify closest intervals.
#'
#' @param x [ivl_df]
#' @param y [ivl_df]
#' @param overlap report overlapping intervals
#' @param suffix colname suffixes in output
#'
#' @template groups
#'
#' @return
#' [ivl_df] with additional columns:
#'   - `.dist` distance to closest interval. Negative distances
#'     denote upstream intervals.
#'   - `.overlap` overlap with closest interval
#'
#' @family multiple set operations
#'
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/closest.html}
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 100,    125
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 25,     50,
#'   'chr1', 140,    175
#' )
#'
#' bed_glyph(bed_closest(x, y))
#'
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 500,    600,
#'   "chr2", 5000,   6000
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 100,    200,
#'   "chr1", 150,    200,
#'   "chr1", 550,    580,
#'   "chr2", 7000,   8500
#' )
#'
#' bed_closest(x, y)
#'
#' bed_closest(x, y, overlap = FALSE)
#'
#' # Report distance based on strand
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 10,	   20,   "a",   1,      "-"
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 8,	     9,	   "b",   1,      "+",
#'   "chr1", 21,	   22,	 "b",   1,      "-"
#' )
#'
#' res <- bed_closest(x, y)
#'
#' # convert distance based on strand
#' res$.dist_strand <- ifelse(res$strand.x == "+", res$.dist, -(res$.dist))
#' res
#'
#' # report absolute distances
#' res$.abs_dist <- abs(res$.dist)
#' res
#'
#' @export
bed_closest <- function(x, y, overlap = TRUE,
                        suffix = c(".x", ".y")) {
  x <- check_interval(x)
  y <- check_interval(y)

  check_suffix(suffix)

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

  suffix <- list(x = suffix[1], y = suffix[2])

  grp_indexes <- shared_group_indexes(x, y)

  res <- closest_impl(x, y,
                      grp_indexes$x,
                      grp_indexes$y,
                      suffix$x,
                      suffix$y)

  if (!overlap) {
    res <- filter(res, .overlap < 1)
    res <- select(res, -.overlap)
  }

  res
}
