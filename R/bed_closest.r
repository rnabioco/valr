#' Identify closest intervals.
#'
#' @param x [tbl_interval()]
#' @param y [tbl_interval()]
#' @param overlap report overlapping intervals
#' @param suffix colname suffixes in output
#'
#' @template groups
#'
#' @return
#' [tbl_interval()] with additional columns:
#'   - `.dist` distance to closest interval. Negative distances
#'     denote upstream intervals.
#'   - `.overlap` overlap with closest interval
#'
#' @family multiple set operations
#'
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/closest.html}
#'
#' @examples
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1', 100,    125
#' )
#'
#' y <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1', 25,     50,
#'   'chr1', 140,    175
#' )
#'
#' bed_glyph(bed_closest(x, y))
#'
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   "chr1", 500,    600,
#'   "chr2", 5000,   6000
#' )
#'
#' y <- trbl_interval(
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
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 10,	   20,   "a",   1,      "-"
#' )
#'
#' y <- trbl_interval(
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
                        suffix = c(".x", ".y")){

  if (!is.tbl_interval(x)) x <- as.tbl_interval(x)
  if (!is.tbl_interval(y)) y <- as.tbl_interval(y)

  check_suffix(suffix)

  x <- bed_sort(x)
  x <- group_by(x, chrom, add = TRUE)

  y <- bed_sort(y)
  y <- group_by(y, chrom, add = TRUE)

  suffix <- list(x = suffix[1], y = suffix[2])

  res <- closest_impl(x, y, suffix$x, suffix$y)

  if (!overlap){
    res <- filter(res, .overlap < 1)
    res <- select(res, -.overlap)
  }

  res
}
