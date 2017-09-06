#' Subtract two sets of intervals.
#'
#' Subtract `y` intervals from `x` intervals.
#'
#' @param x [tbl_interval()]
#' @param y [tbl_interval()]
#' @param any remove any `x` intervals that overlap `y`
#'
#' @template groups
#'
#' @family multiple set operations
#'
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/subtract.html}
#'
#' @examples
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1', 1,      100
#' )
#'
#' y <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1', 50,     75
#' )
#'
#' bed_glyph(bed_subtract(x, y))
#'
#' x <- trbl_interval(
#'  ~chrom, ~start, ~end,
#'  'chr1', 100,    200,
#'  'chr1', 250,    400,
#'  'chr1', 500,    600,
#'  'chr1', 1000,   1200,
#'  'chr1', 1300,   1500
#' )
#'
#' y <- trbl_interval(
#'  ~chrom, ~start, ~end,
#'  'chr1', 150,    175,
#'  'chr1', 510,    525,
#'  'chr1', 550,    575,
#'  'chr1', 900,    1050,
#'  'chr1', 1150,   1250,
#'  'chr1', 1299,   1501
#' )
#'
#' bed_subtract(x, y)
#'
#' bed_subtract(x, y, any = TRUE)
#'
#' @export
bed_subtract <- function(x, y, any = FALSE) {
  if (!is.tbl_interval(x)) x <- as.tbl_interval(x)
  if (!is.tbl_interval(y)) y <- as.tbl_interval(y)

  x <- group_by(x, chrom, add = TRUE)
  y <- group_by(y, chrom, add = TRUE)

  # find groups not in y
  not_y_grps <- setdiff(get_labels(x), get_labels(y))
  # keep x ivls from groups not found in y
  res_no_y <- semi_join(x, not_y_grps, by = colnames(not_y_grps))

  if (any) {
    # collect and return x intervals without overlaps
    res <- intersect_impl(x, y, invert = TRUE)
    anti <- filter(res, is.na(.overlap))
    anti <- select(anti, chrom, start = start.x, end = end.x)

    return(anti)
  }

  res <- subtract_impl(x, y)
  res <- bind_rows(res, res_no_y)
  res <- bed_sort(res)

  res
}
