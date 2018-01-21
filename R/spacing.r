#' Calculate interval spacing.
#'
#' Spacing for the first interval of each chromosome is undefined (`NA`). The
#' leading interval of an overlapping interval pair has a negative value.
#'
#' @param x [tbl_interval()]
#'
#' @return [tbl_interval()] with `.spacing` column.
#'
#' @family utilities
#'
#' @examples
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1', 1,      100,
#'   'chr1', 150,    200,
#'   'chr2', 200,    300
#' )
#'
#' interval_spacing(x)
#'
#' @export
interval_spacing <- function(x) {
  if (!is.tbl_interval(x)) x <- as.tbl_interval(x)

  res <- bed_sort(x)

  gx <- groups(x)

  res <- group_by(res, chrom)
  res <- mutate(res, .spacing = start - lag(end))

  res <- group_by(res, !!! gx)

  res
}
