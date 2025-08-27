#' Calculate interval spacing.
#'
#' Spacing for the first interval of each chromosome is undefined (`NA`). The
#' leading interval of an overlapping interval pair has a negative value.
#'
#' @param x [ivl_df]
#'
#' @return [ivl_df] with `.spacing` column.
#'
#' @family utilities
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 1,      100,
#'   "chr1", 150,    200,
#'   "chr2", 200,    300
#' )
#'
#' interval_spacing(x)
#'
#' @export
interval_spacing <- function(x) {
  x <- check_interval(x)

  res <- bed_sort(x)

  gx <- groups(x)

  res <- group_by(res, .data[["chrom"]])
  res <- mutate(res, .spacing = .data[["start"]] - lag(.data[["end"]]))

  res <- group_by(res, !!!gx)

  res
}
