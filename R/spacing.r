#' Calculate interval spacing.
#'
#' Overlapping intervals are merged. Spacing for the first interval of each
#' chromosome is undefined (\code{NA}).
#'
#' @param x \code{\link{tbl_interval}}
#'
#' @return \code{\link{tbl_interval}} with \code{.spacing} column.
#'
#' @examples
#' x <- tibble::tribble(
#' ~chrom, ~start, ~end,
#' 'chr1', 1,      100,
#' 'chr1', 150,    200,
#' 'chr2', 200,    300
#' )
#'
#' interval_spacing(x)
#'
#' @export
interval_spacing <- function(x) {

  if (!is.tbl_interval(x)) x <- tbl_interval(x)

  res <- bed_merge(x)

  groups_x <- groups(x)

  res <- group_by(res, chrom)
  res <- mutate(res, .spacing = start - lag(end))

  res <- group_by_(res, .dots = groups_x)

  res
}
