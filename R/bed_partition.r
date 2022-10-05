#' Partition intervals into elemental intervals
#'
#' Convert a set of intervals into elemental intervals that contain each start
#' and end position in the set.
#'
#' Summary operations, such as [min()] or [count()] can be performed
#' on elemental intervals by specifying name-value pairs.
#'
#' This function is useful for calculating summaries across overlapping intervals
#' without merging the intervals.
#'
#' @param x [ivl_df]
#' @param ... name-value pairs specifying column names and expressions to apply
#'
#' @template groups
#'
#' @return [ivl_df()]
#'
#' @family single set operations
#'
#' @seealso
#' \url{https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html#partition-p-partition}
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~value, ~strand,
#'  'chr1', 100,    500,  10, "+",
#'  'chr1', 200,    400,  20, "-",
#'  'chr1', 300,    550,  30, "+",
#'  'chr1', 550,    575,   2, "+",
#'  'chr1', 800,    900,   5, "+" )
#'
#'
#' bed_glyph(bed_partition(x))
#' bed_glyph(bed_partition(x, value = sum(value)), label = "value")
#'
#' bed_partition(x)
#'
#' # compute summary over each elemental interval
#' bed_partition(x, value = sum(value))
#'
#' # partition and compute summaries based on group
#' x <- dplyr::group_by(x, strand)
#' bed_partition(x, value = sum(value))
#'
#' # combine values across multiple tibbles
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~value, ~strand,
#'  'chr1', 10,     500,  100, "+",
#'  'chr1', 250,    420,  200, "-",
#'  'chr1', 350,    550,  300, "+",
#'  'chr1', 550,    555,   20, "+",
#'  'chr1', 800,    900,   50, "+" )
#'
#' x <- dplyr::bind_rows(x, y)
#' bed_partition(x, value = sum(value))
#'
#' @export
bed_partition <- function(x, ...) {

  check_required(x)
  x <- check_interval(x)

  groups_df <- group_vars(x)
  x <- bed_sort(x)

  groups <- rlang::syms(unique(c("chrom", groups_df)))
  x <- group_by(x, !!! groups)

  res <- partition_impl(x)

  res <- tibble::as_tibble(res)

  # drop non-grouped cols as values no longer match ivls
  res <- select(res, chrom, start, end, one_of(groups_df))

  # if dots are passed then map values
  if (!is.null(substitute(...))) {
    res <- group_by(res, !!! syms(groups_df))
    res <- bed_map(res, x, ...)
  }
  res
}
