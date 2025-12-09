#' Merge overlapping intervals.
#'
#' Operations can be performed on merged intervals by specifying name-value
#' pairs. Default `max_dist` of `0` means book-ended intervals are
#' merged.
#'
#' @param x [ivl_df]
#' @param max_dist maximum distance between intervals to merge
#' @param ... name-value pairs that specify operations on merged intervals
#'
#' @template groups
#'
#' @return [ivl_df]
#'
#' @family single set operations
#'
#' @seealso
#' \url{https://bedtools.readthedocs.io/en/latest/content/tools/merge.html}
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 1, 50,
#'   "chr1", 10, 75,
#'   "chr1", 100, 120
#' )
#'
#' bed_glyph(bed_merge(x))
#'
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~value, ~strand,
#'   "chr1", 1,      50,   1,      "+",
#'   "chr1", 100,    200,  2,      "+",
#'   "chr1", 150,    250,  3,      "-",
#'   "chr2", 1,      25,   4,      "+",
#'   "chr2", 200,    400,  5,      "-",
#'   "chr2", 400,    500,  6,      "+",
#'   "chr2", 450,    550,  7,      "+"
#' )
#'
#' bed_merge(x)
#'
#' bed_merge(x, max_dist = 100)
#'
#' # merge intervals on same strand
#' bed_merge(dplyr::group_by(x, strand))
#'
#' bed_merge(x, .value = sum(value))
#'
#' @export
bed_merge <- function(x, max_dist = 0, ...) {
  check_required(x)

  x <- check_interval(x)

  if (max_dist < 0) {
    cli::cli_abort("{.var max_dist} must be 0 or greater")
  }

  groups_x <- group_vars(x)

  res <- bed_sort(x)

  group_vars <- rlang::syms(unique(c("chrom", groups_x)))
  res <- group_by(res, !!!group_vars)

  # if no dots are passed then use fast internal merge
  if (!is.null(substitute(...))) {
    res <- merge_impl(res, max_dist, collapse = FALSE)
    res <- tibble::as_tibble(res)
    group_vars <- rlang::syms(unique(c("chrom", ".id_merge", groups_x)))
    res <- group_by(res, !!!group_vars)

    res <- summarize(
      res,
      !!!rlang::quos(
        .start = min(.data[["start"]]),
        .end = max(.data[["end"]]),
        ...
      )
    )
    res <- select(
      res,
      everything(),
      start = all_of(".start"),
      end = all_of(".end")
    )

    res <- ungroup(res)
    res <- select(res, !!quo(-one_of(".id_merge")))
  } else {
    res <- merge_impl(res, max_dist, collapse = TRUE)
    res <- tibble::as_tibble(res)
    res <- select(res, !!!rlang::syms(c("chrom", "start", "end", groups_x)))
  }
  res <- ungroup(res)
  # restore original grouping
  if (!is.null(groups_x)) {
    res <- group_by(res, !!!rlang::syms(groups_x))
  }

  res <- reorder_names(res, x)

  res
}
