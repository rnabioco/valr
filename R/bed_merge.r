#' Merge overlapping intervals.
#'
#' Operations can be performed on merged intervals by specifying name-value
#' pairs. Default `max_dist` of `0` means book-ended intervals are
#' merged.
#'
#' @param x [tbl_interval()]
#' @param max_dist maximum distance between intervals to merge
#' @param ... name-value pairs that specify operations on merged intervals
#'
#' @template groups
#'
#' @return [tbl_interval()]
#'
#' @family single set operations
#'
#' @seealso
#' \url{http://bedtools.readthedocs.org/en/latest/content/tools/merge.html}
#'
#' @examples
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1',  1,      50,
#'   'chr1',  10,     75,
#'   'chr1',  100,    120
#' )
#'
#' bed_glyph(bed_merge(x))
#'
#' x <- trbl_interval(
#'  ~chrom, ~start, ~end, ~value, ~strand,
#'  "chr1", 1,      50,   1,      '+',
#'  "chr1", 100,    200,  2,      '+',
#'  "chr1", 150,    250,  3,      '-',
#'  "chr2", 1,      25,   4,      '+',
#'  "chr2", 200,    400,  5,      '-',
#'  "chr2", 400,    500,  6,      '+',
#'  "chr2", 450,    550,  7,      '+'
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
  if (!is.tbl_interval(x)) x <- as.tbl_interval(x)

  if (max_dist < 0) {
    stop("max_dist must be positive", call. = FALSE)
  }

  if (is_merged(x)) return(x)

  groups_x <- groups(x)

  res <- bed_sort(x)

  res <- group_by(res, chrom, add = TRUE)

  # if no dots are passed then use fast internal merge
  if (!is.null(substitute(...))) {
    res <- merge_impl(res, max_dist, collapse = FALSE)
    group_vars <- rlang::syms(c("chrom", ".id_merge", groups_x))
    res <- group_by(res, !!! group_vars, add = TRUE)

    res <- summarize(res, !!! rlang::quos(
      .start = min(start),
      .end = max(end),
      ...
    ))
    res <- select(res, everything(), start = .start, end = .end)

    res <- ungroup(res)
    res <- select(res, !! quo(-one_of(".id_merge")))
  } else {
    res <- merge_impl(res, max_dist, collapse = TRUE)
    res <- select(res, !!! rlang::syms(c("chrom", "start", "end", groups_x)))
  }

  # restore original grouping
  res <- group_by(res, !!! rlang::syms(groups_x))
  res <- reorder_names(res, x)

  attr(res, "merged") <- TRUE

  res
}

#' Ask whether a tbl is already merged.
#'
#' @param x tbl of intervals
#' @noRd
is_merged <- function(x) {
  merged_attr <- attr(x, "merged")

  if (is.null(merged_attr) || !merged_attr) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}
