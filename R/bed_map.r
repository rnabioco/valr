#' Calculate summaries and statistics from overlapping intervals.
#'
#' Used to apply functions like [min()], [count()] or [concat()] to intersecting
#' intervals. Book-ended intervals are not reported by default, but can be
#' included by setting `min_overlap = 0`.
#'
#' @param x [tbl_interval()]
#' @param y  [tbl_interval()]
#' @param invert report `x` intervals not in `y`
#' @param suffix colname suffixes in output
#' @param min_overlap minimum overlap for intervals.
#' @param ... name-value pairs specifying colnames and expressions to apply
#'
#' @template groups
#'
#' @return [tbl_interval()]
#'
#' @family multiple set operations
#'
#' @seealso
#' \url{http://bedtools.readthedocs.io/en/latest/content/tools/map.html}
#'
#' @examples
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1',      1,      100
#' )
#'
#' y <- trbl_interval(
#'   ~chrom, ~start, ~end, ~value,
#'   'chr1', 1,      20,   10,
#'   'chr1', 30,     50,   20,
#'   'chr1', 90,     120,  30
#' )
#'
#' bed_glyph(bed_map(x, y, value = sum(value)), label = 'value')
#'
#' x <- trbl_interval(
#'  ~chrom, ~start, ~end,
#'  'chr1', 100,    250,
#'  'chr2', 250,    500
#' )
#'
#' y <- trbl_interval(
#'  ~chrom, ~start, ~end, ~value,
#'  'chr1', 100,    250,  10,
#'  'chr1', 150,    250,  20,
#'  'chr2', 250,    500,  500
#' )
#'
#' # also mean, median, sd etc
#' bed_map(x, y, .sum = sum(value))
#'
#' bed_map(x, y, .min = min(value), .max = max(value))
#'
#' bed_map(x, y, .concat = concat(value))
#'
#' # create a list-col
#' bed_map(x, y, .values = list(value))
#'
#' # can also use `nth` family from dplyr
#' bed_map(x, y, .first = dplyr::first(value))
#'
#' bed_map(x, y, .last = dplyr::last(value))
#'
#' bed_map(x, y, .absmax = abs(max(value)))
#'
#' bed_map(x, y, .absmin = abs(min(value)))
#'
#' bed_map(x, y, .count = length(value))
#'
#' bed_map(x, y, .count_distinct = length(unique(value)))
#'
#' bed_map(x, y, .vals = values(value))
#'
#' bed_map(x, y, .vals.unique = values_unique(value))
#'
#' @export
bed_map <- function(x, y, ..., invert = FALSE,
                    suffix = c(".x", ".y"),
                    min_overlap = 1) {

  if (!is.tbl_interval(x)) x <- as.tbl_interval(x)
  if (!is.tbl_interval(y)) y <- as.tbl_interval(y)

  groups_x <- groups(x)

  # used only to get the `x` suffix; `y` suffix is ignored`
  suffix <- list(x = suffix[1], y = suffix[2])

  # `x` names are suffixed to use for grouping later
  x_names <- colnames(x)[!colnames(x) %in% "chrom"]
  x_names_suffix <- stringr::str_c(x_names, suffix$x)

  # note that `y` columns have no suffix so can be referred to by the original
  # names
  res <- bed_intersect(x, y, invert = invert, suffix = c(suffix$x, ""))

  res <- filter(res, .overlap >= min_overlap)

  ##  map supplied functions to each set of intervals
  res <- group_by(res, !!! rlang::syms(c("chrom", x_names_suffix)))
  res <- summarize(res, !!! rlang::quos(...))
  res <- ungroup(res)
  ## remove x suffix, but don't pattern match with '.' regex
  names_no_x <- stringr::str_replace(names(res), stringr::fixed(suffix$x), "")
  names(res) <- names_no_x

  # find rows of `x` that did not intersect
  x_not <- anti_join(x, res, by = c("chrom", x_names))

  res <- bind_rows(res, x_not)
  res <- bed_sort(res)

  # reassign original `x` groups. `y` groups are gone at this point
  res <- group_by(res, !!! rlang::syms(c(groups_x)))

  res
}

#' @export
#'
#' @param .data data
#' @param sep separator character
#'
#' @rdname bed_map
concat <- function(.data, sep = ",") {
  paste0(.data, collapse = sep)
}

#' @export
#' @rdname bed_map
values_unique <- function(.data, sep = ",") {
  concat(unique(.data), sep = sep)
}

#' @export
#' @rdname bed_map
values <- function(.data, sep = ",") {
  concat(rle(.data)$values, sep = sep)
}
