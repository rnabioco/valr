#' Calculate summaries from overlapping intervals.
#'
#' Apply functions like [min()] and [count()] to intersecting intervals.
#' [bed_map()] uses [bed_intersect()] to identify intersecting intervals, so
#' output columns will be suffixed with `.x` and `.y`. Expressions that refer to
#' input columns from `x` and `y` columns must take these suffixes into account.
#'
#' Book-ended intervals can be included by setting `min_overlap = 0`.
#' Non-intersecting intervals from `x` are included in the result with `NA`
#' values
#'
#' @param x [tbl_interval()]
#' @param y  [tbl_interval()]
#' @param ... name-value pairs specifying column names and expressions to apply
#' @param min_overlap minimum overlap for intervals.
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
#' @example inst/example/bed_map.r
#'
#' @export
bed_map <- function(x, y, ..., min_overlap = 1) {
  if (!is.tbl_interval(x)) x <- as.tbl_interval(x)
  if (!is.tbl_interval(y)) y <- as.tbl_interval(y)

  ## add suffixes to all x columns except `chrom`
  x_nms <- str_c(names(x)[!names(x) %in% "chrom"], ".x")

  ## add chrom as a group
  x <- group_by(x, chrom, add = TRUE)
  y <- group_by(y, chrom, add = TRUE)

  res <- intersect_impl(x, y, invert = TRUE)

  ## filter for rows that don't intersect. The `distinct` call is required
  ## because book-ended book-ended intervals in the intersect_impl result can
  ## book-end multiple `y` intervals, causing them to be duplicated after the
  ## `select`
  res_noint <- filter(res, is.na(.overlap) | .overlap < min_overlap)
  res_noint <- select(res_noint, chrom, ends_with('.x'))
  res_noint <- distinct(res_noint)

  ## map supplied functions to each set of intervals
  res_int <- filter(res, !is.na(.overlap) & .overlap >= min_overlap)
  res_int <- group_by(res_int, !!! syms(c("chrom", x_nms)))
  res_int <- summarize(res_int, !!! quos(...))
  res_int <- ungroup(res_int)

  res <- bind_rows(res_int, res_noint)
  # can't use bed_sort because of the `.x` suffixes
  res <- arrange(res, chrom, start.x, end.x)

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
