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
#' @param x [ivl_df]
#' @param y  [ivl_df]
#' @param ... name-value pairs specifying column names and expressions to apply
#' @param min_overlap minimum overlap for intervals.
#'
#' @template groups
#'
#' @return [ivl_df]
#'
#' @family multiple set operations
#'
#' @seealso
#' \url{https://bedtools.readthedocs.io/en/latest/content/tools/map.html}
#'
#' @example inst/example/bed_map.r
#'
#' @export
bed_map <- function(x, y, ..., min_overlap = 1) {
  x <- check_interval(x)
  y <- check_interval(y)

  ## add suffixes to all x columns except `chrom`
  x_nms <- str_c(names(x)[!names(x) %in% "chrom"], ".x")

  ## add integer .id to track each input x ivl
  ## don't use mutate, in order to ignore input grouping
  ## don't overwrite .id if it exists
  .id_col <- ".id"
  if (.id_col %in% names(x)) {
    .id_col <- str_c(.id_col, ".x")
  }
  x[[.id_col]] <- seq_len(nrow(x))

  .id_col_out <- str_c(.id_col, ".x")

  # establish grouping with shared groups (and chrom)
  groups_xy <- shared_groups(x, y)
  groups_xy <- unique(as.character(c("chrom", groups_xy)))
  groups_vars <- rlang::syms(groups_xy)

  # type convert grouping factors to characters if necessary and ungroup
  x <- convert_factors(x, groups_xy)
  y <- convert_factors(y, groups_xy)

  x <- group_by(x, !!! groups_vars)
  y <- group_by(y, !!! groups_vars)

  grp_indexes <- shared_group_indexes(x, y)

  res <- intersect_impl(x, y,
                        grp_indexes$x,
                        grp_indexes$y,
                        invert = TRUE,
                        suffix_x = ".x", suffix_y = "")

  ## filter for rows that don't intersect. The `duplicated` call is required
  ## because book-ended intervals in the intersect_impl result can
  ## book-end multiple `y` intervals, causing them to be duplicated after the
  ## `select`. base::duplicated is ~10x faster than dplyr::distinct
  res_noint <- filter(res, is.na(.overlap) | .overlap < min_overlap)
  res_noint <- select(res_noint, chrom, ends_with(".x"))
  res_noint <- res_noint[!duplicated(res_noint[[.id_col_out]]), ]

  ## identify intersecting intervals
  res_int <- filter(res, !is.na(.overlap) & .overlap >= min_overlap)

  ## drop non-intersecting intervals that are also in the intersecting set
  ## this prevents duplicate reporting of an x interval if it both bookends
  ## and overlaps y intervals. Using base R logical indexing here is ~ 7x faster
  ## than dplyr::anti_join()
  res_noint <- res_noint[!res_noint[[.id_col_out]] %in% res_int[[.id_col_out]], ]
  res_noint <- select(res_noint, -contains(.id_col_out))

  ## map supplied functions to each set of intersecting intervals
  ## group_by .id_col_out to ensure that duplicated input x intervals are reported
 # res_int <- group_by(res_int, !!! syms(c("chrom", x_nms, .id_col_out)))
  res_int <- group_by(res_int, !!! syms(.id_col_out))
  res_int <- summarize(res_int, !!! quos(...))
  res_int <- ungroup(res_int)

  ## join summarize data with intervals based on .id
  ## avoids grouping by many unnecessary columns
  join_grps <- .id_col_out
  names(join_grps) <- .id_col
  x <- ungroup(x)
  res_int <- inner_join(x, res_int, by = join_grps)
  res_int <- select(res_int, -contains(.id_col))

  ## rename to match input columns
  colnames(res_noint) <- stringr::str_replace(colnames(res_noint), "[.]x$", "")

  res <- bind_rows(res_int, res_noint)

  res <- bed_sort(res)

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
