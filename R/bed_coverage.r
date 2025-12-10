#' Compute coverage of intervals.
#'
#' @param x [ivl_df]
#' @param y [ivl_df]
#' @param ... extra arguments (not used)
#'
#' @template min_overlap
#'
#' @template groups
#'
#' @family multiple set operations
#'
#' @return
#' [ivl_df] with the following additional columns:
#'
#'   - `.ints` number of `x` intersections
#'   - `.cov` per-base coverage of `x` intervals
#'   - `.len` total length of `y` intervals covered by `x` intervals
#'   - `.frac` `.len` scaled by the number of `y` intervals
#
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~strand,
#'   "chr1", 100,    500,  "+",
#'   "chr2", 200,    400,  "+",
#'   "chr2", 300,    500,  "-",
#'   "chr2", 800,    900,  "-"
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~value, ~strand,
#'   "chr1", 150,    400,  100,    "+",
#'   "chr1", 500,    550,  100,    "+",
#'   "chr2", 230,    430,  200,    "-",
#'   "chr2", 350,    430,  300,    "-"
#' )
#'
#' bed_coverage(x, y)
#'
#' @seealso \url{https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html}
#'
#' @export
bed_coverage <- function(x, y, ..., min_overlap = NULL) {
  check_required(x)
  check_required(y)

  # Handle min_overlap deprecation: default to 0 (legacy) but warn
  if (is.null(min_overlap)) {
    lifecycle::deprecate_warn(
      when = "0.8.0",
      what = "bed_coverage(min_overlap)",
      details = c(
        "The default will change from 0 (book-ended intervals overlap) to 1 (strict overlap) in a future version.",
        "Set `min_overlap = 0L` to keep the legacy behavior, or `min_overlap = 1L` for bedtools-compatible behavior."
      )
    )
    min_overlap <- 0L
  }

  x <- check_interval(x)
  y <- check_interval(y)

  x <- bed_sort(x)
  y <- bed_sort(y)

  ## add integer .id to track each input x ivl
  .id_col <- get_id_col(x)
  x[[.id_col]] <- seq_len(nrow(x))

  # establish grouping with shared groups (and chrom)
  groups_xy <- shared_groups(x, y)
  groups_xy <- unique(as.character(c("chrom", groups_xy)))
  groups_vars <- rlang::syms(groups_xy)

  # type convert grouping factors to characters if necessary and ungroup
  x <- convert_factors(x, groups_xy)
  y <- convert_factors(y, groups_xy)

  x <- group_by(x, !!!groups_vars)
  y <- group_by(y, !!!groups_vars)

  grp_indexes <- shared_group_indexes(x, y)
  res <- coverage_impl(
    x,
    y,
    grp_indexes$x,
    grp_indexes$y,
    min_overlap
  )
  res <- tibble::as_tibble(res)

  # get x ivls from groups not found in y
  mi <- get_missing_ivls(res, x, .id_col)
  res <- bind_rows(res, mi)

  # reorder by index
  res <- res[order(res[[.id_col]]), ]
  res[[.id_col]] <- NULL
  res
}

get_missing_ivls <- function(res, x, .id_col) {
  x <- ungroup(x)
  mi <- x[!x[[.id_col]] %in% res[[.id_col]], ]
  mi[[".ints"]] <- 0L
  mi[[".cov"]] <- 0L
  mi[[".len"]] <- mi[["end"]] - mi[["start"]]
  mi[[".frac"]] <- 0
  mi
}
