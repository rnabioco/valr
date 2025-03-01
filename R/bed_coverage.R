#' Compute coverage of intervals.
#'
#' @param x [ivl_df]
#' @param y [ivl_df]
#' @param ... extra arguments (not used)
#'
#' @note Book-ended intervals are included in coverage calculations.
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
bed_coverage <- function(x, y, ...) {
  check_required(x)
  check_required(y)

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
    grp_indexes$y
  )

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
