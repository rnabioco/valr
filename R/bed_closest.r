#' Identify closest intervals.
#'
#' @param x [ivl_df]
#' @param y [ivl_df]
#' @param overlap report overlapping intervals
#' @param suffix colname suffixes in output
#'
#' @note For each interval in x `bed_closest()` returns overlapping intervals from y
#' and the closest non-intersecting y interval. Setting `overlap = FALSE` will
#' report the closest non-intersecting y intervals, ignoring any overlapping y
#' intervals.
#'
#' @template groups
#'
#' @return
#' [ivl_df] with additional columns:
#'   - `.overlap` amount of overlap with overlapping interval. Non-overlapping
#'   or adjacent intervals have an overlap of 0. `.overlap` will not be included
#'   in the output if `overlap = FALSE`.
#'   - `.dist` distance to closest interval. Negative distances
#'   denote upstream intervals. Book-ended intervals have a distance of 1.
#'
#' @family multiple set operations
#'
#' @seealso \url{https://bedtools.readthedocs.io/en/latest/content/tools/closest.html}
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 100,    125
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 25,     50,
#'   "chr1", 140,    175
#' )
#'
#' bed_glyph(bed_closest(x, y))
#'
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 500,    600,
#'   "chr2", 5000,   6000
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 100,    200,
#'   "chr1", 150,    200,
#'   "chr1", 550,    580,
#'   "chr2", 7000,   8500
#' )
#'
#' bed_closest(x, y)
#'
#' bed_closest(x, y, overlap = FALSE)
#'
#' # Report distance based on strand
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 10, 20, "a", 1, "-"
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 8, 9, "b", 1, "+",
#'   "chr1", 21, 22, "b", 1, "-"
#' )
#'
#' res <- bed_closest(x, y)
#'
#' # convert distance based on strand
#' res$.dist_strand <- ifelse(res$strand.x == "+", res$.dist, -(res$.dist))
#' res
#'
#' # report absolute distances
#' res$.abs_dist <- abs(res$.dist)
#' res
#'
#' @export
bed_closest <- function(x, y, overlap = TRUE, suffix = c(".x", ".y")) {
  x <- check_interval(x)
  y <- check_interval(y)

  check_suffix(suffix)

  x <- bed_sort(x)
  y <- bed_sort(y)

  id <- get_id_col(x)
  x[[id]] <- seq_len(nrow(x))
  x_id_out <- paste0(id, suffix[1])

  # establish grouping with shared groups (and chrom)
  groups_xy <- shared_groups(x, y)
  groups_xy <- unique(as.character(c("chrom", groups_xy)))
  groups_vars <- rlang::syms(groups_xy)

  # type convert grouping factors to characters if necessary and ungroup
  x <- convert_factors(x, groups_xy)
  y <- convert_factors(y, groups_xy)

  x <- group_by(x, !!!groups_vars)
  y <- group_by(y, !!!groups_vars)

  # TODO: revisit when min_overlap default changes to 1L

  ol_ivls <- bed_intersect(x, y, suffix = suffix, min_overlap = 0L)

  grp_indexes <- shared_group_indexes(x, y)

  res <- closest_impl(
    x,
    y,
    grp_indexes$x,
    grp_indexes$y,
    suffix[1],
    suffix[2]
  )

  # Column names with suffixes for use in mutate
  end_y <- paste0("end", suffix[2])
  start_x <- paste0("start", suffix[1])

  res$.overlap <- 0L
  ol_ivls <- mutate(
    ol_ivls,
    .dist = case_when(
      .data[[".overlap"]] > 0 ~ 0L,
      .data[[end_y]] <= .data[[start_x]] ~ -1L,
      TRUE ~ 1L
    )
  )

  res <- res[colnames(ol_ivls)]
  res <- bind_rows(ol_ivls, res)

  # get x ivls from groups not found in y
  x <- add_colname_suffix(x, suffix[1])
  mi <- get_no_group_ivls(x, res, x_id_out)
  if (nrow(mi) > 0) {
    res <- bind_rows(res, mi)
  }

  if (!overlap) {
    res <- res[which(res$.overlap <= 0 | is.na(res$.overlap)), ]
    res[[".overlap"]] <- NULL
  }

  # reorder by input x ivls
  res <- res[order(res[[x_id_out]], method = "radix"), ]
  res[[x_id_out]] <- NULL
  res
}

add_colname_suffix <- function(x, s) {
  cidx <- which(names(x) != "chrom")
  new_ids <- str_c(colnames(x), s)
  colnames(x)[cidx] <- new_ids[cidx]
  x
}

get_no_group_ivls <- function(x, y, cid) {
  x <- ungroup(x)
  mi <- x[!x[[cid]] %in% y[[cid]], ]
  mi
}
