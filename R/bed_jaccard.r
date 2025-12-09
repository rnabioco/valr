#' Calculate the Jaccard statistic for two sets of intervals.
#'
#' Quantifies the extent of overlap between to sets of intervals in terms of
#' base-pairs. Groups that are shared between input are used to calculate the statistic
#' for subsets of data.
#'
#' @details The Jaccard statistic takes values of `[0,1]` and is measured as:
#'
#' \deqn{ J(x,y) = \frac{\mid x \bigcap y \mid}
#'                      {\mid x \bigcup y \mid} =
#'                 \frac{\mid x \bigcap y \mid}
#'                      {\mid x \mid + \mid y \mid -
#'                       \mid x \bigcap y \mid} }
#'
#' @param x [ivl_df]
#' @param y [ivl_df]
#'
#' @template stats
#'
#' @family interval statistics
#'
#' @return
#' tibble with the following columns:
#'
#'   - `len_i` length of the intersection in base-pairs
#'   - `len_u` length of the union in base-pairs
#'   - `jaccard` value of jaccard statistic
#'   - `n_int` number of intersecting intervals between `x` and `y`
#'
#' If inputs are grouped, the return value will contain one set of values per group.
#'
#' @seealso
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html}
#'
#' @examples
#' genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))
#'
#' x <- bed_random(genome, seed = 1010486)
#' y <- bed_random(genome, seed = 9203911)
#'
#' bed_jaccard(x, y)
#'
#' # calculate jaccard per chromosome
#' bed_jaccard(
#'   dplyr::group_by(x, chrom),
#'   dplyr::group_by(y, chrom)
#' )
#'
#' @export
bed_jaccard <- function(x, y) {
  check_required(x)
  check_required(y)

  x <- check_interval(x)
  y <- check_interval(y)

  groups_shared <- shared_groups(x, y)

  x <- bed_merge(x)
  y <- bed_merge(y)

  # TODO: revisit when min_overlap default changes to 1L
  res_intersect <- bed_intersect(x, y, min_overlap = 0L)

  if (!is.null(groups_shared)) {
    x <- group_by(x, !!!syms(groups_shared))
    y <- group_by(y, !!!syms(groups_shared))

    res_intersect <- group_by(res_intersect, !!!syms(groups_shared))
  }

  res_intersect <- summarize(
    res_intersect,
    sum_overlap = sum(as.numeric(.data[[".overlap"]])),
    n_int = as.numeric(n())
  )

  res_x <- mutate(x, .size = .data[["end"]] - .data[["start"]])
  res_x <- summarize(res_x, sum_x = sum(as.numeric(.data[[".size"]])))

  res_y <- mutate(y, .size = .data[["end"]] - .data[["start"]])
  res_y <- summarize(res_y, sum_y = sum(as.numeric(.data[[".size"]])))

  if (!is.null(groups_shared)) {
    res <- left_join(res_intersect, res_x, by = as.character(groups_shared))
    res <- left_join(res, res_y, by = as.character(groups_shared))

    res <- mutate(res, sum_xy = .data[["sum_x"]] + .data[["sum_y"]])
    group_cols <- select(res, !!!syms(groups_shared))

    res <- transmute(
      res,
      len_i = .data[["sum_overlap"]],
      len_u = .data[["sum_xy"]],
      jaccard = .data[["sum_overlap"]] /
        (.data[["sum_xy"]] - .data[["sum_overlap"]]),
      n = .data[["n_int"]]
    )

    res <- bind_cols(group_cols, res)
  } else {
    n_i <- res_intersect$sum_overlap
    n_u <- res_x$sum_x + res_y$sum_y

    jaccard <- n_i / (n_u - n_i)

    res <- tibble(
      len_i = n_i,
      len_u = n_u,
      jaccard = jaccard,
      n = res_intersect$n_int
    )
  }

  res
}
