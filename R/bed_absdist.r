#' Compute absolute distances between intervals.
#'
#' Computes the absolute distance between the midpoint of each `x` interval and
#' the midpoints of each closest `y` interval.
#'
#' @details Absolute distances are scaled by the inter-reference gap for the
#'   chromosome as follows. For `Q` query points and `R` reference
#'   points on a chromosome, scale the distance for each query point `i` to
#'   the closest reference point by the inter-reference gap for each chromosome.
#'   If an `x` interval has no matching `y` chromosome,
#'   `.absdist` is `NA`.
#'
#'   \deqn{d_i(x,y) = min_k(|q_i - r_k|)\frac{R}{Length\ of\ chromosome}}
#'
#'   Both absolute and scaled distances are reported as `.absdist` and
#'   `.absdist_scaled`.
#'
#' @param x [tbl_interval()]
#' @param y [tbl_interval()]
#' @param genome [tbl_genome()]
#'
#' @return
#' [tbl_interval()] with `.absdist` and `.absdist_scaled` columns.
#'
#' @template stats
#'
#' @family interval statistics
#'
#' @seealso
#' \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002529}
#'
#' @examples
#' genome <- read_genome(valr_example('hg19.chrom.sizes.gz'))
#'
#' x <- bed_random(genome, seed = 1010486)
#' y <- bed_random(genome, seed = 9203911)
#'
#' bed_absdist(x, y, genome)
#'
#' @export
bed_absdist <- function(x, y, genome) {
  if (!is.tbl_interval(x)) x <- as.tbl_interval(x)
  if (!is.tbl_interval(y)) y <- as.tbl_interval(y)
  if (!is.tbl_genome(genome)) genome <- as.tbl_genome(genome)

  # establish grouping with shared groups (and chrom)
  groups_xy <- shared_groups(x, y)
  groups_xy <- unique(as.character(c("chrom", groups_xy)))
  groups_vars <- rlang::syms(groups_xy)

  # type convert grouping factors to characters if necessary and ungroup
  x <- convert_factors(x, groups_xy)
  y <- convert_factors(y, groups_xy)

  x <- group_by(x, !!! groups_vars)
  y <- group_by(y, !!! groups_vars)

  if (utils::packageVersion("dplyr") < "0.7.99.9000"){
    x_cpp <- update_groups(x)
    y_cpp <- update_groups(y)
    grp_indexes <- shared_group_indexes(x_cpp, y_cpp)
    res <- dist_impl(x_cpp, y_cpp,
                     grp_indexes$x,
                     grp_indexes$y,
                     distcalc = "absdist")
  } else {
    grp_indexes <- shared_group_indexes(x, y)
    res <- dist_impl(x, y,
                     grp_indexes$x, grp_indexes$y,
                     distcalc = "absdist")
  }

  # convert groups_xy to character vector
  if (!is.null(groups_xy)) {
    groups_xy <- as.character(groups_xy)
  }

  # calculate reference sizes
  genome <- filter(genome, genome$chrom %in% res$chrom)
  genome <- inner_join(genome, get_labels(y), by = c("chrom"))

  ref_points <- summarize(y, .ref_points = n())
  genome <- inner_join(genome, ref_points, by = c("chrom", groups_xy))

  genome <- mutate(genome, .ref_gap = .ref_points / size)
  genome <- select(genome, -size, -.ref_points)

  # calculate scaled reference sizes
  res <- full_join(res, genome, by = c("chrom", groups_xy))
  res <- mutate(res, .absdist_scaled = .absdist * .ref_gap)
  res <- select(res, -.ref_gap)

  # report back original x intervals not found
  x_missing <- anti_join(x, res, by = c("chrom", groups_xy))
  x_missing <- ungroup(x_missing)
  x_missing <- mutate(x_missing, .absdist = NA, .absdist_scaled = NA)
  res <- bind_rows(res, x_missing)

  res <- bed_sort(res)

  res
}
