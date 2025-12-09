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
#' @param x [ivl_df]
#' @param y [ivl_df]
#' @param genome [genome_df]
#'
#' @return
#' [ivl_df] with `.absdist` and `.absdist_scaled` columns.
#'
#' @template stats
#'
#' @family interval statistics
#'
#' @seealso
#' \url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002529}
#'
#' @examples
#' genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))
#'
#' x <- bed_random(genome, seed = 1010486)
#' y <- bed_random(genome, seed = 9203911)
#'
#' bed_absdist(x, y, genome)
#'
#' @export
bed_absdist <- function(x, y, genome) {
  check_required(x)
  check_required(y)
  check_required(genome)

  x <- check_interval(x)
  y <- check_interval(y)
  genome <- check_genome(genome)

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
  res <- dist_impl(x, y, grp_indexes$x, grp_indexes$y, distcalc = "absdist")
  res <- tibble::as_tibble(res)

  # convert groups_xy to character vector
  if (!is.null(groups_xy)) {
    groups_xy <- as.character(groups_xy)
  }

  # calculate reference sizes
  genome <- filter(genome, .data[["chrom"]] %in% res$chrom)

  if (utils::packageVersion("dplyr") > "1.0.10") {
    genome <- inner_join(
      genome,
      get_labels(y),
      by = c("chrom"),
      multiple = "all"
    )
  } else {
    genome <- inner_join(genome, get_labels(y), by = c("chrom"))
  }

  ref_points <- summarize(y, .ref_points = n())
  genome <- inner_join(genome, ref_points, by = c(groups_xy))

  genome <- mutate(genome, .ref_gap = .data[[".ref_points"]] / .data[["size"]])
  genome <- select(genome, -all_of(c("size", ".ref_points")))

  # calculate scaled reference sizes
  res <- full_join(res, genome, by = c(groups_xy))
  res <- mutate(
    res,
    .absdist_scaled = .data[[".absdist"]] * .data[[".ref_gap"]]
  )
  res <- select(res, -all_of(".ref_gap"))

  # report back original x intervals not found
  x_missing <- anti_join(x, res, by = c(groups_xy))
  x_missing <- ungroup(x_missing)
  x_missing <- mutate(x_missing, .absdist = NA, .absdist_scaled = NA)
  res <- bind_rows(res, x_missing)

  res <- bed_sort(res)

  res
}
