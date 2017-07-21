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

  # find minimum shared groups
  groups_xy <- shared_groups(y, x)
  groups_vars <- rlang::syms(c("chrom", groups_xy))

  x <- group_by(x, !!! groups_vars)
  y <- group_by(y, !!! groups_vars)

  res <- dist_impl(x, y, distcalc = "absdist")

  # convert groups_xy to character vector
  if (!is.null(groups_xy)){
    groups_xy <- as.character(groups_xy)
  }

  # calculate reference sizes
  genome <- filter(genome, genome$chrom %in% res$chrom)
  genome <- inner_join(genome, attributes(y)$labels, by = c("chrom"))

  ref_points <- summarize(y, .ref_points = n())
  genome <- inner_join(genome, ref_points, by = c("chrom", groups_xy))

  genome <- mutate(genome, .ref_gap = .ref_points / size)
  genome <- select(genome, -size, -.ref_points)

  #calculate scaled reference sizes
  res <- full_join(res, genome, by = c("chrom", groups_xy))
  res <- mutate(res, .absdist_scaled = .absdist * .ref_gap)
  res <- select(res, -.ref_gap)

  #report back original x intervals not found
  x_missing <- anti_join(x, res, by = c("chrom", groups_xy))
  x_missing <- ungroup(x_missing)
  x_missing <- mutate(x_missing, .absdist = NA, .absdist_scaled = NA)
  res <- bind_rows(res, x_missing)

  res <- bed_sort(res)

  res
}
