#' Compute absolute distances between intervals.
#'
#' Computes the absolute distance between the midpoints of `x` intervals and
#' the closest midpoints of `y` intervals.
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
#' @return [tbl_interval()] with `.absdist` and `.absdist_scaled`
#'   columns.
#'
#' @template stats
#'
#' @family interval statistics
#'
#' @seealso
#' \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002529}
#'
#' @examples
#' x <- trbl_interval(
#'   ~chrom,   ~start,    ~end,
#'   "chr1",    75,       125
#' )
#'
#' y <- trbl_interval(
#'   ~chrom,   ~start,    ~end,
#'   "chr1",    50,       100,
#'   "chr1",    100,       150
#' )
#'
#' genome <- trbl_genome(
#'   ~chrom, ~size,
#'   "chr1", 500,
#'   "chr2", 1000
#' )
#'
#' bed_absdist(x, y, genome)
#'
#' @export
bed_absdist <- function(x, y, genome) {

  if (!is.tbl_interval(x)) x <- tbl_interval(x)
  if (!is.tbl_interval(y)) y <- tbl_interval(y)
  if (!is.tbl_genome(genome)) genome <- tbl_genome(genome)

  # find minimum shared groups
  groups_xy <- shared_groups(y, x)

  x <- group_by_(x, .dots = c("chrom", groups_xy))
  y <- group_by_(y, .dots = c("chrom", groups_xy))

  res <- absdist_impl(x, y)

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

  res <- arrange(res, chrom, start)

  res
}
