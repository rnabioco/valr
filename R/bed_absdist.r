#' Compute absolute distances between intervals.
#' 
#' @details \code{bed_absdist()} computes the absolute distance between the 
#'   midpoint of query intervals and the closest midpoints of a set of reference
#'   intervals.
#'   
#'   Absolute distances are scaled by the inter-reference gap for the chromosome
#'   as follows. For \code{Q} total query points and \code{R} reference points
#'   on a chromosome, scale the distance for each query point \code{i} to the
#'   closest reference point by the inter-reference gap for each chromosome. If
#'   the chromosome for a supplied  \code{x} interval has no matching  \code{y}
#'   chromosome, the \code{absdist} will be reported as an \code{NA}.
#'   
#'   \deqn{d_i(x,y) = min_k(|q_i - r_k|)\frac{R}{Length\ of\ chromosome}}
#'   
#'   By default both absolute and scaled distances are reported as \code{.absdist} and
#'   \code{.absdist_scaled} respectively.
#'   
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param genome genome tbl
#'   
#' @return \code{data_frame} with \code{.absdist} and \code{.absdist_scaled}
#'   columns.
#'   
#' @family interval-stats
#' @seealso 
#' \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002529}
#' 
#' @examples
#' x <- tibble::frame_data(
#' ~chrom,   ~start,    ~end,
#' "chr1",    75,       125
#'   )
#' 
#' y <- tibble::frame_data(
#'   ~chrom,   ~start,    ~end,
#'   "chr1",    50,       100,
#'   "chr1",    100,       150
#'   )
#'   
#' genome <- tibble::frame_data(
#'   ~chrom, ~size,
#'   "chr1", 500,
#'   "chr2", 1000
#' )
#' 
#' bed_absdist(x, y, genome)
#' 
#' @export
bed_absdist <- function(x, y, genome) {
  
  # find minimum shared groups
  groups_xy <- shared_groups(y, x)
  
  x <- group_by_(x, .dots = c("chrom", groups_xy))
  y <- group_by_(y, .dots = c("chrom", groups_xy))

  res <- absdist_impl(x, y)
  
  # convert groups_xy to character vector
  if (!is.null(groups_xy)){
    groups_xy <- purrr::map_chr(groups_xy, as.character)
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
