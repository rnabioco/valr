#' Fisher's test on number of shared and unique intervals.
#' 
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param genome tbl of chrom sizes
#' 
#' @family interval-stats
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/fisher.html}
#' 
#' @return \code{data_frame}
#' 
#' @examples 
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 10,     20,
#'   "chr1", 30,     40,
#'   "chr1", 51,     52
#' )
#' 
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 15,     25,
#'   "chr1", 51,     52
#' )
#' 
#' genome <- tibble::tribble(
#'   ~chrom, ~size,
#'   "chr1", 500
#' )
#' 
#' bed_fisher(x, y, genome)
#' 
#' @export
bed_fisher <- function(x, y, genome) {

  # number of intervals 
  n_x <- nrow(x) 
  n_y <- nrow(y)
 
  # union of intervals (i.e. total bases covered)
  union_x <- interval_union(x)
  union_y <- interval_union(y)
  
  # mean interval sizes
  mean_x <- union_x / n_x
  mean_y <- union_y / n_y
  
  # heuristic from bedtools fisher.cpp
  mean_total <- mean_x + mean_y
  
  # number of intersections (`n11` in fisher.cpp)
  isect <- bed_intersect(x, y)
  n_i <- nrow(isect)

  # x, not y (`n12` in fisher.cpp)
  n_x_only <- max(0, n_x - n_i)
  # y, not x (`n21` in fisher.cpp)
  n_y_only <- max(0, n_y - n_i)

  genome_size = sum(as.numeric(genome$size))
  
  # estimated total intervals (`n22_full`) 
  total_est <- round(max(n_i + n_x_only + n_y_only,
                     genome_size / mean_total ))
  
  # estimate n for neither x nor y (`n22`)
  not_est <- total_est - n_i - n_x_only - n_y_only
  
  fisher_mat <- matrix(c(n_i, n_x_only, n_y_only, not_est),
                       nrow = 2,
                       dimnames = list('in y?' = c('yes', 'no'),
                                       'in x?' = c('yes', 'no')))
 
  stat <- stats::fisher.test(fisher_mat)
  broom::tidy(stat)
  
}

#' @noRd
interval_union <- function(x) {
  res <- bed_merge(x)
  res <- mutate(res, .size = end - start)
  
  sum(res$.size) 
}
