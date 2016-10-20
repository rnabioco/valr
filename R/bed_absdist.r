#' Compute the absolute distance between query and reference intervals
#' 
#' @details \code{bed_absdist()} computes the absolute distance between the midpoint
#'   of query intervals and the closest midpoints of a set of reference
#'   intervals. Absolute distances are scaled by the inter-reference gap for the
#'   chromosome as follows. For \code{Q} total query points and \code{R} reference points
#'   on a chromosome, scale the distance for each query point \code{i} to the closest
#'   reference point by the inter-reference gap for each chromosome.
#'
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param genome genome tbl
#' 
#' @return \code{data_frame}
#'
#' @family interval-stats
#' @seealso \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002529}
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
#' 
#' @export

bed_absdist <- function(x, y, genome) {
  
  x <- group_by(x, chrom, add = TRUE)
  y <- group_by(y, chrom, add = TRUE)
  
  # make sure that y tbl has same grouping as x tbl
  y <- set_groups(y, x)
  
  res <- absdist_impl(x, y)
  
  # calculate reference sizes
  y_chroms <- unique(y$chrom)
  y_groups <- purrr::map_chr(groups(y), as.character)
  
  genome <- filter(genome, genome$chrom %in% y_chroms)
  genome <- inner_join(genome, attributes(y)$labels, by = c("chrom"))
  
  ref_points <- summarize(y, ref_points = n())
  genome <- inner_join(genome, ref_points, by = c("chrom", y_groups))
  
  genome <- mutate(genome, 
                   ref_gap = ref_points / size)
  genome <- select(genome, -size, -ref_points)
  
  #calculate scaled reference sizes
  res <- full_join(res, genome, by = c("chrom", y_groups))
  res <- mutate(res, scaled_absdist =  absdist * ref_gap)
  res <- select(res, -ref_gap)
  res
  
}


