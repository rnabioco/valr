#' Compute the absolute distance between query and reference intervals
#' 
#' @details \code{bed_absdist()} computes the absolute distance between the midpoint
#'   of query intervals and the closest midpoints of a set of reference
#'   intervals. Absolute distances are scaled by the inter-reference gap for the
#'   chromosome as follows. For \code{Q} total query points and \code{R} reference points
#'   on a chromosome, scale the distance for each query point \code{i} to the closest
#'   reference point by the inter-reference gap for each chromosome. If the chromosome for
#'    a supplied  \code{x} interval has no matching  \code{y} chromosome, the \code{absdist} will
#'    be reported as an \code{NA}.
#'
#'  
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
  
  # find minimum shared groups
  xy_groups <- shared_groups(y, x)
  
  x <- group_by_(x, .dots = c("chrom", xy_groups))
  y <- group_by_(y, .dots = c("chrom", xy_groups))

  res <- absdist_impl(x, y)
  
  # calculate reference sizes
  if (!is.null(xy_groups)){
    xy_groups <- purrr::map_chr(xy_groups, as.character)
  }

  genome <- filter(genome, genome$chrom %in% res$chrom)
  genome <- inner_join(genome, attributes(y)$labels, by = c("chrom"))
  
  ref_points <- summarize(y, ref_points = n())
  genome <- inner_join(genome, ref_points, by = c("chrom", xy_groups))
  
  genome <- mutate(genome, 
                   ref_gap = ref_points / size)
  genome <- select(genome, -size, -ref_points)
  
  #calculate scaled reference sizes
  res <- full_join(res, genome, by = c("chrom", xy_groups))
  res <- mutate(res, scaled_absdist =  absdist * ref_gap)
  res <- select(res, -ref_gap)
  
  #report back original x intervals not found 
  missing_x <- anti_join(x, res, by = c("chrom", xy_groups))
  missing_x <- ungroup(missing_x)
  missing_x <- mutate(missing_x, 
                absdist = NA,
                scaled_absdist = NA)
  res <- bind_rows(res, missing_x)
  res
  
}


