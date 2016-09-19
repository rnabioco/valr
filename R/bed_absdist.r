#' Compute the absolute distance between query intervals and reference intervals
#' 
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param genome genome tbl
#' 
#' @return \code{data_frame}
#' 
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
  
  if ( ! is_sorted(x) )
    x <- bed_sort(x)
  if ( ! is_sorted(y) )
    y <- bed_sort(y)
  
  x <- group_by(x, chrom, add = TRUE)
  y <- group_by(y, chrom, add = TRUE)
  
  res <- absdist_impl(x, y)
  
  # calculate reference sizes
  y_chroms <- unique(y$chrom)
  genome <- filter(genome, genome$chrom %in% y_chroms)
  genome <- mutate(genome, 
                   ref_gap = group_size(y),
                   ref_gap = ref_gap / size)
  genome <- select(genome, -size)
  
  #calculate scaled reference sizes
  res <- full_join(res, genome, by = c("chrom"))
  res <- mutate(res, scaled_absdist =  absdist * ref_gap)
  res <- select(res, -ref_gap)
  res
  
}


