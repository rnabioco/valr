#' Fisher's test on number of shared and unique intervals.
#' 
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param genome tbl of genome intervals
#' @param strand group intervals by strand
#' 
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/fisher.html}
#' 
#' @return \code{data_frame}
#' 
#' @examples 
#' x <- tibble::frame_data(
#'   ~chrom, ~start, ~end,
#'   "chr1", 10,     20,
#'   "chr1", 30,     40,
#'   "chr1", 51,     52,
#'   "chr2", 10,     40
#' )
#' 
#' y <- tibble::frame_data(
#'   ~chrom, ~start, ~end,
#'   "chr1", 15,     25,
#'   "chr1", 51,     52,
#'   "chr2", 35,     60
#' )
#' 
#' genome <- tibble::frame_data(
#'   ~chrom, ~size,
#'   "chr1", 500,
#'   "chr2", 1000
#' )
#' 
#' bed_fisher(x, y, genome)
#' 
#' @export
bed_fisher <- function(x, y, genome, strand = FALSE) {

  x <- group_by(x, chrom)
  y <- group_by(y, chrom)
  
  if (strand) {
    if (! 'strand' %in% colnames(x) || ! 'strand' %in% colnames(y))
      stop('expect `strand` column in input', call. = FALSE)
    
    x <- group_by(x, strand, add = TRUE)
    y <- group_by(y, strand, add = TRUE)
  }
  
  # number of intervals 
  n_x <- summarize(x, n_x = n())
  n_y <- summarize(y, n_y = n()) 

  # number of intersections 
  n_i <- bed_intersect(x, y)
  n_i <- group_by(n_i, chrom)
  n_i <- unique(n_i)
  n_i <- summarize(n_i, n_i = n())
  
  # estimate number of possible intervals as chrom_size / mean(interval_size)
  xy <- bind_rows(list(x, y))
  n_xy <- mutate(xy, .size = end - start)
  n_xy <- group_by(n_xy, chrom)
  n_xy <- summarize(n_xy, size_mean = mean(.size))
  n_xy <- left_join(n_xy, genome)
  n_xy <- mutate(n_xy, int_est = size / size_mean)

  res <- left_join(n_x, n_y)
  res <- left_join(res, n_i)
  res <- left_join(res, n_xy)
  res <- group_by(res, chrom)
  res <- mutate(res,
                p.value = phyper(n_i - (n_x - n_y),
                                 n_x + n_y,
                                 int_est - (n_x + n_y),
                                 n_i, lower.tail = FALSE))
  
  res
  
}
