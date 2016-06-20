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
    assert_that(strand %in% colnames(x) & strand %in% colnames(y))
    x <- group_by(x, strand, add = TRUE)
    y <- group_by(y, strand, add = TRUE)
  }
  
  # number of intervals 
  n_x <- summarize(x, n_x = n())
  n_y <- summarize(y, n_y = n()) 

  # number of intersections 
  n_i <- bed_intersect(x, y) %>%
    group_by(chrom) %>%
    unique() %>%
    summarize(n_i = n())
  
  # estimate number of possible intervals as chrom_size / mean(interval_size)
  xy <- bind_rows(list(x, y))
  n_xy <- suppressMessages(
    xy %>%
    mutate(.size = end - start) %>%
    group_by(chrom) %>%
    summarize(size_mean = mean(.size)) %>%
    left_join(genome) %>%
    mutate(int_est = size / size_mean))

  res <- suppressMessages(
    n_x %>% 
    left_join(n_y) %>% 
    left_join(n_i) %>% 
    left_join(n_xy) %>% 
    group_by(chrom) %>%
    mutate(p.value = phyper(n_i - (n_x - n_y),
                            n_x + n_y,
                            int_est - (n_x + n_y),
                            n_i, lower.tail = FALSE)))
  
  res
  
}
