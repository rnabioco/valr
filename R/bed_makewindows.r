#' nDivide intervals into new intervals with labels
#' 
#' @param bed_df BED data in \code{dplyr::tbl_df} format
#' @param genome genome file with chromosome sizes
#' @param win_size divide intervals into fixed-size windows
#' @param step_size size to step before next window
#' @param num_windows divide intervals to fixed number of windows
#' @param reverse reverse window numbers?
#' 
#' @param win_names one of 'name', 'num', 'namenum'
#' 
#' @return \code{data.frame} with \code{win_num} column
#' 
#' @examples 
#' genome <- dplyr::tibble(
#'  ~chrom, ~size,
#'  "chr1", 5000,
#'  "chr2", 400
#' )
#' 
#' bed_df <- dplyr::tibble(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 100,    200,  'A',   '.',    '+',
#'   "chr2", 300,    350,  'B',   '.',    '-'
#' ) 
#' 
#' # Fixed number of windows 
#' bed_makewindows(bed_df, genome, num_windows = 10)
#' 
#' # Fixed window size
#' bed_makewindows(bed_df, genome, win_size = 10)
#' 
#' # Fixed window size with overlaps
#' bed_makewindows(bed_df, genome, win_size = 10, step_size = 5)
#' 
#' # Reversed window numbering
#' bed_makewindows(bed_df, genome, win_size = 10, reverse = TRUE)
#'  
#' @export
bed_makewindows <- function(bed_df, genome, win_size = 0,
                            step_size = 0, num_windows = 0,
                            reverse = FALSE) {
 
  assert_that(win_size > 0 || num_windows > 0)
  assert_that(step_size >= 0)
  
  res <- bed_df %>%
    by_row(bed_makewindows_, genome, win_size,
           step_size, num_windows,
           reverse, .collate = 'rows',
           .labels = FALSE) %>%
    select(-.row)
    
  res 
}

#' @rdname bed_makewindows
#' @export
bed_makewindows_ <- function(interval, genome, win_size, step_size,
                             num_windows, reverse) {
  
  # get size of chrom for coord check later
  cur_chrom <- interval$chrom
  chrom_size <- genome[genome$chrom == cur_chrom,]$size
  
  if (num_windows > 0) {
    win_size <- round((interval$end - interval$start) / num_windows)
  } 
   
  res <- interval %>%
    transform(.start = seq(from = start, to = end,
                           by = win_size - step_size)) %>%
    mutate(.end = ifelse(.start + win_size < end,
                         .start + win_size, end),
           win_num = row_number()) %>%
    filter(.start != .end & .end <= chrom_size) %>%
    mutate(start = .start, end = .end) %>%
    select(-.start, -.end)
  
  if (reverse) {
    res <- res %>%
       mutate(win_num = rank(-win_num)) 
  }
  
  res
}
