#' Divide intervals into new sub-intervals ("windows").
#' 
#' @param x tbl of intervals
#' @param genome genome file with chromosome sizes
#' @param win_size divide intervals into fixed-size windows
#' @param step_size size to step before next window
#' @param num_win divide intervals to fixed number of windows
#' @param reverse reverse window numbers
#'   
#' @note The \code{name} and \code{win_id} columns can be used to create new 
#'   interval names (see 'namenum' example below) or in subsequent 
#'   \code{group_by} operations (see vignette).
#'   
#' @return \code{data_frame} with \code{win_id} column that contains a numeric 
#'   identifier for the window.
#'   
#' @examples 
#' genome <- tibble::frame_data(
#'  ~chrom, ~size,
#'  "chr1", 5000,
#'  "chr2", 400
#' )
#' 
#' x <- tibble::frame_data(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 100,    200,  'A',   '.',    '+',
#'   "chr2", 300,    350,  'B',   '.',    '-'
#' ) 
#' 
#' 
#' # Fixed number of windows 
#' bed_makewindows(x, genome, num_win = 10)
#' 
#' # Fixed window size
#' bed_makewindows(x, genome, win_size = 10)
#' 
#' # Fixed window size with overlaps
#' bed_makewindows(x, genome, win_size = 10, step_size = 5)
#' 
#' # reverse win_id
#' bed_makewindows(x, genome, win_size = 10, reverse = TRUE)
#' 
#' # bedtools 'namenum'
#' wins <- bed_makewindows(x, genome, win_size = 10)
#' dplyr::mutate(wins, namenum = stringr::str_c(name, '_', win_id))
#' 
#' @export
bed_makewindows <- function(x, genome, win_size = 0,
                            step_size = 0, num_win = 0,
                            reverse = FALSE) {
 
  assert_that(win_size > 0 || num_win > 0)
  assert_that(step_size >= 0)
  
  res <- x %>%
    by_row(split_interval, genome, win_size,
           step_size, num_win, 
           reverse, .collate = 'rows',
           .labels = FALSE) %>%
    select(-.row)
    
  res 
}

#' @describeIn bed_makewindows Helper to divide interval into labeled sub-intervals.
#' 
#' @param interval row of data frame
#' 
split_interval <- function(interval, genome, win_size, step_size,
                           num_win, reverse) {
  
  # get size of chrom for coord check later
  cur_chrom <- interval$chrom
  chrom_size <- genome[genome$chrom == cur_chrom,]$size
  
  if (num_win > 0) {
    win_size <- round((interval$end - interval$start) / num_win)
  } 
   
  res <- interval %>%
    transform(.start = seq(from = start, to = end,
                           by = win_size - step_size)) %>%
    mutate(.end = ifelse(.start + win_size < end,
                         .start + win_size, end),
           .win_num = row_number()) %>%
    filter(.start != .end & .end <= chrom_size) %>%
    mutate(start = .start, end = .end) %>%
    select(-.start, -.end)
 
  # add win_id 
  if (reverse) {
    res <- mutate(res, win_id = rank(-.win_num))
  } else {
    res <- mutate(res, win_id = rank(.win_num))
  }
  
  res <- select(res, -.win_num)
  
  res
}
