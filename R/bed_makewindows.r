#' Divide intervals into new sub-intervals ("windows").
#' 
#' @param x tbl of intervals
#' @param genome genome file with chromosome sizes
#' @param win_size divide intervals into fixed-size windows
#' @param step_size size to step before next window
#' @param num_win divide intervals to fixed number of windows
#' @param reverse reverse window numbers
#'   
#' @note The \code{name} and \code{.win_id} columns can be used to create new 
#'   interval names (see 'namenum' example below) or in subsequent 
#'   \code{group_by} operations (see vignette).
#' 
#' @family utils  
#' @return \code{data_frame} with \code{.win_id} column that contains a numeric 
#'   identifier for the window.
#'   
#' @examples 
#' genome <- tibble::tribble(
#'  ~chrom, ~size,
#'  "chr1", 200
#' )
#' 
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 100,    200,  'A',   '.',    '+'
#' )
#' 
#' bed_glyph(bed_makewindows(x, genome, num_win = 10), label = '.win_id') 
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
#' dplyr::mutate(wins, namenum = stringr::str_c(name, '_', .win_id))
#' 
#' @export
bed_makewindows <- function(x, genome, win_size = 0,
                            step_size = 0, num_win = 0,
                            reverse = FALSE) {

  if (win_size == 0 && num_win == 0)
    stop('specify either `win_size` or `num_win`', call. = FALSE)
  
  res <- purrr::by_row(x, split_interval, genome, win_size,
           step_size, num_win, 
           reverse, .collate = 'rows',
           .labels = FALSE)
  
  res <- select(res, -.row)
    
  res 
}

#' @param interval row of data frame
#' @noRd
split_interval <- function(interval, genome, win_size, step_size,
                           num_win, reverse) {
  
  # get size of chrom for coord check later
  cur_chrom <- interval$chrom
  chrom_size <- genome[genome$chrom == cur_chrom,]$size
  
  if (num_win > 0) {
    win_size <- round((interval$end - interval$start) / num_win)
  } 
   
  res <- transform(interval, .start = seq(from = start, to = end,
                   by = win_size - step_size))
  
  res <- mutate(res, .end = ifelse(.start + win_size < end,
                                   .start + win_size, end),
                     .win_num = row_number())
  res <- filter(res, .start != .end & .end <= chrom_size)
  res <- mutate(res, start = .start, end = .end)
  res <- select(res, -.start, -.end)
 
  # add .win_id column
  if (reverse) {
    res <- mutate(res, .win_id = rank(-.win_num))
  } else {
    res <- mutate(res, .win_id = rank(.win_num))
  }
  
  res <- select(res, -.win_num)
  
  res
}
