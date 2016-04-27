#' Divide intervals into new intervals with labels.
#' 
#' @param bed_df BED data in \code{dplyr::tbl_df} format
#' @param genome genome file with chromosome sizes
#' @param win_size divide intervals into fixed-size windows
#' @param step_size size to step before next window
#' @param num_windows divide intervals to fixed number of windows
#' @param reverse reverse window numbers
#' 
#' @param win_id one of \code{name}, \code{num}, \code{namenum} (default \code{name})
#' 
#' @return \code{data_frame} with \code{.win_id} column
#' 
#' @examples 
#' genome <- tibble::frame_data(
#'  ~chrom, ~size,
#'  "chr1", 5000,
#'  "chr2", 400
#' )
#' 
#' bed_df <- tibble::frame_data(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 100,    200,  'A',   '.',    '+',
#'   "chr2", 300,    350,  'B',   '.',    '-'
#' ) 
#' 
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
#' # named intervals (name)
#' bed_makewindows(bed_df, genome, win_size = 10, win_id = 'name')
#' 
#' # named intervals (num)
#' bed_makewindows(bed_df, genome, win_size = 10, win_id = 'num')
#' 
#' # named intervals (reversed num)
#' bed_makewindows(bed_df, genome, win_size = 10, win_id = 'num', reverse = TRUE)
#' 
#' # named intervals (namenum)
#' bed_makewindows(bed_df, genome, win_size = 10, win_id = 'namenum', TRUE)
#' 
#' # named intervals (reversed namenum)
#' bed_makewindows(bed_df, genome, win_size = 10, win_id = 'namenum', reverse = TRUE)
#' 
#' small_bed_df <- tibble::frame_data(
#'   ~chrom, ~start, ~end,
#'   "chr1", 100,    200, 
#'   "chr2", 300,    350
#' )
#'
#' # named intervals (for BED3 tbls)
#' bed_makewindows(small_bed_df, genome, win_size = 10, win_id = 'name')
#' 
#' @export
bed_makewindows <- function(bed_df, genome, win_size = 0,
                            step_size = 0, num_windows = 0,
                            reverse = FALSE, win_id = NULL) {
 
  assert_that(win_size > 0 || num_windows > 0)
  assert_that(step_size >= 0)
  
  win_id <- match.arg(win_id, c('name', 'num', 'namenum'))
  
  res <- bed_df %>%
    by_row(split_interval, genome, win_size,
           step_size, num_windows, win_id,
           reverse, .collate = 'rows',
           .labels = FALSE) %>%
    select(-.row)
    
  res 
}

split_interval <- function(interval, genome, win_size, step_size,
                           num_windows, win_id, reverse) {
  
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
           .win_num = row_number()) %>%
    filter(.start != .end & .end <= chrom_size) %>%
    mutate(start = .start, end = .end) %>%
    select(-.start, -.end)
 
  # add win_id 
  if (reverse) {
    res <- mutate(res, .win_num = rank(-.win_num))
  } else {
    res <- mutate(res, .win_num = rank(.win_num))
  } 
  
  if (win_id == 'name') {
    if ( ! 'name' %in% colnames(res) ) {
      res <- mutate(res, .win_id = str_c(chrom, ':', start, '-', end)) %>% select(-.win_num) 
    } else {
      res <- mutate(res, .win_id = name) %>% select(-.win_num)
    }
  } else if (win_id == 'num') {
    res <- rename(res, .win_id = .win_num)
  } else if (win_id == 'namenum') {
    res <- mutate(res, .win_id = str_c(name, "_", .win_num)) %>% select(-.win_num)
  }
  
  res
}
