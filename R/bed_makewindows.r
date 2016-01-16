#'
#' bed_makewindows
#' 
#' Divide intervals into new intervals ("windows") with labels
#' 
#' @param bed_df BED data in \code{dplyr::tbl_df} format
#' @param genome genome file with chromosome sizes
#' @param win_size divide intervals into fixed-size windows
#' @param step_size size to step before next window
#' @param num_wins divide intervals to fixed number of windows
#' @param reverse reverse window numbers?
#' 
#' @param win_names one of 'name', 'num', 'namenum'
#' 
#' @return \code{dplyr::tbl_df} with `win.num` column
#' 
#' @examples
#' 
#' genome <- dplyr::tibble(
#'  ~chrom, ~size,
#'  "chr1", 5000,
#'  "chr2", 400
#' )
#' 
#'  bed_df <- dplyr::tibble(
#'    ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'    "chr1", 100, 200, 'A', '.', '+',
#'    "chr2", 300, 350, 'B', '.', '-'
#'  ) 
#'  
#' @export
bed_makewindows <- function(bed_df, genome, win_size = 0,
                            step_size = 0, num_windows = 0,
                            reverse = FALSE) {
  
  assert_that(win_size >= 0)
  assert_that(step_size >= 0)
  assert_that(num_windows >= 0)
  
  # Warning messages:
  # 1: In bind_rows_(x, .id) : Unequal factor levels: coercing to character
  # 2: In bind_rows_(x, .id) : Unequal factor levels: coercing to character
  
  res <- bed_df %>%
    rowwise() %>%
    do(calculate_intervals(., genome, win_size,
                           step_size, num_windows
                           reverse)) %>%
    ungroup()

  res 
}

calculate_intervals <- function(df, genome, win_size, step_size, num_windows) {
 
  # get size of chrom for coord check later
  cur_chrom <- df$chrom
  chrom_size <- genome[genome$chrom == cur_chrom,]$size
  
  if (num_windows > 0) {
    win_size <- round((df$end - df$start) / num_windows)
  } 
  
  res <- df %>%
    transform(.start = seq(from = start,
                           to = end,
                           by = win_size - step_size)) %>%
    mutate(.end = .start + win_size,
           win_num = row_number()) %>%
    filter(.end <= end & .end <= chrom_size) %>%
    mutate(start = .start,
           end = .end) %>%
    select(-.start, -.end)
  
  if (reverse) {
    res <- res %>%
       mutate(win_num = rank(-win_num)) 
  }
  
  res
}