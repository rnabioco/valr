#' Create flanks from input intervals.
#' 
#' @param x tbl of intervals
#' @param genome tbl of chrom sizes
#' @param both number of bases on both sizes 
#' @param left number of bases on left side
#' @param right number of bases on right side
#' @param strand define \code{left} and \code{right} based on strand
#' @param fraction define flanks based on fraction of interval length
#' @param trim adjust coordinates for out-of-bounds intervals
#' 
#' @return \code{data_frame}
#' 
#' @seealso
#'   \url{http://bedtools.readthedocs.org/en/latest/content/tools/flank.html}
#' 
#' @examples 
#' genome <- tibble::frame_data(
#'  ~chrom, ~size,
#'  "chr1", 5000
#' )
#' 
#' x <- tibble::frame_data(
#'  ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'  "chr1", 500,    1000, '.',   '.',    '+',
#'  "chr1", 1000,   1500, '.',   '.',    '-'
#' )
#' 
#' bed_flank(x, genome, left = 100)
#' bed_flank(x, genome, right = 100)
#' bed_flank(x, genome, both = 100)
#'
#' bed_flank(x, genome, both = 0.5, fraction=TRUE)
#' 
#' @export
bed_flank <- function(x, genome, both = 0, left = 0,
                      right = 0, fraction = FALSE,
                      strand = FALSE, trim = FALSE) {

  assert_that(both > 0 || left > 0 || right > 0)
  
  if (both != 0 && (left != 0 || right != 0)) {
    stop('ambiguous side spec for bed_flank')
  } 
  
  if (both) {
    left <- both 
    right <- both
  }
  
  if (strand) {
    if (fraction) {
      res <- x %>%
        mutate(.interval_size = end - start,
               left_start = ifelse(strand == '+', 
                                start - round( left * .interval_size ),
                                end),
               left_end = ifelse(strand == '+',
                              start,
                              end + round( left * .interval_size )),
               right_start = ifelse(strand == '+',
                                end,
                                start - round( right * .interval_size )),
               right_end = ifelse(strand == '+',
                              end + round( right * .interval_size ),
                              start)) %>%
        select(-start, -end, -.interval_size) 
      
    } else {
      res <- x %>%
        mutate(left_start = ifelse(strand == '+', 
                                start - left,
                                end),
               left_end = ifelse(strand == '+',
                              start,
                              end + left),
               right_start = ifelse(strand == '+',
                                end,
                                start - right),
               right_end = ifelse(strand == '+',
                              end + right,
                              start)) %>%
        select(-start, -end)
    }
    
  } else {
    if (fraction) {
      res <- x %>%
        mutate(.interval_size = end - start,
               right_start =  start - round( left * .interval_size ), 
               right_end = start, 
               left_start = end, 
               left_end = end + round( right * .interval_size )) %>%
        select(-start, -end, -.interval_size) 
      
    } else {
      res <- x %>%
        mutate(right_start = start - left, 
               right_end = start, 
               left_start = end, 
               left_end = end + right) %>%
        select(-start, -end)
    }
  }
  
  if (right && !left) {
    res <- res %>%
      mutate(start = right_start,
             end = right_end) %>%
      select(chrom, start, end, everything(), 
             -left_start, -left_end, -right_start, -right_end)
    
  } else if (left && !right) {
    res <- res %>%
      mutate(start = left_start,
             end = left_end) %>%
      select(chrom, start, end, everything(), 
             -left_start, -left_end, -right_start, -right_end)
    
  } else {
    res <- res %>%
      gather(key, value, right_start, right_end, left_start, left_end) %>% 
      separate(key, c('key', 'pos'), sep = '_') %>% 
      spread(key, value) %>%
      select(chrom, start, end, everything(), -pos) 
  }   
  
  if (trim) {
    res <- res %>%
      bound_intervals(genome, trim = T) %>%
      bed_sort()
    
  } else {
    res <- res %>%
      bound_intervals(genome, trim = F) %>%
      bed_sort()
  }
  
  res
}
 
