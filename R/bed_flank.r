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
  #assert_that(fraction >= 0 && fraction <= 1)
  
  if (both != 0 && (left != 0 || right != 0)) {
    stop('ambiguous side spec for bed_flank')
  } 
  
  if (both) {
    if (fraction) {
      res <- x %>%
        mutate(.interval_size = end - start,
               start_r = start - round(both * .interval_size), end_r = start, 
               start_l = end, end_l = end + round(both * .interval_size)) %>%
        select(-start, -end, -.interval_size) 
      
    } else {
      res <- x %>%
        mutate(start_r = start - both, end_r = start, 
               start_l = end, end_l = end + both) %>%
        select(-start, -end)
    }
   
    res <- res %>%
      gather(key, value, start_r, end_r, start_l, end_l) %>% 
      separate(key, c("key", "pos"), sep = "_") %>% 
      spread(key, value) %>%
      select(chrom, start, end, everything(), -pos)
  } 
  
  # not `both`
  if (!strand) {
    if (left) {
      if (fraction) {
        res <- x %>%
          mutate(.start = start,
                 .interval_size = end - start,
                 start = start - round(left * .interval_size),
                 end = .start) %>%
          select(-.start, -.interval_size)
        
      } else {
        res <- x %>%
          mutate(.start = start,
                 start = start - left,
                 end = .start) %>%
          select(-.start)
      }
    
    } else if (right) {
      if (fraction) {
        res <- x %>%
          mutate(.interval_size = end - start,
                 start = end,
                 end = end + round(right * .interval_size)) %>%
          select(-.interval_size)
        
      } else {
        res <- x %>%
          mutate(start = end,
                 end = end + right)
      }
    }
  
  } else {
    # calc left and right based on strand
    if (left) {
      if (fraction) {
        res <- x %>%
          mutate(.interval_size = end - start,
                 start = ifelse(strand == '+',
                                start - round(left * .interval_size),
                                end),
                 end = ifelse(strand == '+',
                              start + round(left * .interval_size),
                              end + round(left * .interval_size))) %>%
          select(-.interval_size)
        
      } else {
        res <- x %>%
          mutate(start = ifelse(strand == '+',
                                start - left,
                                end),
                 end = ifelse(strand == '+',
                              start + left,
                              end + left))
      }
      
      # WHY DOESN'T THIS WORK?
      #res <- x %>%
      #  mutate(.interval_size = end - start,
      #         .shift_amt = ifelse(fraction == T,
      #                             left * .interval_size,
      #                             left),
      #         start = ifelse(strand == '+',
      #                        start - .shift_amt,
      #                        end),
      #         end = ifelse(strand == '+',
      #                      start + .shift_amt,
      #                      end + .shift_amt)) 
     
    } else if (right) {
      if (fraction) {
        res <- x %>%
          mutate(.interval_size = end - start,
                 start = ifelse(strand == '+',
                                end,
                                start - round(right * .interval_size)),
                 end = ifelse(strand == '+',
                              end + round(right * .interval_size),
                              start + round(right * .interval_size))) %>%
          select(-.interval_size)
        
      } else {
        res <- x %>%
          mutate(start = ifelse(strand == '+',
                                end,
                                start - right),
                 end = ifelse(strand == '+',
                              end + right,
                              start + right))
      }
    }
  }    
    
  res <- res %>%
    bound_intervals(genome, trim) %>%
    bed_sort()
  
  res
}
