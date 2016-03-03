#' create flanking intervals from input intervals
#' 
#' @param bed_tbl tbl of intervals
#' @param both number of bases on both sizes 
#' @param left number of bases on left side
#' @param right number of bases on right side
#' @param strand define \code{left} and \code{right} based on strand
#' @param fraction define flanks based on fraction of interval length
#' 
#' @return \code{data.frame}
#' 
#' @examples 
#' genome <- dplyr::tibble(
#'  ~chrom, ~size,
#'  "chr1", 5000
#' )
#' 
#' bed_tbl <- dplyr::tibble(
#'  ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'  "chr1", 500,    1000, '.',   '.',    '+',
#'  "chr1", 1000,   1500, '.',   '.',    '-'
#' )
#' 
#' bed_flank(bed_tbl, left = 100)
#' bed_flank(bed_tbl, right = 100)
#' bed_flank(bed_tbl, both = 100)
#'
#' bed_flank(bed_tbl, both = 0.5, fraction=TRUE)
#' 
#' @export
bed_flank <- function(bed_tbl, both = 0, left = 0,
                      right = 0, fraction = FALSE,
                      strand = FALSE) {

  assert_that(both > 0 || left > 0 || right > 0)
  assert_that(is.flag(strand) && 'strand' %in% colnames(bed_df))
  assert_that(!is.flag(both))
  assert_that(fraction >= 0 && fraction <= 1)
  
  if (both != 0 && (left != 0 || right != 0)) {
    stop('ambiguous side spec for bed_flank')
  } 
  
  if (both) {
    if (fraction) {
      flank_result <- bed_df %>%
        mutate(.interval_size = end - start,
               .starts = list(start - (fraction * .interval_size), end),
               .ends = list(start, end + (fraction * .interval_size))) %>%
        select(-.interval_size)       
    } else {
      flank_result <- bed_df %>%
        mutate(.starts = list(start - both, end),
               .ends = list(start, end + both))
    }
   
    # XXX figure out how to put start, end in original position, they come out the end
    flank_result <- flank_result %>% 
      tidyr::unnest() %>%
      select(-start, -end) %>%
      rename(start = .starts, end = .ends)
    
    flank_result
  } 
  
  # not `both`
  if (!strand) {
    if (left) {
      flank_result <- bed_df %>%
        mutate(.start = start,
               start = start - left,
               end = .start) %>%
        select(-.start) 
    } else if (right) {
      flank_result <- bed_df %>%
        mutate(start = end,
               end = end + right)
    } 
  } else {
    if (left) {
    # calc left and right based on strand
      flank_result <- bed_df %>%
        mutate(start = ifelse(strand == '+',
                              start - left,
                              end),
               end = ifelse(strand == '+',
                            start,
                            end + left))
    } else if (right) {
       flank_result <- bed_df %>%
        mutate(start = ifelse(strand == '+',
                              end,
                              start - right),
               end = ifelse(strand == '+',
                            end + right,
                            start))
    }
  }    
    
  flank_result
}
