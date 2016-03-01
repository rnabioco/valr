#' increase the size of input intervals
#'
#' @inheritParams bed_flank
#' 
#' @return \code{data.frame}
#' 
#' @examples 
#' genome <- dplyr::tibble(
#'  ~chrom, ~size,
#'  "chr1", 5000
#' )
#' 
#' bed_df <- dplyr::tibble(
#'  ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'  "chr1", 500,    1000, '.',   '.',     '+',
#'  "chr1", 1000,   1500, '.',   '.',     '-'
#' )
#' 
#' bed_df %>% bed_slop(left = 100, genome = genome)
#' bed_df %>% bed_slop(right = 100, genome)
#' bed_df %>% bed_slop(both = 100, genome)
#'
#' bed_df %>% bed_slop(both = 0.5, fraction=TRUE)
#' 
#' @export
bed_slop <- function(bed_df, genome, both = 0, left = 0,
                     right = 0, fraction = TRUE,
                     strand = FALSE) {

  assert_that(is.flag(strand) && 'strand' %in% colnames(bed_df))
  assert_that(!is.flag(both))
  
  if (both != 0 && (left != 0 || right != 0)) {
    stop('ambiguous side spec for bed_slop')
  } 
  
  if (!genome) {
    warning('genome file not specified. computed intervals may be out-of-bounds.')
  }
 
  if (both) {
    if (fraction) {
      slop_result <- bed_df %>%
        mutate(.interval_size = end - start,
               start = start ,
               end = end + both)
    } else {
      slop_result <- bed_df %>%
        mutate(start = start - both,
               end = end + both)
    }
    return(slop_result)
  } 
  
  if (!strand) {
    if (left) {
      slop_result <- bed_df %>%
        mutate(start = start - left)
    } else if (right) {
      slop_result <- bed_df %>%
        mutate(end = end + right)
    } 
  } else {
    # calc left and rigth based on strand
    slop_result <- bed_df %>%
      mutate(start = ifelse(strand == '+', start, end))
  }    
    
  slop_result
}