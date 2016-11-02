#' Flip strands in intervals.
#' 
#' Flips positive (\code{+}) stranded intervals to negative (\code{-}) strands,
#' and vice-versa. Facilitates comparisons among intervals on opposing strands.
#' 
#' @param x tbl of intervals
#'   
#' @examples 
#' x <- tibble::tribble(
#' ~chrom, ~start, ~end, ~strand,
#' 'chr1', 1,      100,  '+',
#' 'chr2', 1,      100,  '-' 
#' )
#' 
#' flip_strands(x)
#'  
#' @export
flip_strands <- function(x) {
  
  if (! 'strand' %in% colnames(x))
    stop('`strand` column not found in `x`', call. = FALSE)

  # remove existing groups
  groups_x <- groups(x)
  res <- ungroup(x)
  
  res <- mutate(res, .strand = ifelse(strand == '+', '-', '+'))
  res <- select(res, -strand)
  res <- rename(res, strand = .strand)
 
  res <- group_by_(res, .dots = groups_x) 
  res
}
