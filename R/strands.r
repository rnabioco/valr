#' Flip strands in intervals.
#' 
#' Flips \code{+} stranded intervals to \code{-} and vice-versa.
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

  x <- mutate(x, .strand = ifelse(strand == '+', '-', '+'))
  x <- select(x, -strand)
  x <- rename(x, strand = .strand)
  
  x
}
