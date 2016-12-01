#' Subtract intervals.
#' 
#' Subtract \code{y} intervals from \code{x} intervals.
#' 
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param any remove any \code{x} intervals that overlap \code{y}
#' 
#' @template groups
#' 
#' @family multi-set-ops 
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/subtract.html}
#' 
#' @examples
#' x <- tibble::tribble(
#' ~chrom, ~start, ~end,
#' 'chr1', 1,      100
#' )
#'  
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 50,     75
#' )
#'  
#' bed_glyph(bed_subtract(x, y))
#' 
#' x <- tibble::tribble(
#'  ~chrom, ~start, ~end,
#'  "chr1", 100,    200,
#'  "chr1", 250,    400,
#'  "chr1", 500,    600,
#'  "chr1", 1000,   1200,
#'  "chr1", 1300,   1500
#' )
#' 
#' y <- tibble::tribble(
#'  ~chrom, ~start, ~end,
#'  "chr1", 150,    175,
#'  "chr1", 510,    525,
#'  "chr1", 550,    575,
#'  "chr1", 900,    1050,
#'  "chr1", 1150,   1250,
#'  "chr1", 1299,   1501
#' )
#' 
#' bed_subtract(x, y)
#' 
#' bed_subtract(x, y, any = TRUE)
#'  
#' @export
bed_subtract <- function(x, y, any = FALSE) {
  
  x <- group_by(x, chrom, add = TRUE)
  y <- bed_merge(y)
  y <- group_by(y, chrom, add = TRUE)
    
  if (any) {
    # collect and return x intervals without overlaps 
    res <- bed_intersect(x, y)
    colspec <- c('chrom', 'start' = 'start.x', 'end' = 'end.x')
    anti <- anti_join(x, res, by = colspec)
   
    return(anti)
  }

  res <- subtract_impl(x, y)
  res <- arrange(res, chrom, start)

  res  
}
