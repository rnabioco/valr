#' Subtract intervals.
#' 
#' Subtract \code{y} intervals from \code{x} intervals.
#' 
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param any remove any \code{x} intervals that overlap \code{y}
#' @param strand subtract intervals on same strand
#' @param strand_opp subtract intervals on opposite strand
#' 
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/subtract.html}
#' 
#' @examples
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
#'  # dangling left
#'  "chr1", 900,    1050,
#'  # dangling right
#'  "chr1", 1150,   1250,
#'  # full containment
#'  "chr1", 1299,   1501
#' )
#' 
#' bed_subtract(x, y)
#' bed_subtract(x, y, any = TRUE)
#'  
#' @export
bed_subtract <- function(x, y, any = FALSE,
                         strand = FALSE, strand_opp = FALSE) {

  x <- group_by(x, chrom, add = TRUE)
  y <- bed_merge(y, strand = (strand_opp || strand))
  y <- group_by(y, chrom, add = TRUE)
    

  if (any) {
    # if `any` then only return x intervals without overlaps 
    res <- bed_intersect(x, y, strand, strand_opp)
    # collect x intervals with no overlaps 
    colspec <- c('chrom', 'start' = 'start.x', 'end' = 'end.x')
    anti <- anti_join(x, res, by = colspec)
   
    return(anti)
  }

  # otherwise return the subtracted set - this includes x intervals
  # without overlaps.

  if (!strand && !strand_opp) {
    res <- subtract_impl(x, y)
  } 
  
  if (strand) {
    x_pos <- filter(x, strand == "+") 
    y_pos <- filter(y, strand == "+") 
    x_neg <- filter(x, strand == "-") 
    y_neg <- filter(y, strand == "-") 
    res_pos <- subtract_impl(x_pos, y_pos)
    res_neg <- subtract_impl(x_neg, y_neg)
    res <- bind_rows(res_pos, res_neg)
  } else if (strand_opp) {
    x_pos <- filter(x, strand == "+") 
    y_pos <- filter(y, strand == "+") 
    x_neg <- filter(x, strand == "-") 
    y_neg <- filter(y, strand == "-") 
    res_pos <- subtract_impl(x_pos, y_neg)
    res_neg <- subtract_impl(x_neg, y_pos)
    res <- bind_rows(res_pos, res_neg)
  }
  res <- bed_sort(res)

  res  
}
