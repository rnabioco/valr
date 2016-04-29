#' Identify intersecting intervals within a specified distance
#' 
#' @param x BED intervals 
#' @param y BED intervals 
#' @param both add basepairs upstream and downstream to x features
#' @param left add left basepairs to x features 
#' @param right add right basepairs to x features 
#' @param fraction define both, left, and right distance based on fraction of x interval length
#' @param sw define left and right based on strand 
#' @param strand intersect intervals on same strand
#' @param strand_opp intersect intervals on opposite strand
#' @param trim adjust coordinates for out-of-bounds intervals
#' 
#' @examples 
#' x <- tibble::frame_data(
#' ~chrom, ~start, ~end,
#' "chr1", 10,    100,
#' "chr2", 200,    400,
#' "chr2", 300,    500,
#' "chr2", 800,    900
#' )
#' 
#' y <- tibble::frame_data(
#' ~chrom, ~start, ~end,
#' "chr1", 150,    400,
#' "chr2", 230,    430,
#' "chr2", 350,    430
#' )
#' 
#' bed_intersect(x, y)
#' bed_window(x, y, w = 100)
#' 
#'  
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/window.html}
#'  
#' @export
bed_window <- function(x, y, genome, both = 0, left = 0, right = 0,
                       fraction = FALSE, sw = FALSE, strand = FALSE, 
                       strand_opp = FALSE, trim = FALSE){
  
  x <- mutate(x, start.org = start,
                 end.org = end)
  
  x_slop <- bed_slop(x, genome, both = both, left = left,
           right = right, fraction = fraction,
           strand = sw, trim = trim)
  
  res <- bed_intersect(x_slop, y, strand = strand, strand_opp = strand_opp)
  
  # reassign original x bed_df start and end positions
  res <- mutate(res, start.x = start.org.x)
  res <- mutate(res, end.x = end.org.x)
  res <- select(res, -start.org.x, -end.org.x)

  res

}
