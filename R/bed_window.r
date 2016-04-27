#' Identify intersecting intervals within a specified distance
#' 
#' @param x BED intervals 
#' @param y BED intervals 
#' @param w basepairs upstream and downstream to add to x features
#' @param strand intersect intervals on same strand
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
bed_window <- function(x, 
                       y, 
                       w = 1000,
                       l = 0,
                       r = 0,
                       strand = FALSE, 
                       strand_opp = FALSE){
  
  if (strand && x$strand == "-"){
    if (l > 0 || r > 0)
      x$start.org <- x$start
      x$end.org <- x$end
      x$start <- x$start - r  
      x$end <- x$end + l 
  } else {
  if (l > 0 || r > 0) {
    x$start.org <- x$start
    x$end.org <- x$end
    x$start <- x$start - l  
    x$end <- x$end + r 
  } else {
    x$start.org <- x$start
    x$end.org <- x$end
    x$start <- x$start - w  
    x$end <- x$end + w  
  }
  }  
  # need genome info to avoid having intervals past chrom end
  
  res <- bed_intersect(x, y, strand, strand_opp)
  
  res <- select(res, -start.x, -end.x)

  colnames(res)[colnames(res) == "start.org.x"] <- "start.x"
  colnames(res)[colnames(res) == "end.org.x"] <- "end.x"
  
  res
}
