#' Identify closest intervals.
#' 
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param overlap report overlapping intervals 
#' @param suffix colname suffixes in output
#' @param dist format for distance to nearest interval
#'              
#' @template groups
#' 
#' @details \code{dist} can take one of the following values:
#'   \itemize{
#'     \item{\code{genome}}{ negative distances signify upstream intervals (default)}
#'     \item{\code{strand}}{ upstream defined based on strand}
#'     \item{\code{abs}}{ absolute value of distance}}
#' 
#' @return \code{data_frame} with additional columns:
#'   \itemize{
#'     \item{\code{.dist}}{ distance to closest interval}
#'     \item{\code{.overlap}}{ overlap with closest interval}
#'   }
#'   
#' @family multi-set-ops
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/closest.html}
#' 
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1',      100,     125
#' )
#' 
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1',      25,      50,
#'   'chr1',     140,     175
#' )
#'  
#' bed_glyph(bed_closest(x, y))
#' 
#' x <- tibble::tribble(
#' ~chrom, ~start, ~end,
#' "chr1", 500,    600,
#' "chr2", 5000,   6000
#' ) 
#' 
#' y <- tibble::tribble(
#' ~chrom, ~start, ~end,
#' "chr1", 100,    200,
#' "chr1", 150,    200,
#' "chr1", 550,    580,
#' "chr2", 7000,   8500
#' ) 
#' 
#' bed_closest(x, y)
#' 
#' bed_closest(x, y, overlap = FALSE)
#' 
#' @export
bed_closest <- function(x, y, overlap = TRUE,
                        suffix = c('.x', '.y'),
                        dist = c("genome", "strand", "abs")) {
  
  check_suffix(suffix) 
 
  if ( ! is_sorted(x) )
    x <- bed_sort(x)
  if ( ! is_sorted(y) )
    y <- bed_sort(y)
  
  x <- group_by(x, chrom, add = TRUE)
  y <- group_by(y, chrom, add = TRUE)
 
  suffix <- list(x = suffix[1], y = suffix[2])
  
  res <- closest_impl(x, y, suffix$x, suffix$y)

  dist <- match.arg(dist, c("genome", "strand", "abs"))
  
  # modify distance output based on user input 
  # genome type reporting is default output from closest_impl())
  if (dist == "strand") {
    res$.dist <- ifelse(res$strand.x == "+", res$.dist, -(res$.dist))
  } else if (dist == "abs") {
    res$.dist <- abs(res$.dist)
  } 
  
  # remove negative overlap information 
  # not necessary to keep, due to distance column
  res$.overlap <- ifelse(res$.overlap < 0, 0, res$.overlap )
  
  if (!overlap){
    res <- filter(res, .overlap < 1)
    res <- select(res, -.overlap)
  }
    
  res
}
