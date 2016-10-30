#' Identify closest intervals.
#' 
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param overlap include overlapping intervals and report overlap
#' @param suffix colname suffixes in output
#' @param distance_type reporting method for distance to 
#' nearest interval (
#' genome (default) = use negative distances to report upstream intervals,
#'           strand = define upstream based on strand,
#'              abs = report absolute value of distance)
#' 
#' 
#' @return \code{data_frame}
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
#' bed_closest(x, y, overlap = FALSE)
#' 
#' 
#' @export

bed_closest <- function(x, y, overlap = TRUE,
                        strand = FALSE, strand_opp = FALSE, suffix = c('.x', '.y'),
                        distance_type = c("genome", "strand", "abs")) {
  
  if (strand && !('strand' %in% colnames(x) && 'strand' %in% colnames(y)))
    stop("`strand` specified on unstranded data_frame", call. = FALSE)
  
  check_suffix(suffix) 
 
  if ( ! is_sorted(x) )
    x <- bed_sort(x)
  if ( ! is_sorted(y) )
    y <- bed_sort(y)
  
  x <- group_by(x, chrom, add = TRUE)
  y <- group_by(y, chrom, add = TRUE)
 
  suffix <- list(x = suffix[1], y = suffix[2])
  
  if (!strand && !strand_opp){
    res <- closest_impl(x, y, suffix$x, suffix$y)
  }

  strand.x = paste(strand, suffix[1], sep = '')
  strand.y = paste(strand, suffix[2], sep = '')
  
  distance_type <- match.arg(distance_type, c("genome", "strand", "abs"))
  
  if (strand) {
    x_pos <- filter(x, strand == "+") 
    y_pos <- filter(y, strand == "+") 
    x_neg <- filter(x, strand == "-") 
    y_neg <- filter(y, strand == "-") 
    res_pos <- closest_impl(x_pos, y_pos, suffix$x, suffix$y)
    res_neg <- closest_impl(x_neg, y_neg, suffix$x, suffix$y)
    res <- bind_rows(res_pos, res_neg)
    res <- filter(res, strand.x == strand.y) 
  } else if (strand_opp) {
    x_pos <- filter(x, strand == "+") 
    y_pos <- filter(y, strand == "+") 
    x_neg <- filter(x, strand == "-") 
    y_neg <- filter(y, strand == "-") 
    res_pos <- closest_impl(x_pos, y_neg, suffix$x, suffix$y)
    res_neg <- closest_impl(x_neg, y_pos, suffix$x, suffix$y)
    res <- bind_rows(res_pos, res_neg)
    res <- filter(res, strand.x != strand.y) 
  }
  
  # modify distance output based on user input 
  # genome type reporting is default output from closest_impl())
  if (distance_type == "strand" &&  "strand" %in% colnames(x)) {
    res$.distance <- ifelse(res$strand.x == "+", 
                            res$.distance,
                            -(res$.distance))
  } else if (distance_type == "abs") {
    res$.distance <- abs(res$.distance)
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
