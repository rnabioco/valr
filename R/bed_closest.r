#' Identify closest intervals.
#' 
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param overlap include overlapping intervals and report overlap
#' @param strand intersect intervals on same strand
#' @param strand_opp intersect intervals on opposite strands
#' @param suffix colname suffixes in output
#' @param distance_type reporting method for distance to 
#' nearest interval (
#' genome (default) = use negative distances to report upstream intervals,
#'           strand = define upstream based on strand,
#'              abs = report absolute value of distance)
#' 
#' 
#' @return \code{data_frame}
#' 
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/closest.html}
#' 
#' @examples
#' x <- tibble::frame_data(
#' ~chrom, ~start, ~end,
#' "chr1", 500,    600,
#' "chr2", 5000,   6000
#' ) 
#' 
#' y <- tibble::frame_data(
#' ~chrom, ~start, ~end,
#' "chr1", 100,    200,
#' "chr1", 150,    200,
#' "chr1", 550,    580,
#' "chr2", 7000,   8500
#' ) 
#' 
#' bed_closest(x, y)
#' bed_closest(x, y, overlap = F)
#' 
#' 
#' @export

bed_closest <- function(x, y, overlap = TRUE,
                        strand = FALSE, strand_opp = FALSE, suffix = c('.x', '.y'),
                        distance_type = c("genome", "strand", "abs")) {
  if ( ! is_sorted(x) )
    x <- bed_sort(x)
  if ( ! is_sorted(y) )
    y <- bed_sort(y)
  
  if (is.null(groups(x)) || groups(x) != "chrom")
    x <- group_by(x, chrom)
  if (is.null(groups(y)) || groups(y) != "chrom")
    y <- group_by(y, chrom)
  
  # dplyr::check_suffix
  if (!is.character(suffix) || length(suffix) != 2) {
    stop("`suffix` must be a character vector of length 2.", call. = FALSE)
  }
  
  suffix <- list(x = suffix[1], y = suffix[2])
  
  res <- closest_impl(x, y, suffix$x, suffix$y)

  strand.x = paste(strand, suffix[1], sep = '')
  strand.y = paste(strand, suffix[2], sep = '')
  
  distance_type <- match.arg(distance_type, c("genome", "strand", "abs"))
  
  if (strand) {
    if (!('strand' %in% colnames(x) && 'strand' %in% colnames(y))){
      stop("`strand` specified on unstranded data_frame", call. = FALSE)
    }
    res <- filter(res, strand.x == strand.y) 
  } else if (strand_opp) {
    res <- filter(res, strand.x != strand.y) 
  }
  
  # modify distance output based on user input 
  # genome type reporting is default output from closest_impl())
  if (distance_type == "strand" &&  "strand" %in% colnames(x)) {
    res$.distance <- ifelse(strand.x == "+", 
                            .distance,
                            -(.distance))
  } else if (distance_type == "abs") {
    res$.distance <- abs(res$.distance)
  } 
  
  # remove negative overlap information 
  # not necessary to keep, due to distance column
  res$.overlap <- ifelse(res$.overlap < 0, 0, res$.overlap )
  
  if (!overlap){
    res <- res %>% 
      filter(.overlap < 1) %>% 
      select(-.overlap)
  }
    
  res
}
