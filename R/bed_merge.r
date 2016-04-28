#' Merge overlapping intervals.
#'
#' @param x tbl of intervals 
#' @param max_dist maximum distance between intervals to merge
#' @param .keep keep dot columns (\code{.merge_id}, \code{.overlap})
#' @param ... name-value pairs that specify merging operations
#'
#' @return \code{data_frame}
#'  
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/merge.html}
#' 
#' @examples 
#' x <- tibble::frame_data(
#'  ~chrom, ~start, ~end, ~value, ~strand,
#'  "chr1", 1,      50,   1,      '+',
#'  "chr1", 100,    200,  2,      '+',
#'  "chr1", 150,    250,  3,      '-',
#'  "chr2", 1,      25,   4,      '+',
#'  "chr2", 200,    400,  5,      '-',
#'  "chr2", 400,    500,  6,      '+',
#'  "chr2", 450,    550,  7,      '+'
#' )
#' 
#' bed_merge(x)
#' bed_merge(x, max_dist = 100)
#' bed_merge(x, strand = TRUE)
#' bed_merge(x, .value = sum(value))
#' 
#' @export
bed_merge <- function(x, max_dist = 0, strand = FALSE, ...) {
 
  assert_that(max_dist >= 0)
  
  if ( ! is_sorted(x) ) {
    res <- bed_sort(x)
  }
 
  res <- group_by(res, chrom)
  
  if (strand)
    res <- group_by(res, strand, add = TRUE)

  res <- merge_impl(res, max_dist)
  
  dots <- list(.start = ~min(start), .end = ~max(end))
  dots <- c(dots, lazyeval::lazy_dots(...))
  
  res <- res %>% 
    group_by(chrom, .merge_id) %>%
    summarize_(.dots = dots) %>%
    rename(start = .start, end = .end) %>%
    select(-.merge_id)
  
  attr(res, 'merged') <- TRUE

  res  
}

#' Ask whether a tbl is already merged.
#' 
#' @param x tbl of intervals
#' 
#' @export
is_merged <- function(x) {
  
  merged_attr <- attr(x, "merged")
  
  if ( is.null(merged_attr) || ! merged_attr ) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}

