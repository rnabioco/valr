#' Merge overlapping intervals
#' 
#' @param max_dist maximum distance between intervals to merge
#' 
#' @note TODO: implement strand (group_by?), max_dist
#' 
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/merge.html}
#' 
#' @examples 
#' bed_tbl <- dplyr::tibble(
#'  ~chrom, ~start, ~end,
#'  "chr1", 1,      50,
#'  "chr1", 100,    200,
#'  "chr1", 150,    250,
#'  "chr2", 1,      25,
#'  "chr2", 200,    400,
#'  "chr2", 400,    500,
#'  "chr2", 450,    550
#' )
#' 
#' bed_merge(bed_tbl)
#' 
#' bed_merge(bed_tbl, max_dist)
#' 
#' @export
bed_merge <- function(bed_tbl, max_dist = 0) {
 
  assert_that(max_dist >= 0)
  
  if ( ! is_sorted(bed_tbl) ) {
    res <- bed_sort(bed_tbl)
  }
  
  # res must be sorted *again*
  res <- merge_impl(res, max_dist) %>% bed_sort
  
  # add `merged` attribute. this attribute can be tested to determine whether a 
  # merge needs to be done
  attr(res, 'merged') <- TRUE

  res  
}

#' determine whether tbl has been previously merged
#' 
#' @export
is_merged <- function(bed_tbl) {
  
  merged_attr <- attr(bed_tbl, "merged")
  
  if ( is.null(merged_attr) || ! merged_attr ) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}

