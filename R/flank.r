#' flank
#' 
#' make flanks from input BED regions
#' 
#' @param strand define -l and -r based on strand
#' @param percent define flanks based on fraction of feature's length
#' @param left number of bases on left side
#' @param right number of bases on right side
#' 
#' @return \code{dplyr::tbl_df}
#' 
#' @export
bedtools_flank <- function(bed_df, size = 0, both = FALSE, left = FALSE,
                           right = FALSE, strand = FALSE, fraction = 0) {

  if (strand && !'strand' %in% colnames(bed_df)) {
    stop('missing strand column for stranded flank calculation')
  } 
  if (both && (left || right)) {
    stop('ambiguous side spec for flank')
  } 
  
  if (both) {
    bed_df <- bed_df %>%
      dplyr::mutate(start = start - both,
                    end = end + both)
    return(bed_df)
  } 
  
  if (!strand) {
    if (left) {
      bed_df <- bed_df %>%
        dplyr::mutate(start = start - left)
      
    } else if (right) {
      bed_df <- bed_df %>%
        dplyr::mutate(end = end + right)
    }      
  } else {
    # calc left and rigth based on strand
  }    
    
  return(bed_df)
}