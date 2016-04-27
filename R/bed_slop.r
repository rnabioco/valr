#' Increase the size of input intervals.
#'
#' @inheritParams bed_flank
#' 
#' @return \code{data_frame}
#' 
#' @seealso
#'   \url{http://bedtools.readthedocs.org/en/latest/content/tools/slop.html}
#'   
#' @examples 
#' genome <- tibble::frame_data(
#'  ~chrom, ~size,
#'  "chr1", 5000
#' )
#' 
#' bed_df <- tibble::frame_data(
#'  ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'  "chr1", 500,    1000, '.',   '.',     '+',
#'  "chr1", 1000,   1500, '.',   '.',     '-'
#' )
#' 
#' bed_slop(bed_df, genome, left = 100)
#' bed_slop(bed_df, genome, right = 100)
#' bed_slop(bed_df, genome, both = 100)
#'
#' bed_slop(bed_df, genome, both = 0.5, fraction=TRUE)
#' 
#' @export
bed_slop <- function(bed_df, genome, both = 0, left = 0,
                     right = 0, fraction = FALSE,
                     strand = FALSE) {

  assert_that(is.flag(strand) && 'strand' %in% colnames(bed_df))
  
  if (both != 0 && (left != 0 || right != 0)) {
    stop('ambiguous side spec for bed_slop')
  } 

  if (fraction) {
    bed_df <- mutate(bed_df, .interval_size = end - start)
  }  
  
  if (both != 0) {
    if (fraction) {
      out <- bed_df %>%
        mutate(start = start - both * .interval_size,
               end = end + both * .interval_size)
    } else {
      out <- bed_df %>%
        mutate(start = start - both,
               end = end + both)
    }
  } else {
    # calc left and rigth based on strand
    if (strand) {
      if (fraction) {
        out <- bed_df %>%
          mutate(start = ifelse(strand == '+',
                                start - left * .interval_size,
                                start - right * .interval_size),
                 end = ifelse(strand == '+',
                              end + right * .interval_size,
                              end + left * .interval_size))
      } else {
        out <- bed_df %>%
          mutate(start = ifelse(strand == '+',
                                start - left,
                                end - right),
                 end = ifelse(strand == '+',
                              end + right,
                              end + left))
      }
    } else {
      if ( fraction ) {
        out <- bed_df %>%
          mutate(start = start - left * .interval_size,
                 end = end + right * .interval_size)
      } else {
        out <- bed_df %>%
          mutate(start = start - left,
                 end = end + right)
      }
    }    
  }
   
  if ( fraction ) out <- select(out, -.interval_size) 
  
  out <- out %>%
    bound_intervals(genome) %>%
    bed_sort()

  out

}
