#' Increase the size of input intervals.
#'
#' @inheritParams bed_flank
#' @param trim adjust coordinates for res-of-bounds intervals
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
#' x <- tibble::frame_data(
#'  ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'  "chr1", 500,    1000, '.',   '.',     '+',
#'  "chr1", 1000,   1500, '.',   '.',     '-'
#' )
#' 
#' bed_slop(x, genome, left = 100)
#' bed_slop(x, genome, right = 100)
#' bed_slop(x, genome, both = 100)
#'
#' bed_slop(x, genome, both = 0.5, fraction=TRUE)
#' 
#' @export
bed_slop <- function(x, genome, both = 0, left = 0,
                     right = 0, fraction = FALSE,
                     strand = FALSE, trim = FALSE) {

  assert_that(is.flag(strand))
  
  if (strand){
    assert_that('strand' %in% colnames(x))
  }

  if (both != 0 && (left != 0 || right != 0)) {
    stop('ambiguous side spec for bed_slop')
  } 

  if (fraction) {
    x <- mutate(x, .interval_size = end - start)
  }  
  
  if (both != 0) {
    if (fraction) {
      res <- mutate(x, start = start - round( both * .interval_size ),
               end = end + round( both * .interval_size ))
    } else {
      res <- mutate(x, start = start - both,
               end = end + both)
    }
  } else {
    # calc left and right based on strand
    if (strand) {
      if (fraction) {
        res <- mutate(x, start = ifelse(strand == '+',
                                start - round( left * .interval_size ),
                                start - round( right * .interval_size )),
                         end = ifelse(strand == '+',
                              end + round( right * .interval_size ),
                              end + round( left * .interval_size )))
      } else {
        res <- mutate(x, start = ifelse(strand == '+',
                                start - left,
                                start - right),
                         end = ifelse(strand == '+',
                              end + right,
                              end + left))
      }
    } else {
      if ( fraction ) {
        res <- mutate(x, start = start - round( left * .interval_size ),
                 end = end + round( right * .interval_size ))
      } else {
        res <- mutate(x, start = start - left,
                 end = end + right)
      }
    }    
  }
   
  if ( fraction ) res <- select(res, -.interval_size) 
  
  res <- bound_intervals(res, genome, trim)
  res <- bed_sort(res)

  res
}
