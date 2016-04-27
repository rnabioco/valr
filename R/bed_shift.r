#' Adjust intervals by a fixed size.
#'
#' Out-of-bounds intervals are removed by default.
#'  
#' @param bed_tbl tbl of intervals 
#' @param genome chromosome sizes
#' @param size number of bases to shift. postive numbers shift right, negative shift left. 
#' @param strand logical indicating shift according to +/- in strand column
#' @param fraction define \code{size} as a fraction of interval
#' @param trim adjust coordinates for out-of-bounds intervals
#' 
#' @return \code{data.frame}
#' 
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/shift.html}
#' 
#' @examples
#' bed_tbl <- tibble::frame_data(
#'    ~chrom, ~start, ~end, ~strand,
#'    "chr1", 100, 150, "+",
#'    "chr1", 200, 250, "+",
#'    "chr2", 300, 350, "+",
#'    "chr2", 400, 450, "-",
#'    "chr3", 500, 550, "-",
#'    "chr3", 600, 650, "-" 
#'    )
#' 
#' genome <- tibble::frame_data(
#'    ~chrom, ~size,
#'    "chr1", 1000,
#'    "chr2", 2000,
#'    "chr3", 3000
#'    )
#' 
#' bed_shift(bed_tbl, genome, 100)
#' bed_shift(bed_tbl, genome, fraction = 0.5)
#' bed_shift(bed_tbl, genome, 100, strand = TRUE)
#' bed_shift(bed_tbl, genome, strand = TRUE, fraction = 0.5)
#' 
#' @export
bed_shift <- function(bed_tbl, genome, size = 0, strand = FALSE, fraction = 0, trim = FALSE){

  # shift invervals
  if (!strand && !fraction) {
    res <- bed_tbl %>%
      mutate(start = start + size,
             end = end + size)
  }
  
  # shift by percent of interval size
  if (!strand && fraction){
    res <- bed_tbl %>%
      mutate(.size = end - start) %>%
      mutate(start = start + round(.size * fraction),
             end = end + round(.size * fraction)) %>%
      select(-.size)
  }
  
  # shift by strand
  if (strand && !fraction){
    res <- bed_tbl %>%
      mutate(start = ifelse(strand == '+',
                            start + size,
                            start - size),
             end = ifelse(strand == '+',
                          end + size,
                          end - size))
  }
  
  # shift by strand and percent
  if (strand && fraction){
    res <- bed_tbl %>%
      mutate(.size = end - start) %>%
      mutate(start = ifelse(strand == "+",
                            start + round(.size * fraction),
                            start - round(.size * fraction)),
             end = ifelse(strand == "+",
                          end + round(.size * fraction),
                          end - round(.size * fraction))) %>%
      select(-.size)
  }
  
  res <- bound_intervals(res, genome, trim)

  res  
}
