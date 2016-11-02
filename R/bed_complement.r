#' Identify intervals in a genome not covered by a query.
#' 
#' @param x tbl of intervals
#' @param genome chrom sizes
#' 
#' @family single-set-ops
#' 
#' @return \code{data_frame}
#' 
#' @examples 
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1',      1,      10,
#'   'chr1',      75,    100    
#' )
#' 
#' genome <- tibble::tribble(
#'   ~chrom, ~size, 
#'   'chr1', 200
#' )
#' 
#' bed_glyph(bed_complement(x, genome))
#' 
#' genome <- tibble::tribble(
#'    ~chrom,  ~size,
#'    "chr1", 500,
#'    "chr2", 600,
#'    "chr3", 800
#' ) 
#' 
#' x <- tibble::tribble(
#'    ~chrom, ~start, ~end,
#'    "chr1", 100,    300,
#'    "chr1", 200,    400,
#'    "chr2", 1,      100,
#'    "chr2", 200,    400,
#'    "chr3", 500,    600
#' )
#' 
#' # intervals not covered by x
#' bed_complement(x, genome)

#' @export
bed_complement <- function(x, genome) {

  if ( ! is_merged(x) ) {
    res <- bed_merge(x)
  } 
 
  # non-overlapping chroms
  chroms_no_overlaps <- anti_join(genome, res, by = 'chrom')
  chroms_no_overlaps <- mutate(chroms_no_overlaps, start = 1)
  chroms_no_overlaps <- select(chroms_no_overlaps, chrom, start, end = size)
 
  # remove rows from x that are not in genome
  res <- semi_join(res, genome, by = 'chrom')
  
  res <- group_by(res, chrom)
  
  res <- complement_impl(res, genome)
  res <- as_data_frame(res)

  res <- bind_rows(res, chroms_no_overlaps)
  res <- bed_sort(res)
 
  res 
}
