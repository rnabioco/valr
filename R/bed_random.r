#' generate randomly placed intervals on a genome
#' 
#' @param genome genome tbl
#' @param length length of intervals
#' @param n number of intervals to generate
#' @param seed seed RNG for reproducible intervals
#' 
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/random.html}
#'
#' @examples
#' genome <- tibble::frame_data(
#'   ~chrom,  ~size,
#'   "chr1",  10000000,
#'   "chr2",  50000000,
#'   "chr3",  60000000,
#'   "chrX",  5000000
#' ) 
#' 
#' # random intervals (unsorted)
#' bed_random(genome)
#'
#' # 500 random intervals of length 500 
#' bed_random(genome, length = 500, n = 500)
#' 
#' # reproducible random intervals
#' bed_random(genome, seed = 2016)
#' 
#' @export
bed_random <- function(genome, length = 1000, n = 1e6, seed = 0) {
  out <- random_impl(genome, length, n, seed) %>% as_data_frame
  out  
}
