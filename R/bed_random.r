#' generate randomly placed intervals
#' 
#' @param genome genome tbl
#' @param length length of intervals
#' @param n number of intervals to generate
#' @param seed seed RNG for reproducible intervals
#' 
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/random.html}
#'
#' @example
#' genome <- dplyr::tibble(
#'   ~chrom,  ~size,
#'   "chr1",  10000000,
#'   "chr2",  50000000,
#'   "chr3",  60000000,
#'   "chrX",  5000000
#' ) 
#' 
#' bed_random(genome)
#' 
#' @export
bed_random <- function(genome, length = 100, n = 100, seed = NULL) {

  # randomly select chrom, sizes from genome
  set.seed(seed)
  rows <- genome[sample(nrow(genome), n, replace = TRUE), ]
  
  set.seed(seed)
  res <- rows %>%
    rowwise() %>%
    mutate(start = as.double(sample(size, 1))) %>%
    mutate(end = start + length) %>%
    select(-size) %>% 
    ungroup()
    
  res 
}
