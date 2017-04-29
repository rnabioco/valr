#' Generate randomly placed intervals on a genome.
#'
#' @param genome [tbl_genome()]
#' @param length length of intervals
#' @param n number of intervals to generate
#' @param sort_by sorting variables
#' @param seed seed RNG for reproducible intervals
#'
#' @details Sorting can be suppressed with `sort_by = NULL`.
#'
#' @return [tbl_interval()]
#'
#' @family randomizing operations
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/random.html}
#'
#' @examples
#' genome <- trbl_genome(
#'   ~chrom,  ~size,
#'   "chr1",  10000000,
#'   "chr2",  50000000,
#'   "chr3",  60000000,
#'   "chrX",  5000000
#' )
#'
#' bed_random(genome, seed = 10104)
#'
#' # sorting can be suppressed
#' bed_random(genome, sort_by = NULL, seed = 10104)
#'
#' # 500 random intervals of length 500
#' bed_random(genome, length = 500, n = 500, seed = 10104)
#'
#' @export
bed_random <- function(genome, length = 1000, n = 1e6, sort_by = c('chrom', 'start'), seed = 0) {

  if (!is.tbl_genome(genome)) genome <- tbl_genome(genome)

  if(!all(genome$size > length))
    stop('`length` must be greater than all chrom sizes', call. = FALSE)

  out <- random_impl(genome, length, n, seed)

  if (!is.null(sort_by) && length(sort_by) > 0)
    out <- arrange_(out, .dots = sort_by)

  out <- tibble::as_tibble(out)
  out
}
