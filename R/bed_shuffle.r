#' Shuffle input intervals.
#'
#' @param x [tbl_interval()]
#' @param genome [tbl_genome()]
#' @param incl [tbl_interval()] of included intervals
#' @param excl [tbl_interval()] of excluded intervals
#' @param max_tries maximum tries to identify a bounded interval
#' @param within shuffle within chromosomes
#' @param seed seed for reproducible intervals
#'
#' @return [tbl_interval()]
#' @family randomizing operations
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/shuffle.html}
#'
#' @examples
#' genome <- trbl_genome(
#'  ~chrom, ~size,
#'  "chr1", 1e6,
#'  "chr2", 2e6,
#'  "chr3", 4e6
#' )
#'
#' x <- bed_random(genome)
#' bed_shuffle(x, genome)
#'
#' @export
bed_shuffle <- function(x, genome, incl = NULL, excl = NULL,
                        max_tries = 1000, within = FALSE, seed = 0) {

  if (!is.tbl_interval(x)) x <- tbl_interval(x)
  if (!is.tbl_genome(genome)) genome <- tbl_genome(genome)

  # flatten incl and excl
  if (!is.null(incl))
      incl <- bed_merge(incl)
  if (!is.null(excl))
      excl <- bed_merge(excl)

  # make genome into an interval tbl
  genome_incl <- mutate(genome, start = 1, end = size)
  genome_incl <- select(genome_incl, chrom, start, end)

  # find the included intervals bounds. case where only incl intervals are
  # defined is not evaluated explicitly, so is the default
  if (is.null(incl) && is.null(excl)) {
    incl <- genome_incl
  } else if (is.null(incl) && !is.null(excl)) {
    incl <- bed_subtract(genome_incl, excl)
  } else if (!is.null(incl) && !is.null(excl)) {
    incl <- bed_subtract(incl, excl)
  }

  if (nrow(incl) == 0 || is.null(incl))
    stop('no intervals to sample from', call. = FALSE)

  #shuffle_impl will drop all columns except chrom, start, and end
  res <- shuffle_impl(x, incl, within, max_tries, seed)
  res <- as_data_frame(res)

  # by default pass all original x column data to
  # result (except chrom, start, end) which are shuffled
  # see issue # 81 in github

  res <- bind_cols(res, x[, !colnames(x) %in% colnames(res)])

  res
}
