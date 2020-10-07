#' Shuffle input intervals.
#'
#' @param x [ivl_df]
#' @param genome [genome_df]
#' @param incl [ivl_df] of included intervals
#' @param excl [ivl_df] of excluded intervals
#' @param max_tries maximum tries to identify a bounded interval
#' @param within shuffle within chromosomes
#' @param seed seed for reproducible intervals
#'
#' @return [ivl_df]
#'
#' @family randomizing operations
#'
#' @seealso \url{https://bedtools.readthedocs.io/en/latest/content/tools/shuffle.html}
#'
#' @examples
#' genome <- tibble::tribble(
#'  ~chrom, ~size,
#'  "chr1", 1e6,
#'  "chr2", 2e6,
#'  "chr3", 4e6
#' )
#'
#' x <- bed_random(genome, seed = 1010486)
#'
#' bed_shuffle(x, genome, seed = 9830491)
#'
#' @export
bed_shuffle <- function(x, genome, incl = NULL, excl = NULL,
                        max_tries = 1000, within = FALSE, seed = 0) {
  x <- check_interval(x)
  genome <- check_genome(genome)

  # flatten incl and excl
  if (!is.null(incl)) {
    incl <- bed_merge(incl)
  }
  if (!is.null(excl)) {
    excl <- bed_merge(excl)
  }

  # make genome into an interval tbl
  genome_incl <- mutate(genome, start = 0, end = size)
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

  if (nrow(incl) == 0 || is.null(incl)) {
    stop("no intervals to sample from", call. = FALSE)
  }

  res <- shuffle_impl(x, incl, within, max_tries, seed)

  # bind original x column data to result (#81)
  res <- bind_cols(res, as_tibble(x[, !colnames(x) %in% colnames(res)]))

  as_tibble(res)
}
