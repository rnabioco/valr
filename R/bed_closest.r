#' Identify closest intervals.
#'
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param overlap report overlapping intervals
#' @param suffix colname suffixes in output
#'
#' @template groups
#'
#' @return \code{data_frame} with additional columns:
#'   \itemize{
#'     \item{\code{.dist}}{ distance to closest interval, negative distances denote upstream intervals}
#'     \item{\code{.overlap}}{ overlap with closest interval}
#'   }
#'
#' @family multi-set-ops
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/closest.html}
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1',      100,     125
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1',      25,      50,
#'   'chr1',     140,     175
#' )
#'
#' bed_glyph(bed_closest(x, y))
#'
#' x <- tibble::tribble(
#' ~chrom, ~start, ~end,
#' "chr1", 500,    600,
#' "chr2", 5000,   6000
#' )
#'
#' y <- tibble::tribble(
#' ~chrom, ~start, ~end,
#' "chr1", 100,    200,
#' "chr1", 150,    200,
#' "chr1", 550,    580,
#' "chr2", 7000,   8500
#' )
#'
#' bed_closest(x, y)
#'
#' bed_closest(x, y, overlap = FALSE)
#'
#' # Report distance based on strand
#' x  <- tibble::tribble(
#' ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
#' "chr1",	10,	20,	"a",	1,	"-"
#' )
#'
#' y <- tibble::tribble(
#'  ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
#'  "chr1",	8,	9,	"b",	1,	"+",
#'  "chr1",	21,	22,	"b",	1, "-"
#' )
#'
#' res <- bed_closest(x, y)
#' res$.strand_dist <- ifelse(res$strand.x == "+", res$.dist, -(res$.dist))
#'
#' #Report absolute distances
#' res$.abs_dist <- abs(res$.dist)
#'
#' @export
bed_closest <- function(x, y, overlap = TRUE,
                        suffix = c('.x', '.y')){

  check_suffix(suffix)

  x <- arrange(x, chrom, start)
  x <- group_by(x, chrom, add = TRUE)

  y <- arrange(y, chrom, start)
  y <- group_by(y, chrom, add = TRUE)

  suffix <- list(x = suffix[1], y = suffix[2])

  res <- closest_impl(x, y, suffix$x, suffix$y)

  # remove negative overlap information
  # not necessary to keep, due to distance column
  res$.overlap <- ifelse(res$.overlap < 0, 0, res$.overlap )

  if (!overlap){
    res <- filter(res, .overlap < 1)
    res <- select(res, -.overlap)
  }

  res
}
