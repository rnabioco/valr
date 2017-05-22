#' Compute relative distances between intervals.
#'
#' @param x [tbl_interval()]
#' @param y [tbl_interval()]
#' @param detail report relative distances for each `x` interval.
#'
#' @family interval statistics
#'
#' @return If `detail = FALSE`, a [tbl_interval()] that summarizes
#'  calculated `.reldist` values with the following columns:
#'   - `.reldist` relative distance metric
#'   - `.counts` number of metric observations
#'   - `.total` total observations
#'   - `.freq` frequency of observation
#'
#'   If `detail = TRUE`, a new `.reldist` column reports the relative
#'   distance for each input `x` interval.
#'
#' @template stats
#'
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/reldist.html}
#'
#' @examples
#' x <- trbl_interval(
#'   ~chrom,   ~start,    ~end,
#'   "chr1",    75,       125
#' )
#'
#' y <- trbl_interval(
#'   ~chrom,   ~start,    ~end,
#'   "chr1",    50,       100,
#'   "chr1",    100,       150
#' )
#'
#' bed_reldist(x, y)
#'
#' bed_reldist(x, y, detail = TRUE)
#'
#' @export
bed_reldist <- function(x, y, detail = FALSE) {

  if (!is.tbl_interval(x)) x <- tbl_interval(x)
  if (!is.tbl_interval(y)) y <- tbl_interval(y)

  x <- group_by(x, chrom, add = TRUE)
  y <- group_by(y, chrom, add = TRUE)

  res <- reldist_impl(x, y)

  if (detail) return(res)

  res[['.reldist']] <- floor(res[['.reldist']] * 100) / 100
  nr <- nrow(res)
  res <- group_by(res, .reldist)
  res <- summarize(res,
                   .counts = n(),
                   .total = nr,
                   .freq = .counts / .total)
  res
}


