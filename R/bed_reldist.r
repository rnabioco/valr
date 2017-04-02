#' Compute relative distances between intervals.
#'
#' @param x \code{\link{tbl_interval}}
#' @param y \code{\link{tbl_interval}}
#' @param detail report relative distances for each \code{x} interval.
#'
#' @family interval statistics
#'
#' @return If \code{detail = FALSE}, a \code{\link{tbl_interval}} that summarizes
#'  calcuclated \code{.reldist} values with the following columns:
#'   \itemize{
#'     \item{\code{.reldist}}{ relative distance metric}
#'     \item{\code{.counts}}{ number of metric observations}
#'     \item{\code{.total}}{ total observations}
#'     \item{\code{.freq}}{ frequency of observation}}
#'
#'     If \code{detail = TRUE}, a new \code{.reldist} column reports the relative
#'     distance for each input \code{x} interval.
#'
#' @template stats
#'
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/reldist.html}
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom,   ~start,    ~end,
#'   "chr1",    75,       125
#' )
#'
#' y <- tibble::tribble(
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


