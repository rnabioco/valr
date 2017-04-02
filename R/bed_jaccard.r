#' Calculate jaccard statistics on two sets of intervals.
#'
#' @details \code{bed_jaccard()} quantifies the extent of overlap between to sets of
#' intervals. The Jaccard statistic takes values of \code{[0,1]} and is measured as:
#'
#' \deqn{ J(x,y) = \frac{\mid x \bigcap y \mid}
#'                      {\mid x \bigcup y \mid} =
#'                 \frac{\mid x \bigcap y \mid}
#'                      {\mid x \mid + \mid y \mid - \mid x \bigcap y \mid} }
#'
#' @param x \code{\link{tbl_interval}}
#' @param y \code{\link{tbl_interval}}
#'
#' @template stats
#'
#' @family interval statistics
#' @return \code{\link{tbl_interval}} with the following columns:
#'   \itemize{
#'     \item{\code{len_i}}{ length of the intersection}
#'     \item{\code{len_u}}{ length of the union}
#'     \item{\code{jaccard}}{ jaccard statistic}
#'     \item{\code{n_int}}{ number of intersecting intervals between x and y}
#'     }
#'
#' @seealso
#'   \url{http://bedtools.readthedocs.org/en/latest/content/tools/jaccard.html}
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 10,     20,
#'   "chr1", 30,     40
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 15,     20
#' )
#'
#' bed_jaccard(x, y)
#'
#' @export
bed_jaccard <- function(x, y) {

  if (!is.tbl_interval(x)) x <- tbl_interval(x)
  if (!is.tbl_interval(y)) y <- tbl_interval(y)

  x <- bed_merge(x)
  y <- bed_merge(y)

  res_intersect <- bed_intersect(x, y)
  res_intersect <- summarize(res_intersect,
                             sum_overlap = sum(as.numeric(.overlap)),
                             n_int = as.numeric(n()))

  res_x <- mutate(x, .size = end - start)
  res_x <- summarize(res_x, sum_x = sum(as.numeric(.size)))

  res_y <- mutate(y, .size = end - start)
  res_y <- summarize(res_y, sum_y = sum(as.numeric(.size)))

  n_i <- res_intersect$sum_overlap
  n <- res_intersect$n_int

  n_x <- res_x$sum_x
  n_y <- res_y$sum_y

  n_u <- n_x + n_y

  jaccard <- n_i / (n_u - n_i)

  res <- tibble::tribble(
    ~len_i, ~len_u, ~jaccard, ~n,
    n_i,    n_u,    jaccard,  n
  )

  res
}
