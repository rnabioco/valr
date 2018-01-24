#' Sort a set of intervals.
#'
#' @param x [tbl_interval()]
#' @param by_size sort by interval size
#' @param by_chrom sort within chromosome
#' @param reverse reverse sort order
#'
#' @seealso
#' \url{http://bedtools.readthedocs.org/en/latest/content/tools/sort.html}
#'
#' @examples
#' x <- trbl_interval(
#'    ~chrom, ~start, ~end,
#'    "chr8", 500,    1000,
#'    "chr8", 1000,   5000,
#'    "chr8", 100,    200,
#'    "chr1", 100,    300,
#'    "chr1", 100,    200
#' )
#'
#' # sort by chrom and start
#' bed_sort(x)
#'
#' # reverse sort order
#' bed_sort(x, reverse = TRUE)
#'
#' # sort by interval size
#' bed_sort(x, by_size = TRUE)
#'
#' # sort by decreasing interval size
#' bed_sort(x, by_size = TRUE, reverse = TRUE)
#'
#' # sort by interval size within chrom
#' bed_sort(x, by_size = TRUE, by_chrom = TRUE)
#'
#' @export
bed_sort <- function(x, by_size = FALSE, by_chrom = FALSE, reverse = FALSE) {
  if (!is.tbl_interval(x)) x <- tbl_interval(x)

  if (by_size) {
    res <- mutate(x, .size = end - start)

    if (by_chrom) {
      if (reverse) {
        res <- arrange(res, chrom, desc(.size))
      } else {
        res <- arrange(res, chrom, .size)
      }
    } else {
      if (reverse) {
        res <- arrange(res, desc(.size))
      } else {
        res <- arrange(res, .size)
      }
    }



    # remove .size column and groups in result
    res <- select(res, -.size)
  } else {
    if (by_chrom) {
      res <- group_by(x, chrom)
    }

    # sort by coordinate
    if (reverse) {
      res <- arrange(x, chrom, desc(start))
    } else {
      res <- arrange(x, chrom, start, end)
    }
  }

  # add `sorted` attribute
  attr(res, "sorted") <- TRUE

  res <- as.tbl_interval(res)

  res
}
