#' Divide intervals into new sub-intervals ("windows").
#'
#' @param x [tbl_interval()]
#' @param win_size divide intervals into fixed-size windows
#' @param step_size size to step before next window
#' @param num_win divide intervals to fixed number of windows
#' @param reverse reverse window numbers
#' @param genome this argument has been deprecated and is not used
#'
#' @note The `name` and `.win_id` columns can be used to create new
#'   interval names (see 'namenum' example below) or in subsequent
#'   `group_by` operations (see vignette).
#'
#' @family utilities
#'
#' @return [tbl_interval()] with `.win_id` column that contains a numeric
#'   identifier for the window.
#'
#' @examples
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 100,    200,  'A',   '.',    '+'
#' )
#'
#' bed_glyph(bed_makewindows(x, num_win = 10), label = '.win_id')
#'
#' # Fixed number of windows
#' bed_makewindows(x, num_win = 10)
#'
#' # Fixed window size
#' bed_makewindows(x, win_size = 10)
#'
#' # Fixed window size with overlaps
#' bed_makewindows(x, win_size = 10, step_size = 5)
#'
#' # reverse win_id
#' bed_makewindows(x, win_size = 10, reverse = TRUE)
#'
#' # bedtools 'namenum'
#' wins <- bed_makewindows(x, win_size = 10)
#' dplyr::mutate(wins, namenum = stringr::str_c(name, '_', .win_id))
#'
#' @export
bed_makewindows <- function(x,
                            genome = NULL,
                            win_size = 0,
                            step_size = 0,
                            num_win = 0,
                            reverse = FALSE) {

  # handle deprecated genome argument
  if (!is.null(genome)) {
    if (is.list(genome)) {
      warning("genome argument has been deprecated, ignoring",
        call. = FALSE
      )
    } else if (is.numeric(genome)) {
      # if win_size is passed as positional argument it will be genome
      # reassign to win_size
      win_size <- genome
    }
  }

  if (!is.tbl_interval(x)) x <- as.tbl_interval(x)

  if (win_size == 0 && num_win == 0) {
    stop("specify either `win_size` or `num_win`", call. = FALSE)
  }

  if (win_size < 0 || num_win < 0) {
    stop("`win_size` or `num_win` must be positive values", call. = FALSE)
  }

  # dummy win_ids
  x <- mutate(x, .win_id = 0)
  res <- makewindows_impl(x, win_size, num_win, step_size, reverse)
  res <- as.tbl_interval(res)

  res
}
