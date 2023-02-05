#' Divide intervals into new sub-intervals ("windows").
#'
#' @param x [ivl_df]
#' @param win_size divide intervals into fixed-size windows
#' @param step_size size to step before next window
#' @param num_win divide intervals to fixed number of windows
#' @param reverse reverse window numbers
#'
#' @note The `name` and `.win_id` columns can be used to create new
#'   interval names (see 'namenum' example below) or in subsequent
#'   `group_by` operations (see vignette).
#'
#' @family utilities
#'
#' @return [ivl_df] with `.win_id` column that contains a numeric
#'   identifier for the window.
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 100,    200,  "A",   ".",    "+"
#' )
#'
#' bed_glyph(bed_makewindows(x, num_win = 10), label = ".win_id")
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
#' dplyr::mutate(wins, namenum = stringr::str_c(name, "_", .win_id))
#'
#' @export
bed_makewindows <- function(x,
                            win_size = 0,
                            step_size = 0,
                            num_win = 0,
                            reverse = FALSE) {
  check_required(x)

  x <- check_interval(x)

  if (win_size == 0 && num_win == 0) {
    cli::cli_abort("specify either {.var win_size} or {.var num_win}")
  }

  if (win_size < 0 || num_win < 0) {
    cli::cli_abort("{.var win_size} and {.var num_win} must be >= 0")
  }

  if (any(x$end - x$start < num_win)) {
    cli::cli_alert_warning("interval lengths < {.var num_win} will be skipped.")
  }

  # dummy win_ids
  x <- mutate(x, .win_id = 0)

  res <- makewindows_impl(x, win_size, num_win, step_size, reverse)

  res
}
