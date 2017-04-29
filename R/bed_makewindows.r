#' Divide intervals into new sub-intervals ("windows").
#'
#' @param x [tbl_interval()]
#' @param genome [tbl_genome()]
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
#' @return [tbl_interval()] with `.win_id` column that contains a numeric
#'   identifier for the window.
#'
#' @examples
#' genome <- trbl_genome(
#'  ~chrom, ~size,
#'  "chr1", 200
#' )
#'
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 100,    200,  'A',   '.',    '+'
#' )
#'
#' bed_glyph(bed_makewindows(x, genome, num_win = 10), label = '.win_id')
#'
#' # Fixed number of windows
#' bed_makewindows(x, genome, num_win = 10)
#'
#' # Fixed window size
#' bed_makewindows(x, genome, win_size = 10)
#'
#' # Fixed window size with overlaps
#' bed_makewindows(x, genome, win_size = 10, step_size = 5)
#'
#' # reverse win_id
#' bed_makewindows(x, genome, win_size = 10, reverse = TRUE)
#'
#' # bedtools 'namenum'
#' wins <- bed_makewindows(x, genome, win_size = 10)
#' dplyr::mutate(wins, namenum = stringr::str_c(name, '_', .win_id))
#'
#' @export
bed_makewindows <- function(x, genome, win_size = 0,
                            step_size = 0, num_win = 0,
                            reverse = FALSE) {

  if (!is.tbl_interval(x)) x <- tbl_interval(x)
  if (!is.tbl_genome(genome)) genome <- tbl_genome(genome)

  if (win_size == 0 && num_win == 0)
    stop('specify either `win_size` or `num_win`', call. = FALSE)

  x <- ungroup(x)
  x <- mutate(x, .row_id = row_number())
  x <- rowwise(x)
  if (num_win > 0) {
    x <- mutate(x,
                 .win_size = round((end - start) / num_win))
  } else {
    x <- mutate(x, .win_size = win_size)
  }

  res <- mutate(x, .start = list(seq(start, end,
                                     by = .win_size - step_size)),
                .win_num = list(seq(1, length(.start))))
  res <- tidyr::unnest(res)
  res <- mutate(res, .end = ifelse(.start + .win_size < end,
                                   .start + .win_size, end))
  res <- filter(res, .start != .end )
  res <- mutate(res, start = .start, end = .end)
  res <- select(res, -.start, -.end, -.win_size)

  # add .win_id column
  res <- group_by(res, .row_id)
  if (reverse) {
    res <- mutate(res, .win_id = rank(-.win_num))
  } else {
    res <- mutate(res, .win_id = rank(.win_num))
  }

  res <- ungroup(res)
  res <- select(res, -.win_num, -.row_id)
  res
}
