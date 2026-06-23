#' Identify intervals within a specified distance.
#'
#' @param x [ivl_df]
#' @param y [ivl_df], or a path or URL to a bigWig (`.bw`) or bigBed (`.bb`)
#'   file. When a file is supplied, only the windowed regions around `x` are
#'   read from it (local files and `http(s)://` URLs are both supported),
#'   avoiding the cost of loading the entire file.
#' @param ... params for bed_slop and bed_intersect
#' @inheritParams bed_slop
#' @inheritParams bed_intersect
#'
#' @template groups
#'
#' @family multiple set operations
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 25,     50,
#'   "chr1", 100,    125
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 60,     75
#' )
#'
#' genome <- tibble::tribble(
#'   ~chrom, ~size,
#'   "chr1", 125
#' )
#'
#' bed_glyph(bed_window(x, y, genome, both = 15))
#'
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 10, 100,
#'   "chr2", 200, 400,
#'   "chr2", 300, 500,
#'   "chr2", 800, 900
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 150,    400,
#'   "chr2", 230,    430,
#'   "chr2", 350,    430
#' )
#'
#' genome <- tibble::tribble(
#'   ~chrom, ~size,
#'   "chr1", 500,
#'   "chr2", 1000
#' )
#'
#' bed_window(x, y, genome, both = 100)
#'
#' # `y` can be a bigWig/bigBed file path or `http(s)://` URL; only the
#' # windowed regions around `x` are read from the file
#' xf <- tibble::tribble(
#'   ~chrom, ~start,   ~end,
#'   "chr1", 4840000, 4841000
#' )
#'
#' gf <- read_genome(valr_example("hg19.chrom.sizes.gz"))
#'
#' bed_window(xf, valr_example("test.bb"), gf, both = 10000)
#'
#' # add a `.dist` column to the output
#' \dontrun{
#' bed_window(x, y, genome, both = 200) |>
#'  mutate(
#'    .dist = case_when(
#'      .overlap == 0 ~ abs(pmax(start.x, start.y) - pmin(end.x, end.y)),
#'      .default = 0
#'    )
#'  )
#' }
#'
#' @seealso \url{https://bedtools.readthedocs.io/en/latest/content/tools/window.html}
#'
#' @export
bed_window <- function(x, y, genome, ...) {
  check_required(x)
  check_required(y)
  check_required(genome)

  x <- check_interval(x)
  # a bigWig/bigBed `y` is resolved below, once the window has been applied
  if (!is_bigwig_path(y)) {
    y <- check_interval(y)
  }
  genome <- check_genome(genome)

  x <- mutate(x, .start = .data[["start"]], .end = .data[["end"]])

  # capture command line args
  cmd_args <- list(...)

  # get arguments for bed_slop and bed_intersect
  slop_arg_names <- names(formals(bed_slop))
  intersect_arg_names <- names(formals(bed_intersect))

  # parse supplied args into those for bed_slop or bed_intersect
  slop_args <- cmd_args[names(cmd_args) %in% slop_arg_names]
  intersect_args <- cmd_args[names(cmd_args) %in% intersect_arg_names]

  # pass new list of args to bed_slop
  slop_x <- do.call(
    bed_slop,
    c(list("x" = x, "genome" = genome), slop_args)
  )

  # for a bigWig/bigBed `y`, query only the windowed (slopped) x regions
  if (is_bigwig_path(y)) {
    y <- read_bigwig_regions(slop_x, y)
  }

  # pass new list of args to bed_intersect
  # TODO: revisit when min_overlap default changes to 1L
  if (!"min_overlap" %in% names(intersect_args)) {
    intersect_args$min_overlap <- 0L
  }
  res <- do.call(
    bed_intersect,
    c(list("x" = slop_x, "y" = y), intersect_args)
  )

  res <- mutate(res, start.x = .data[[".start.x"]], end.x = .data[[".end.x"]])

  res <- ungroup(res)

  res <- select(res, -all_of(c(".start.x", ".end.x")))

  res
}
