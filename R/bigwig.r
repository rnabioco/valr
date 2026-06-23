# Helpers for using bigWig/bigBed files directly as `y` in interval operations.
#
# When a file path or URL is supplied where an interval dataframe is expected,
# the file is queried with only the regions spanned by `x` (via cpp11bigwig's
# multi-range query), avoiding the cost of materializing the whole file. Works
# for local paths and remote (http/https) URLs.

# Lowercased file extension, ignoring any URL query string or fragment.
#' @noRd
bigwig_ext <- function(file) {
  path <- sub("[?#].*$", "", file)
  tolower(sub(".*\\.([^.]+)$", "\\1", path))
}

# Is `y` a single path/URL to a bigWig or bigBed file?
#' @noRd
is_bigwig_path <- function(y) {
  is.character(y) &&
    length(y) == 1L &&
    !is.na(y) &&
    bigwig_ext(y) %in% c("bw", "bigwig", "bb", "bigbed")
}

# Read a bigWig/bigBed file scoped to the intervals in `x`.
#
# Returns a de-duplicated [ivl_df] of file entries overlapping `x`. Overlapping
# query intervals can otherwise return the same entry more than once.
#' @noRd
read_bigwig_regions <- function(x, file) {
  chrom <- as.character(x[["chrom"]])
  start <- as.integer(x[["start"]])
  end <- as.integer(x[["end"]])

  y <- switch(
    bigwig_ext(file),
    bw = ,
    bigwig = read_bigwig(file, chrom = chrom, start = start, end = end),
    bb = ,
    bigbed = read_bigbed(file, chrom = chrom, start = start, end = end),
    cli::cli_abort(c(
      "Can't read {.file {file}} as a bigWig or bigBed file.",
      "i" = "Supported extensions are {.val bw}, {.val bigwig}, {.val bb}, and {.val bigbed}."
    ))
  )

  dplyr::distinct(y)
}
