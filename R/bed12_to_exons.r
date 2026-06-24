#' Convert BED12 to individual exons in BED6.
#'
#' After conversion to BED6 format, the `score` column contains the exon
#' number, with respect to strand (i.e., the first exon for `-` strand
#' genes will have larger start and end coordinates).
#'
#' @param x [ivl_df], or a path or URL to a BED12 bigBed (`.bb`) file, in which
#'   case the whole file is read and validated as BED12 (via its header's
#'   declared field count) before conversion.
#'
#' @family utilities
#'
#' @examples
#' x <- read_bed12(valr_example("mm9.refGene.bed.gz"))
#'
#' bed12_to_exons(x)
#'
#' # BED12 from any source with the standard column order is accepted,
#' # including `cpp11bigwig::read_bigbed()`, which uses UCSC column names.
#' bb <- system.file("extdata", "test.bb", package = "cpp11bigwig")
#' bed12_to_exons(cpp11bigwig::read_bigbed(bb))
#'
#' # a .bb path can also be passed directly; it is read and validated as BED12
#' bed12_to_exons(bb)
#'
#' @export
bed12_to_exons <- function(x) {
  check_required(x)

  # accept a path/URL to a BED12 bigBed file directly, validating its schema
  if (is_bigbed_path(x)) {
    x <- read_bigbed12(x)
  }

  x <- check_interval(x)

  if (!ncol(x) == 12) {
    cli::cli_abort("expected 12 column input")
  }

  # The parser looks up `exon_sizes`/`exon_starts` by name. If they are absent,
  # assume a foreign BED12 source in standard column order (e.g.
  # cpp11bigwig::read_bigbed(), which uses UCSC names like blockSizes/chromStarts)
  # and standardize names positionally. Input that already carries valr's names is
  # left untouched, so reordered-but-named columns still resolve correctly.
  if (!all(c("exon_sizes", "exon_starts") %in% names(x))) {
    names(x) <- names(bed12_coltypes)
  }

  # The parser reads `exon_sizes`/`exon_starts` as comma-separated character
  # values. Validate up front so non-BED12 input fails with a clear message
  # instead of a cryptic error from the C++ backend.
  not_chr <- c("exon_sizes", "exon_starts")[
    !vapply(x[c("exon_sizes", "exon_starts")], is.character, logical(1))
  ]
  if (length(not_chr) > 0) {
    cli::cli_abort(c(
      "{.arg x} does not look like BED12 input.",
      "x" = "{.field {not_chr}} must be {.cls character} ({.val comma,separated,values})."
    ))
  }

  res <- bed12toexons_impl(x)
  res <- tibble::as_tibble(res)
  res <- res[, c("chrom", "start", "end", "name", "score", "strand")]
  res <- bed_sort(res)

  res
}
