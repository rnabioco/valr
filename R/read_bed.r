#' sniff number of fields in a bed file
#' @noRd
sniff_fields <- function(filename) {
  ncol(
    readr::read_tsv(
      filename,
      n_max = 10,
      comment = "#",
      show_col_types = FALSE,
      name_repair = "minimal"
    )
  )
}

#' @title Read BED and related files.
#'
#' @description read functions for BED and related formats. Filenames can be
#'   local file or URLs. The read functions load data into tbls with consistent
#'   `chrom`, `start` and `end` colnames.
#'
#' @param filename file or URL
#' @param col_types column type spec for [readr::read_tsv()]
#' @param sort sort the tbl by chrom and start
#' @param ... options to pass to [readr::read_tsv()]
#' @param n_fields `r lifecycle::badge("deprecated")`
#'
#' @return [ivl_df]
#'
#' @family read functions
#'
#' @details \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}
#'
#' @examples
#' # read_bed assumes 3 field BED format.
#' read_bed(valr_example("3fields.bed.gz"))
#'
#' # result is sorted by chrom and start unless `sort = FALSE`
#' read_bed(valr_example("3fields.bed.gz"), sort = FALSE)
#'
#' @export
read_bed <- function(
  filename,
  col_types = bed12_coltypes,
  sort = TRUE,
  ...,
  n_fields = NULL
) {
  check_required(filename)

  if (!is.null(n_fields)) {
    lifecycle::deprecate_warn(
      "0.6.9",
      "read_bed(n_fields)",
      details = "fields are now determined automatically from the file"
    )
  }

  n_fields <- sniff_fields(filename)

  coltypes <- col_types[1:n_fields]
  colnames <- names(coltypes)

  out <- readr::read_tsv(
    filename,
    col_names = colnames,
    col_types = coltypes,
    ...
  )

  if (sort) out <- bed_sort(out)

  out
}

#' @rdname read_bed
#'
#' @details \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}
#'
#' @examples
#'
#' read_bed12(valr_example("mm9.refGene.bed.gz"))
#'
#' @export
read_bed12 <- function(filename, ...) {
  check_required(filename)
  n_fields <- sniff_fields(filename)
  if (n_fields != 12) {
    cli::cli_abort("expected 12 columns in bed12")
  }
  bed12_tbl <- read_bed(filename)
  bed12_tbl
}

#' @rdname read_bed
#'
#' @details \url{https://genome.ucsc.edu/goldenPath/help/bedgraph.html}
#'
#' @examples
#'
#' read_bedgraph(valr_example("test.bg.gz"))
#'
#' @export
read_bedgraph <- function(filename, ...) {
  # load as bed4, rename `value` column and covert to double
  check_required(filename)
  n_fields <- sniff_fields(filename)
  if (n_fields != 4) {
    cli::cli_abort("expected 4 columns in bedgraph")
  }
  out <- read_bed(filename, sort = FALSE)
  out <- select(out, everything(), value = name)
  out <- mutate(out, value = as.double(value))
  out
}

#' @rdname read_bed
#'
#' @details \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format12}
#'
#' @examples
#'
#' read_narrowpeak(valr_example("sample.narrowPeak.gz"))
#'
#' @export
read_narrowpeak <- function(filename, ...) {
  check_required(filename)
  n_fields <- sniff_fields(filename)
  colnames <- names(peak_coltypes)
  out <- readr::read_tsv(
    filename,
    col_types = peak_coltypes,
    col_names = colnames
  )
  out
}

#' @rdname read_bed
#'
#' @details \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format13}
#'
#' @examples
#'
#' read_broadpeak(valr_example("sample.broadPeak.gz"))
#'
#' @export
read_broadpeak <- function(filename, ...) {
  check_required(filename)
  n_fields <- sniff_fields(filename)
  coltypes <- peak_coltypes[1:length(peak_coltypes) - 1]
  colnames <- names(coltypes)
  out <- readr::read_tsv(filename, col_names = colnames, col_types = coltypes)
  out
}

peak_coltypes <- list(
  chrom = readr::col_character(),
  start = readr::col_integer(),
  end = readr::col_integer(),
  name = readr::col_character(),
  score = readr::col_integer(),
  strand = readr::col_character(),
  signal = readr::col_double(),
  pvalue = readr::col_double(),
  qvalue = readr::col_double(),
  peak = readr::col_integer()
)

bed12_coltypes <- list(
  chrom = readr::col_character(),
  start = readr::col_integer(),
  end = readr::col_integer(),
  name = readr::col_character(),
  score = readr::col_character(),
  strand = readr::col_character(),
  cds_start = readr::col_integer(),
  cds_end = readr::col_integer(),
  item_rgb = readr::col_character(),
  exon_count = readr::col_integer(),
  exon_sizes = readr::col_character(),
  exon_starts = readr::col_character()
)


#' Read a bigwig file into a valr compatible tbl
#'
#' This function will output a 4 column tibble with
#' zero-based chrom, start, end, value columns.
#'
#' @param path path to bigWig file
#' @param ... params for `cpp11bigwig::read_bigwig()`
#'
#' @examples
#' read_bigwig(valr_example("hg19.dnase1.bw"))
#'
#' read_bigwig(valr_example("hg19.dnase1.bw"), as = "GRanges")
#'
#' @export
read_bigwig <- function(path, ...) {
  cpp11bigwig::read_bigwig(path, ...)
}

#' Import and convert a GTF/GFF file into a valr compatible bed tbl format
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function will output a tibble with the
#' required chrom, start, and end columns, as well as other columns depending
#' on content in GTF/GFF file.
#'
#' @param path path to gtf or gff file
#' @param zero_based if TRUE, convert to zero based
#'
#' @examples
#'
#' \dontrun{
#' gtf <- read_gtf(valr_example("hg19.gencode.gtf.gz"))
#' head(gtf)
#' }
#'
#' @export
read_gtf <- function(path, zero_based = TRUE) {
  lifecycle::deprecate_stop(
    when = "0.8.3",
    what = "read_gtf()",
    details = c(
      x = paste0(
        "read_gtf() was removed because rtracklayer does not pass ",
        "CRAN AddressSantizer checks of the UCSC C-library code vendored ",
        "in rtracklayer."
      ),
      i = paste0(
        "convert GTF to BED, and then `read_bed()`, ",
        "or use `rtracklayer::import()` then `gr_to_bed()`."
      )
    )
  )
}
