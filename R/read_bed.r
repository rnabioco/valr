#' @title Read BED and related files.
#'
#' @description read functions for BED and related formats. Filenames can be
#'   local file or URLs. The read functions load data into tbls with consistent
#'   `chrom`, `start` and `end` colnames.
#'
#' @param filename file or URL
#' @param n_fields number fields in the BED file
#' @param col_types column type spec for [readr::read_tsv()]
#' @param sort sort the tbl by chrom and start
#' @param ... options to pass to [readr::read_tsv()]
#'
#' @return [ivl_df]
#'
#' @family read functions
#'
#' @details \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}
#'
#' @examples
#' # read_bed assumes 3 field BED format.
#' read_bed(valr_example('3fields.bed.gz'))
#'
#' read_bed(valr_example('6fields.bed.gz'), n_fields = 6)
#'
#' # result is sorted by chrom and start unless `sort = FALSE`
#' read_bed(valr_example('3fields.bed.gz'), sort = FALSE)
#'
#' @export
read_bed <- function(filename, n_fields = 3, col_types = bed12_coltypes,
                     sort = TRUE, ...) {

  check_required(filename)

  coltypes <- col_types[1:n_fields]
  colnames <- names(coltypes)

  out <- readr::read_tsv(
    filename,
    col_names = colnames,
    col_types = coltypes, ...
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
#' read_bed12(valr_example('mm9.refGene.bed.gz'))
#'
#' @export
read_bed12 <- function(filename, ...) {
  check_required(filename)
  bed12_tbl <- read_bed(filename, n_fields = 12)
  bed12_tbl
}

#' @rdname read_bed
#'
#' @details \url{https://genome.ucsc.edu/goldenPath/help/bedgraph.html}
#'
#' @examples
#'
#' read_bedgraph(valr_example('test.bg.gz'))
#'
#' @export
read_bedgraph <- function(filename, ...) {
  # load as bed4, rename `value` column and covert to double
  check_required(filename)
  out <- read_bed(filename, n_fields = 4, sort = FALSE)
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
#' read_narrowpeak(valr_example('sample.narrowPeak.gz'))
#'
#' @export
read_narrowpeak <- function(filename, ...) {
  check_required(filename)
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
#' read_broadpeak(valr_example('sample.broadPeak.gz'))
#'
#' @export
read_broadpeak <- function(filename, ...) {
  check_required(filename)
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


#' Import and convert a bigwig file into a valr compatible tbl
#'
#' @description This function will output a 5 column tibble with
#' zero-based chrom, start, end, score, and strand columns.
#'
#' @param path path to bigWig file
#' @param set_strand strand to add to output (defaults to "+")
#'
#' @note This functions uses \code{rtracklayer} to import bigwigs which
#' has unstable support for the windows platform and therefore may error
#' for windows users (particularly for 32 bit window users).
#'
#' @examples
#' \dontrun{
#' if (.Platform$OS.type != "windows") {
#'   bw <- read_bigwig(valr_example('hg19.dnase1.bw'))
#'   head(bw)
#' }
#' }
#' @importFrom rtracklayer import
#' @export
read_bigwig <- function(path, set_strand = "+") {
  check_required(path)
  # note that rtracklayer will produce a one-based GRanges object
  res <- rtracklayer::import(path)
  res <- dplyr::as_tibble(res)
  res <- dplyr::mutate(res,
                       chrom = as.character(seqnames),
                       start = start - 1L,
                       strand = set_strand)
  dplyr::select(res, chrom, start, end, score, strand)
}

#' Import and convert a GTF/GFF file into a valr compatible bed tbl format
#'
#' @description This function will output a tibble with the
#' required chrom, start, and end columns, as well as other columns depending
#' on content in GTF/GFF file.
#'
#' @param path path to gtf or gff file
#' @param zero_based if TRUE, convert to zero based
#'
#' @examples
#'
#' gtf <- read_gtf(valr_example('hg19.gencode.gtf.gz'))
#' head(gtf)
#'
#' @importFrom rtracklayer import
#' @export
read_gtf <- function(path, zero_based = TRUE){
  gtf <- rtracklayer::import(path)
  gtf <- as.data.frame(gtf)
  gtf <- dplyr::mutate_if(gtf, is.factor, as.character)
  res <- dplyr::rename(gtf, chrom = seqnames)

  if(zero_based) {
    res <- dplyr::mutate(res, start = start - 1L)
  }

  tibble::as_tibble(res)
}
