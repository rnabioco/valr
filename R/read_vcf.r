#' Read a VCF file.
#'
#' @param vcf vcf filename
#'
#' @family read functions
#' @return `data_frame`
#'
#' @note return value has `chrom`, `start` and `end` columns.
#'   Interval lengths are the size of the 'REF' field.
#'
#' @examples
#' vcf_file <- valr_example('test.vcf.gz')
#' read_vcf(vcf_file)
#'
#' @export
read_vcf <- function(vcf) {
  res <- suppressMessages(readr::read_tsv(vcf, comment = "##"))
  colnames(res) <- stringr::str_replace(colnames(res), "^#", "")

  res <- mutate(
    res,
    chrom = stringr::str_c("chr", CHROM),
    start = POS,
    end = start + stringr::str_length(REF)
  )

  res
}
