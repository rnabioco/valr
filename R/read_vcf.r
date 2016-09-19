#' Read a VCF file.
#' 
#' @param vcf vcf filename
#'   
#' @return \code{data_frame}
#'   
#' @note return value has \code{chrom}, \code{start} and \code{end} columns.
#'   Interval lengths are the size of the 'REF' field.
#'   
#' @examples
#' v <- system.file('extdata', 'test.vcf.gz', package = 'valr')
#' read_vcf(v)
#' 
#' @export
read_vcf <- function(vcf) {
  
   res <- suppressMessages(readr::read_tsv(vcf, comment = '##'))
   colnames(res) <- stringr::str_replace(colnames(res), '^#', '')
   
   res <- mutate(res,
                        chrom = str_c('chr', CHROM),
                        start = POS,
                        end = start + stringr::str_length(REF))
   
   res
}
