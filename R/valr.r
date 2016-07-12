#' valr: genome interval arithmetic in R
#' 
#' valr provides tools to read and manipulate intervals and signals on a 
#' genome reference. While other analysis suites like BEDtools and BEDOPS enable
#' fast command line analysis, valr enables similar analyses \emph{within
#' R}.
#' 
#' valr focuses on manipulating intervals in BED format and signal in 
#' bedGraph format. Eventually this package will power interactive 
#' visualizations of genome-scale data with \code{shiny}.
#' 
#' To learn more about valr, start with the vignette: 
#' \code{browseVignettes(package = "valr")}
#' 
#' @author Jay Hesselberth <jay.hesselberth@@gmail.com>
#'   
#' @docType package
#' @name valr
#'   
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/index.html}
#' @seealso \url{https://pythonhosted.org/pybedtools/}
#' @seealso \url{http://bedops.readthedocs.org/en/latest/index.html}
#' @seealso \url{https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html}
#'   
#' @useDynLib valr
#' @importFrom tibble frame_data
#' @importFrom readr read_tsv col_integer col_character
#' @importFrom stringr str_replace str_split str_c str_length str_detect
#' @importFrom tidyr unnest gather separate spread
#' @importFrom assertthat assert_that is.flag
#' @importFrom purrr by_row
#' @importFrom lazyeval lazy_eval lazy_dots
#' @importFrom stats phyper
#' @importFrom utils head tail
#' @import dplyr
#' @exportPattern "^[[:alpha:]]+"
NULL

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`
