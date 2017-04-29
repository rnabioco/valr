#' valr: genome interval arithmetic in R
#'
#' valr provides tools to read and manipulate intervals and signals on a genome
#' reference. valr was developed to facilitate interactive analysis of
#' genome-scale data sets, leveraging the power of dplyr and piping.
#'
#' To learn more about valr, start with the vignette:
#' `browseVignettes(package = "valr")`
#'
#' @author Jay Hesselberth <jay.hesselberth@@gmail.com>
#' @author Kent Riemondy <kent.riemondy@@gmail.com>
#'
#' @docType package
#' @name valr
#'
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/index.html}
#' @seealso \url{https://pythonhosted.org/pybedtools/}
#' @seealso \url{http://bedops.readthedocs.org/en/latest/index.html}
#' @seealso
#'   \url{https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html}
#'
#' @seealso \url{https://CRAN.R-project.org/package=bedr}
#'
#' @useDynLib valr, .registration = TRUE
#' @importFrom tibble tribble as_tibble
#' @importFrom readr read_tsv col_integer col_character col_double
#' @importFrom stringr str_replace str_split str_c str_length
#' @importFrom tidyr unnest
#' @importFrom lazyeval lazy_dots
#' @importFrom stats fisher.test
#' @importFrom broom tidy
#' @importFrom utils head tail
#' @import ggplot2
#' @import dplyr
"_PACKAGE"
