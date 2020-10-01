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
#' @seealso Report bugs at \url{https://github.com/rnabioco/valr/issues}
#'
#' @useDynLib valr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom tibble tribble as_tibble is_tibble
#' @importFrom readr read_tsv col_integer col_character col_double
#' @importFrom stringr str_replace str_split str_c str_length fixed
#' @importFrom rlang quos sym syms
#' @importFrom stats fisher.test na.omit
#' @importFrom utils head tail packageVersion
#' @importFrom broom tidy
#' @import ggplot2
#' @import dplyr
"_PACKAGE"
