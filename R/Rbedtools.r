#' Rbedtools: genome arithmetic in R
#' 
#' Rbedtools provides tools to read and manipulate intervals and signals on a 
#' genome reference. While other analysis suites like BEDtools and BEDOPS enable
#' fast command line analysis, Rbedtools aims to enable similar analyses 
#' \emph{within R}.
#' 
#' The package focuses on manipulating intervals in BED format 
#' and signal in bedGraph format. Eventually this package will power interactive
#' visualizations of genome-scale data with \code{shiny}.
#' 
#' To learn more about Rbedtools, start with the vignette: 
#' \code{browseVignettes(package = "Rbedtools")}
#' 
#' @author Jay Hesselberth <jay.hesselberth@gmail.com>
#'   
#' @docType package
#' @name Rbedtools
#'   
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/index.html}
#' @seealso \url{http://bedops.readthedocs.org/en/latest/index.html}
#'   
#' @import dplyr
#' @import stringr
#' @import readr
#' @import tidyr
#' @import assertthat
#' @importFrom purrr by_row
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @useDynLib Rbedtools
NULL
