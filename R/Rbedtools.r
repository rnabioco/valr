#' Rbedtools: genome arithmetic in R
#' 
#' Rbedtools provides tools implemented in the BEDtools suite to 
#' read and manipulate intervals and signals on a genome reference.
#' 
#' To learn more about Rbedtools, start with the vignette:
#' \code{browseVignettes(package = "Rbedtools")}
#' 
#' @author Jay Hesselberth <jay.hesselberth@gmail.com>
#' 
#' @docType package
#' @name Rbedtools
#' 
#' @import dplyr
#' @import readr
#' @import tidyr
#' @importFrom assertthat assert_that
#' @importFrom purrr by_row
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @useDynLib Rbedtools
NULL
