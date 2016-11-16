#' Fetch data from remote databases
#' 
#' Currently \code{db_ucsc} is available for connections.
#' 
#' @param dbname name of database
#' @param host hostname
#' @param user username
#' @param password password
#' @param port MySQL connection port
#' @param ... params for connection
#' 
#' @examples
#' ucsc <- db_ucsc('hg38')
#' 
#' # fetch the `refGene` tbl
#' tbl(ucsc, "refGene")
#' 
#' @export
db_ucsc <- function(dbname, host = 'genome-mysql.cse.ucsc.edu',
                    user = 'genomep', password = 'password',
                    port = 3306, ...) {
  src_mysql(dbname, host, port, user, password, ...)
}

#' Convert a remote tbl to BED format
#'
as_bed <- function(x, n = 6) {
  select(x, chrom, start = txStart, end = txEnd,
         name, strand) 
}
