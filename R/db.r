#' Fetch data from remote databases.
#' 
#' Currently \code{db_ucsc} and \code{db_ensembl} are available for connections.
#' 
#' @name db
#' 
#' @param dbname name of database
#' @param host hostname
#' @param user username
#' @param password password
#' @param port MySQL connection port
#' @param ... params for connection
NULL

#' @rdname db
#' @seealso \url{https://genome.ucsc.edu/goldenpath/help/mysql.html}
#'   
#' @examples
#' \dontrun{
#' ucsc <- db_ucsc('hg38')
#' 
#' # fetch the `refGene` tbl
#' tbl(ucsc, "refGene")
#' }
#' 
#' @export
db_ucsc <- function(dbname, host = 'genome-mysql.cse.ucsc.edu',
                    user = 'genomep', password = 'password',
                    port = 3306, ...) {
  src_mysql(dbname, host, port, user, password, ...)
}

#' @rdname db
#' @seealso \url{http://www.ensembl.org/info/data/mysql.html}
#' 
#' @examples
#' \dontrun{
#' # squirrel genome
#' ensembl <- db_ensembl('spermophilus_tridecemlineatus_core_67_2')
#' 
#' tbl(ensembl, "gene") 
#' }
#' 
#' @export
db_ensembl <- function(dbname, host = 'ensembldb.ensembl.org',
                       user = 'anonymous', password = '',
                       port = 3306, ...) {
  src_mysql(dbname, host, port, user, password, ...)
}
