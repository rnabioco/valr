#' Fetch data from remote databases.
#'
#' Currently `db_ucsc` and `db_ensembl` are available for connections.
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
#' if(require(RMariaDB)) {
#'   ucsc <- db_ucsc('hg38')
#'
#'   # fetch the `refGene` tbl
#'   tbl(ucsc, "refGene")
#'
#'   # the `chromInfo` tbls have size information
#'   tbl(ucsc, "chromInfo")
#' }
#' }
#' @importFrom DBI dbConnect
#' @importFrom RMariaDB MariaDB
#' @export
db_ucsc <- function(dbname, host = "genome-mysql.cse.ucsc.edu",
                    user = "genomep", password = "password",
                    port = 3306, ...) {
  DBI::dbConnect(RMariaDB::MariaDB(),
                 dbname = dbname,
                 user = user,
                 password = password,
                 host = host,
                 post = port, ...) # nocov
}

#' @rdname db
#' @seealso \url{http://www.ensembl.org/info/data/mysql.html}
#'
#' @examples
#' \dontrun{
#' if(require(RMariaDB)) {
#'   # squirrel genome
#'   ensembl <- db_ensembl('spermophilus_tridecemlineatus_core_67_2')
#'
#'   tbl(ensembl, "gene")
#' }
#' }
#'
#' @export
db_ensembl <- function(dbname, host = "ensembldb.ensembl.org",
                       user = "anonymous", password = "",
                       port = 3306, ...) {
  DBI::dbConnect(RMariaDB::MariaDB(),
                 dbname = dbname,
                 user = user,
                 password = password,
                 host = host,
                 post = port, ...) # nocov # nocov
}
