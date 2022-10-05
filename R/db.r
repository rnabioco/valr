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
#'   library(dplyr)
#'   ucsc <- db_ucsc('hg38')
#'
#'   # fetch the `refGene` tbl
#'   tbl(ucsc, "refGene")
#'
#'   # the `chromInfo` tbls have size information
#'   tbl(ucsc, "chromInfo")
#' }
#' }
#' @export
db_ucsc <- function(dbname, host = "genome-mysql.cse.ucsc.edu",
                    user = "genomep", password = "password",
                    port = 3306, ...) {
  check_db_packages()

  DBI::dbConnect(RMariaDB::MariaDB(),
                 dbname = dbname,
                 user = user,
                 password = password,
                 host = host,
                 post = port, ...) # nocov
}

#' @rdname db
#' @seealso \url{https://www.ensembl.org/info/data/mysql.html}
#'
#' @examples
#' \dontrun{
#' if(require(RMariaDB)) {
#'   library(dplyr)
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
  check_db_packages()

  DBI::dbConnect(RMariaDB::MariaDB(),
                 dbname = dbname,
                 user = user,
                 password = password,
                 host = host,
                 post = port, ...) # nocov
}

check_db_packages <- function() {
  db_pkgs <- c("dbplyr", "DBI", "RMariaDB")

  pkgs_found <- sapply(db_pkgs, requireNamespace, quietly = TRUE)

  if (!all(pkgs_found)) {
    missing <- db_pkgs[!pkgs_found]
    n <- length(missing)
    cli::cli_abort("install {n} package{?s} {.pkg {missing}} for db functions.")
  }
}
