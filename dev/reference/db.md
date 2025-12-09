# Fetch data from remote databases.

Currently `db_ucsc` and `db_ensembl` are available for connections.

## Usage

``` r
db_ucsc(
  dbname,
  host = "genome-mysql.cse.ucsc.edu",
  user = "genomep",
  password = "password",
  port = 3306,
  ...
)

db_ensembl(
  dbname,
  host = "ensembldb.ensembl.org",
  user = "anonymous",
  password = "",
  port = 3306,
  ...
)
```

## Arguments

- dbname:

  name of database

- host:

  hostname

- user:

  username

- password:

  password

- port:

  MySQL connection port

- ...:

  params for connection

## See also

<https://genome.ucsc.edu/goldenpath/help/mysql.html>

<https://www.ensembl.org/info/data/mysql.html>

## Examples

``` r
if (FALSE) { # \dontrun{
if (require(RMariaDB)) {
  library(dplyr)
  ucsc <- db_ucsc("hg38")

  # fetch the `refGene` tbl
  tbl(ucsc, "refGene")

  # the `chromInfo` tbls have size information
  tbl(ucsc, "chromInfo")
}
} # }
if (FALSE) { # \dontrun{
if (require(RMariaDB)) {
  library(dplyr)
  # squirrel genome
  ensembl <- db_ensembl("spermophilus_tridecemlineatus_core_67_2")

  tbl(ensembl, "gene")
}
} # }
```
