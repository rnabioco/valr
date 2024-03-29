% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bed_coverage.R
\name{bed_coverage}
\alias{bed_coverage}
\title{Compute coverage of intervals.}
\usage{
bed_coverage(x, y, ...)
}
\arguments{
\item{x}{\link{ivl_df}}

\item{y}{\link{ivl_df}}

\item{...}{extra arguments (not used)}
}
\value{
\link{ivl_df} with the following additional columns:
\itemize{
\item \code{.ints} number of \code{x} intersections
\item \code{.cov} per-base coverage of \code{x} intervals
\item \code{.len} total length of \code{y} intervals covered by \code{x} intervals
\item \code{.frac} \code{.len} scaled by the number of \code{y} intervals
}
}
\description{
Compute coverage of intervals.
}
\details{
input tbls are grouped by \code{chrom} by default, and additional
groups can be added using \code{\link[dplyr:group_by]{dplyr::group_by()}}. For example,
grouping by \code{strand} will constrain analyses to the same strand. To
compare opposing strands across two tbls, strands on the \code{y} tbl can
first be inverted using \code{\link[=flip_strands]{flip_strands()}}.
}
\note{
Book-ended intervals are included in coverage calculations.
}
\examples{
x <- tibble::tribble(
  ~chrom, ~start, ~end, ~strand,
  "chr1", 100,    500,  "+",
  "chr2", 200,    400,  "+",
  "chr2", 300,    500,  "-",
  "chr2", 800,    900,  "-"
)

y <- tibble::tribble(
  ~chrom, ~start, ~end, ~value, ~strand,
  "chr1", 150,    400,  100,    "+",
  "chr1", 500,    550,  100,    "+",
  "chr2", 230,    430,  200,    "-",
  "chr2", 350,    430,  300,    "-"
)

bed_coverage(x, y)

}
\seealso{
\url{https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html}

Other multiple set operations: 
\code{\link{bed_closest}()},
\code{\link{bed_intersect}()},
\code{\link{bed_map}()},
\code{\link{bed_subtract}()},
\code{\link{bed_window}()}
}
\concept{multiple set operations}
