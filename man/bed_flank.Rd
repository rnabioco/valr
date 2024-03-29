% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bed_flank.r
\name{bed_flank}
\alias{bed_flank}
\title{Create flanking intervals from input intervals.}
\usage{
bed_flank(
  x,
  genome,
  both = 0,
  left = 0,
  right = 0,
  fraction = FALSE,
  strand = FALSE,
  trim = FALSE,
  ...
)
}
\arguments{
\item{x}{\link{ivl_df}}

\item{genome}{\link{genome_df}}

\item{both}{number of bases on both sizes}

\item{left}{number of bases on left side}

\item{right}{number of bases on right side}

\item{fraction}{define flanks based on fraction of interval length}

\item{strand}{define \code{left} and \code{right} based on strand}

\item{trim}{adjust coordinates for out-of-bounds intervals}

\item{...}{extra arguments (not used)}
}
\value{
\link{ivl_df}
}
\description{
Create flanking intervals from input intervals.
}
\examples{
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 25, 50,
  "chr1", 100, 125
)

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 130
)

bed_glyph(bed_flank(x, genome, both = 20))

x <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 500,    1000, ".",   ".",    "+",
  "chr1", 1000,   1500, ".",   ".",    "-"
)

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 5000
)

bed_flank(x, genome, left = 100)

bed_flank(x, genome, right = 100)

bed_flank(x, genome, both = 100)

bed_flank(x, genome, both = 0.5, fraction = TRUE)

}
\seealso{
\url{https://bedtools.readthedocs.io/en/latest/content/tools/flank.html}

Other single set operations: 
\code{\link{bed_cluster}()},
\code{\link{bed_complement}()},
\code{\link{bed_genomecov}()},
\code{\link{bed_merge}()},
\code{\link{bed_partition}()},
\code{\link{bed_shift}()},
\code{\link{bed_slop}()}
}
\concept{single set operations}
