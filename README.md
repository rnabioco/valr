
# valr <img src="man/figures/logo.png" align="right" />

[![Build
Status](https://travis-ci.org/rnabioco/valr.svg?branch=master)](https://travis-ci.org/rnabioco/valr)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/rnabioco/valr?branch=master&svg=true)](https://ci.appveyor.com/project/jayhesselberth/valr)
[![Coverage
Status](https://img.shields.io/codecov/c/github/rnabioco/valr/master.svg)](https://codecov.io/github/rnabioco/valr?branch=master)
[![](https://www.r-pkg.org/badges/version/valr)](https://CRAN.R-project.org/package=valr)
[![](https://cranlogs.r-pkg.org/badges/valr?color=FFD700)](https://www.r-pkg.org/pkg/valr)
[![DOI](https://zenodo.org/badge/49370633.svg)](https://zenodo.org/badge/latestdoi/49370633)

**`valr` provides tools to read and manipulate genome intervals and
signals**, similar to the \[`BEDtools`\]\[1\] suite. `valr` enables
analysis in the R/RStudio environment, leveraging modern R tools in the
\[`tidyverse`\]\[14\] for a terse, expressive syntax. Compute-intensive
algorithms are implemented in \[`Rcpp`\]\[3\]/C++, and many methods take
advantage of the speed and grouping capability provided by
\[`dplyr`\]\[2\]. See `vignette(valr)` for more details.

## Installation

The latest stable version can be installed from CRAN:

``` r
install.packages('valr')
```

The latest development version can be installed from github:

``` r
# install.packages("devtools")
devtools::install_github('rnabioco/valr')
```
