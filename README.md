
valr <img src="logo.png" align="right" />
=========================================

[![Build Status](https://travis-ci.org/jayhesselberth/valr.svg?branch=master)](https://travis-ci.org/jayhesselberth/valr) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/jayhesselberth/valr?branch=master&svg=true)](https://ci.appveyor.com/project/jayhesselberth/valr) [![Coverage Status](https://img.shields.io/codecov/c/github/jayhesselberth/valr/master.svg)](https://codecov.io/github/jayhesselberth/valr?branch=master) [![](http://www.r-pkg.org/badges/version/valr)](http://www.r-pkg.org/pkg/valr)

**`valr` provides tools to read and manipulate genome intervals and signals**, similar to the [`BEDtools`](http://bedtools.readthedocs.org/en/latest/) suite. `valr` enables BEDtools-like analysis in the R/RStudio environment, and leverages modern R tools for a terse, expressive syntax. Compute-intensive algorithms are implemented in [`Rcpp`](http://www.rcpp.org/)/C++, and many methods take advantage of the speed and grouping capability provided by [`dplyr`](https://github.com/hadley/dplyr).

Installation
------------

The latest development version can be installed from github:

``` r
# install.packages("devtools")
devtools::install_github('jayhesselberth/valr')
```

Why `valr`?
-----------

**Why another tool set for interval manipulations?** Based on our experience teaching genome analysis, we were motivated to develop interval arithmetic software that faciliates genome analysis in a single environment (RStudio), eliminating the need to master the command-line and exploratory analysis tools.

`valr` can currently be used for analysis of pre-processed data in BED and related formats. We plan to support BAM and VCF files soon via tabix indexes.

### Familiar tools, all within R

The functions in `valr` have similar names to their `BEDtools` counterparts, and so will be familiar to users coming from the `BEDtools` suite. Similar to [`pybedtools`](https://daler.github.io/pybedtools/#why-pybedtools), `valr` has a terse syntax:

``` r
library(valr)
library(dplyr)

snps <- read_bed(valr_example('hg19.snps147.chr22.bed.gz'), n_fields = 6)
genes <- read_bed(valr_example('genes.hg19.chr22.bed.gz'), n_fields = 6)

# find snps in intergenic regions
intergenic <- bed_subtract(snps, genes)
# find distance from intergenic snps to nearest gene
nearby <- bed_closest(intergenic, genes)

nearby %>%
  select(starts_with('name'), .overlap, .dist) %>%
  filter(abs(.dist) < 5000)
```

### Remote databases

Remote databases can be accessed with `db_ucsc()` (to access the UCSC Browser) and `db_ensembl()` (to access Ensembl databases).

``` r
# access the `refGene` tbl on the `hg38` assembly
ucsc <- db_ucsc('hg38')
tbl(ucsc, 'refGene')
```

### Visual documentation

`valr` includes helpful glyphs to illustrate the results of specific operations, similar to those found in the `BEDtools` documentation. This glyph (generated using `bed_glyph()`) shows the result of intersecting `x` and `y` intervals with `bed_intersect()`:

![](img/README-intersect_glyph-1.png)

### Reproducible reports

`valr` can be used in RMarkdown documents to generate reproducible work-flows for data processing. Because `valr` is reasonably fast, it can be for exploratory analysis with `RMarkdown`, and for interactive analysis using `shiny`.

API
---

Function names are similar to their their [BEDtools](http://bedtools.readthedocs.org/en/latest/) counterparts, with some additions.

### Reading data

-   BED and related files are read with `read_bed()`, `read_bed12()`, `read_bedgraph()`, `read_narrowpeak()` and `read_broadpeak()`.

-   Genome files containing chromosome name and size information are loaded with `read_genome()`.

-   VCF files are loaded with `read_vcf()`.

-   Remote databases can be accessed with `db_ucsc()` and `db_ensembl()`.

### Transforming single interval sets

-   Intervals are ordered with `bed_sort()`.

-   Interval coordinates are adjusted with `bed_slop()` and `bed_shift()`, and new flanking intervals are created with `bed_flank()`.

-   Nearby intervals are combined with `bed_merge()` and identified (but not merged) with `bed_cluster()`.

-   Intervals not covered by a query are created with `bed_complement()`.

### Comparing multiple interval sets

-   Find overlaps between two sets of intervals with `bed_intersect()`.

-   Apply functions to selected columns for overlapping intervals with `bed_map()`.

-   Remove intervals based on overlaps between two files with `bed_subtract()`.

-   Find overlapping intervals within a window with `bed_window()`.

-   Find the closest intervals independent of overlaps with `bed_closest()`.

### Randomizing intervals

-   Generate random intervals from an input genome with `bed_random()`.

-   Shuffle the coordinates of input intervals with `bed_shuffle()`.

-   Random sampling of input intervals is done with the `sample_` function family in `dplyr`.

### Interval statistics

-   Calculate significance of overlaps between two sets of intervals with `bed_fisher()` and `bed_projection()`.

-   Quantify relative and absolute distances between sets of intervals with `bed_reldist()` and `bed_absdist()`.

-   Quantify extent of overlap between two sets of intervals with `bed_jaccard()`.

### Utilities

-   Visualize the actions of valr functions with `bed_glyph()`.

-   Constrain intervals to a genome reference with `bound_intervals()`.

-   Subdivide intervals with `bed_makewindows()`.

-   Convert BED12 to BED6 format with `bed12_to_exons()`.

-   Calculate spacing between intervals with `interval_spacing()`.

Related work
------------

-   Command-line tools [BEDtools](http://bedtools.readthedocs.org/en/latest/) and [bedops](http://bedops.readthedocs.org/en/latest/index.html).

-   The Python library [pybedtools](https://pythonhosted.org/pybedtools/) wraps BEDtools.

-   The R packages [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html), [bedr](https://CRAN.R-project.org/package=bedr), [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html) and [GenometriCorr](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002529) provide similar capability with a different philosophy.
