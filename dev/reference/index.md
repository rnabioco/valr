# Package index

## Data types

- [`tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
  [`tibble_row()`](https://tibble.tidyverse.org/reference/tibble.html) :
  Build a data frame (from tibble)
- [`tribble()`](https://tibble.tidyverse.org/reference/tribble.html) :
  Row-wise tibble creation (from tibble)

## Load data

- [`read_bed()`](https://rnabioco.github.io/valr/dev/reference/read_bed.md)
  [`read_bed12()`](https://rnabioco.github.io/valr/dev/reference/read_bed.md)
  [`read_bedgraph()`](https://rnabioco.github.io/valr/dev/reference/read_bed.md)
  [`read_narrowpeak()`](https://rnabioco.github.io/valr/dev/reference/read_bed.md)
  [`read_broadpeak()`](https://rnabioco.github.io/valr/dev/reference/read_bed.md)
  : Read BED and related files.
- [`read_genome()`](https://rnabioco.github.io/valr/dev/reference/read_genome.md)
  : Read genome files.
- [`read_gtf()`](https://rnabioco.github.io/valr/dev/reference/read_gtf.md)
  **\[deprecated\]** : Import and convert a GTF/GFF file into a valr
  compatible bed tbl format
- [`read_vcf()`](https://rnabioco.github.io/valr/dev/reference/read_vcf.md)
  : Read a VCF file.
- [`db_ucsc()`](https://rnabioco.github.io/valr/dev/reference/db.md)
  [`db_ensembl()`](https://rnabioco.github.io/valr/dev/reference/db.md)
  : Fetch data from remote databases.
- [`read_bigwig()`](https://rnabioco.github.io/cpp11bigwig/reference/read_bigwig.html)
  : Read data from bigWig files. (from cpp11bigwig)
- [`read_bigbed()`](https://rnabioco.github.io/cpp11bigwig/reference/read_bigbed.html)
  : Read data from bigBed files. (from cpp11bigwig)

## Single set operations

- [`bed_cluster()`](https://rnabioco.github.io/valr/dev/reference/bed_cluster.md)
  : Cluster neighboring intervals.
- [`bed_complement()`](https://rnabioco.github.io/valr/dev/reference/bed_complement.md)
  : Identify intervals in a genome not covered by a query.
- [`bed_flank()`](https://rnabioco.github.io/valr/dev/reference/bed_flank.md)
  : Create flanking intervals from input intervals.
- [`bed_genomecov()`](https://rnabioco.github.io/valr/dev/reference/bed_genomecov.md)
  : Calculate coverage across a genome
- [`bed_merge()`](https://rnabioco.github.io/valr/dev/reference/bed_merge.md)
  : Merge overlapping intervals.
- [`bed_partition()`](https://rnabioco.github.io/valr/dev/reference/bed_partition.md)
  : Partition intervals into elemental intervals
- [`bed_slop()`](https://rnabioco.github.io/valr/dev/reference/bed_slop.md)
  : Increase the size of input intervals.
- [`bed_shift()`](https://rnabioco.github.io/valr/dev/reference/bed_shift.md)
  : Adjust intervals by a fixed size.
- [`bed_sort()`](https://rnabioco.github.io/valr/dev/reference/bed_sort.md)
  : Sort a set of intervals.

## Multiple set operations

- [`bed_closest()`](https://rnabioco.github.io/valr/dev/reference/bed_closest.md)
  : Identify closest intervals.
- [`bed_coverage()`](https://rnabioco.github.io/valr/dev/reference/bed_coverage.md)
  : Compute coverage of intervals.
- [`bed_intersect()`](https://rnabioco.github.io/valr/dev/reference/bed_intersect.md)
  : Identify intersecting intervals.
- [`bed_map()`](https://rnabioco.github.io/valr/dev/reference/bed_map.md)
  [`concat()`](https://rnabioco.github.io/valr/dev/reference/bed_map.md)
  [`values_unique()`](https://rnabioco.github.io/valr/dev/reference/bed_map.md)
  [`values()`](https://rnabioco.github.io/valr/dev/reference/bed_map.md)
  : Calculate summaries from overlapping intervals.
- [`bed_subtract()`](https://rnabioco.github.io/valr/dev/reference/bed_subtract.md)
  : Subtract two sets of intervals.
- [`bed_window()`](https://rnabioco.github.io/valr/dev/reference/bed_window.md)
  : Identify intervals within a specified distance.

## Randomizing intervals

- [`bed_random()`](https://rnabioco.github.io/valr/dev/reference/bed_random.md)
  : Generate randomly placed intervals on a genome.
- [`bed_shuffle()`](https://rnabioco.github.io/valr/dev/reference/bed_shuffle.md)
  : Shuffle input intervals.

## Interval statistics

- [`bed_absdist()`](https://rnabioco.github.io/valr/dev/reference/bed_absdist.md)
  : Compute absolute distances between intervals.
- [`bed_reldist()`](https://rnabioco.github.io/valr/dev/reference/bed_reldist.md)
  : Compute relative distances between intervals.
- [`bed_fisher()`](https://rnabioco.github.io/valr/dev/reference/bed_fisher.md)
  : Fisher's test to measure overlap between two sets of intervals.
- [`bed_jaccard()`](https://rnabioco.github.io/valr/dev/reference/bed_jaccard.md)
  : Calculate the Jaccard statistic for two sets of intervals.
- [`bed_projection()`](https://rnabioco.github.io/valr/dev/reference/bed_projection.md)
  : Projection test for query interval overlap.

## Utilities

- [`create_introns()`](https://rnabioco.github.io/valr/dev/reference/create_introns.md)
  : Create intron features.
- [`create_tss()`](https://rnabioco.github.io/valr/dev/reference/create_tss.md)
  : Create transcription start site features.
- [`create_utrs3()`](https://rnabioco.github.io/valr/dev/reference/create_utrs3.md)
  : Create 3' UTR features.
- [`create_utrs5()`](https://rnabioco.github.io/valr/dev/reference/create_utrs5.md)
  : Create 5' UTR features.
- [`bed_makewindows()`](https://rnabioco.github.io/valr/dev/reference/bed_makewindows.md)
  : Divide intervals into new sub-intervals ("windows").
- [`bed12_to_exons()`](https://rnabioco.github.io/valr/dev/reference/bed12_to_exons.md)
  : Convert BED12 to individual exons in BED6.
- [`bed_glyph()`](https://rnabioco.github.io/valr/dev/reference/bed_glyph.md)
  : Create example glyphs for valr functions.
- [`bound_intervals()`](https://rnabioco.github.io/valr/dev/reference/bound_intervals.md)
  : Select intervals bounded by a genome.
- [`flip_strands()`](https://rnabioco.github.io/valr/dev/reference/flip_strands.md)
  : Flip strands in intervals.
- [`gr_to_bed()`](https://rnabioco.github.io/valr/dev/reference/gr_to_bed.md)
  : Convert Granges to bed tibble
- [`interval_spacing()`](https://rnabioco.github.io/valr/dev/reference/interval_spacing.md)
  : Calculate interval spacing.
- [`check_interval()`](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)
  [`check_genome()`](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)
  : Bed-like data.frame requirements for valr functions
- [`valr`](https://rnabioco.github.io/valr/dev/reference/valr.md) :
  valr: genome interval arithmetic in R
- [`valr_example()`](https://rnabioco.github.io/valr/dev/reference/valr_example.md)
  : Provide working directory for valr example files.
