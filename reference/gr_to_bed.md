# Convert Granges to bed tibble

Convert Granges to bed tibble

## Usage

``` r
gr_to_bed(x)
```

## Arguments

- x:

  GRanges object to convert to bed tibble.

## Value

[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)

## Examples

``` r
if (FALSE) { # \dontrun{
gr <- GenomicRanges::GRanges(
  seqnames = S4Vectors::Rle(
    c("chr1", "chr2", "chr1", "chr3"),
    c(1, 1, 1, 1)
  ),
  ranges = IRanges::IRanges(
    start = c(1, 10, 50, 100),
    end = c(100, 500, 1000, 2000),
    names = head(letters, 4)
  ),
  strand = S4Vectors::Rle(
    c("-", "+"), c(2, 2)
  )
)

gr_to_bed(gr)

# There are two ways to convert a bed-like data.frame to GRanges:

gr <- GenomicRanges::GRanges(
  seqnames = S4Vectors::Rle(x$chrom),
  ranges = IRanges::IRanges(
    start = x$start + 1,
    end = x$end,
    names = x$name
  ),
  strand = S4Vectors::Rle(x$strand)
)
# or:

gr <- GenomicRanges::makeGRangesFromDataFrame(dplyr::mutate(x, start = start + 1))
} # }
```
