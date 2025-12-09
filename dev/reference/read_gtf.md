# Import and convert a GTF/GFF file into a valr compatible bed tbl format

**\[deprecated\]**

This function will output a tibble with the required chrom, start, and
end columns, as well as other columns depending on content in GTF/GFF
file.

## Usage

``` r
read_gtf(path, zero_based = TRUE)
```

## Arguments

- path:

  path to gtf or gff file

- zero_based:

  if TRUE, convert to zero based

## Examples

``` r
if (FALSE) { # \dontrun{
gtf <- read_gtf(valr_example("hg19.gencode.gtf.gz"))
head(gtf)
} # }
```
