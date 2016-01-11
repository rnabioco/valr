read
====

``` r
library(Rbedtools)

bed_tbl <- read_bed('inst/extdata/3fields.bed.gz')
bedgraph_tbl <- read_bedgraph('inst/extdata/test.bg.gz')

bed_tbl
bedgraph_tbl
```

intersect
=========

``` r
# bed_tbl = "A" file
# bedgraph_tbl = "B" file
# A intersect B == B regions that intersect A intervals
bed_tbl %>%
  intersect(bedgraph_tbl) %>%
  summarize(n_insersects = n())
```
