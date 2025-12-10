# read gtf

    Code
      read_gtf(gtf_path)
    Condition
      Error:
      ! `read_gtf()` was deprecated in valr 0.8.3 and is now defunct.
      x read_gtf() was removed because rtracklayer does not pass CRAN AddressSantizer checks of the UCSC C-library code vendored in rtracklayer.
      i convert GTF to BED, and then `read_bed()`, or use `rtracklayer::import()` then `gr_to_bed()`.

