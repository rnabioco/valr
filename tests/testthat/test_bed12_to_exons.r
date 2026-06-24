test_that("invalid ncol causes an error", {
  x <- tibble::tribble()
  expect_snapshot(bed12_to_exons(x), error = TRUE)
})

test_that("non-BED12 12-column input errors clearly", {
  # 12 columns but the exon columns are not comma-separated character values
  x <- tibble::as_tibble(
    setNames(
      c(
        list(chrom = "chr1", start = 1L, end = 100L),
        rep(list(0L), 9)
      ),
      c("chrom", "start", "end", letters[1:9])
    )
  )
  expect_error(bed12_to_exons(x), "does not look like BED12")
})

test_that("BED12 is parsed correctly", {
  x <- read_bed12(valr_example("mm9.refGene.bed.gz"))
  expect_equal(nrow(bed12_to_exons(x)), 1683)
  expect_equal(ncol(bed12_to_exons(x)), 6)
})

test_that("BED12 with UCSC column names (e.g. read_bigbed) is accepted (#461)", {
  bb <- system.file("extdata", "test.bb", package = "cpp11bigwig")
  x <- cpp11bigwig::read_bigbed(bb)

  res <- bed12_to_exons(x)
  expect_named(res, c("chrom", "start", "end", "name", "score", "strand"))
  expect_equal(ncol(res), 6)
  expect_equal(nrow(res), sum(x$blockCount))

  # positional rename matches valr-named input exactly
  y <- x
  names(y) <- c(
    "chrom", "start", "end", "name", "score", "strand",
    "cds_start", "cds_end", "item_rgb", "exon_count", "exon_sizes", "exon_starts"
  )
  expect_equal(res, bed12_to_exons(y))

  # valr-named columns are looked up by name, so reordering them is preserved
  # (positional renaming is only applied to foreign sources)
  expect_equal(res, bed12_to_exons(y[, rev(names(y))]))
})
