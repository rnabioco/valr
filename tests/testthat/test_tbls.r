context("tbls")

test_that("invalid column names throw error", {
  x <- tibble::tribble(
    ~pork, ~pie, ~hat,
    "chr1", 1, 50
  )
  expect_error(check_interval(x), "expected 3 required names, missing: chrom, start, end")

  # missing 1 only
  x <- tibble::tribble(
    ~chrom, ~start, ~oops,
    "chr1", 1, 50
  )
  expect_error(check_interval(x), "expected 3 required names, missing: end")
})

test_that("invalid column names throw error", {
  genome <- tibble::tribble(
    ~foo, ~bar,
    "chr1", 1e4
  )
  expect_error(check_genome(genome), "expected 2 required names, missing: chrom, size")
})

test_that("duplicate chromosomes refs throw error", {
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 1e4,
    "chr1", 1e4
  )
  expect_error(check_genome(genome), "duplicate chroms in genome: chr1")
})

test_that("gr_to_bed coerces GRanges objects", {
  skip_if_not_installed("GenomicRanges")

  require(GenomicRanges)
  gr <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
    ranges = IRanges(1:10, end = 7:16, names = head(letters, 10)),
    strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2))
  )

  res <- gr_to_bed(gr)
  expect_silent(check_interval(res))
  expect_is(res, "tbl_df")
})
