context("tbls")

test_that("tbl_ivl classes are set", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    'chr1',  1,     50
  )
  x <- tbl_interval(x)
  expect_true(is.tbl_interval(x))
})

test_that("tbl_szs classes are set", {
  genome <- tibble::tribble(
    ~chrom, ~size,
    'chr1',  1e4
  )
  genome <- tbl_sizes(genome)
  expect_true(is.tbl_sizes(genome))
})

test_that("invalid tbl_ivl names throw error", {
  x <- tibble::tribble(
    ~pork, ~pie, ~hat,
    'chr1',  1,  50
  )
  expect_error(tbl_interval(x), 'expected 3 required names, missing: chrom, start, end')

  # missing 1 only
  x <- tibble::tribble(
    ~chrom, ~start, ~oops,
    'chr1',  1,  50
  )
  expect_error(tbl_interval(x), 'expected 3 required names, missing: end')
})

test_that("invalid tbl_szs names throw error", {
  genome <- tibble::tribble(
    ~foo,   ~bar,
    'chr1', 1e4
  )
  expect_error(tbl_sizes(genome), 'expected 2 required names, missing: chrom, size')
})

test_that("duplicate tbl_szs refs throw error", {
  genome <- tibble::tribble(
    ~chrom, ~size,
    'chr1', 1e4,
    'chr1', 1e4
  )
  expect_error(tbl_sizes(genome), 'duplicate chroms in genome: chr1')
})
