context('bed_random')

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 5000,
  "chr2", 1e6
) %>% mutate(chrom = as.character(chrom))

# Seed for reproducible bed_random tests
seed <- 1010486

test_that("returns correct number of intervals", {
  res <- bed_random(genome, n = 1e6, seed = seed)
  expect_equal(nrow(res), 1e6)
})

test_that("returns correctly sized intervals", {
  len <- 1000
  res <- bed_random(genome, length = len, n = 1e6, seed = seed)
  expect_true(all(res$end - res$start == len))
})

test_that("all ends are less or equal to than chrom size", {
  len <- 1000
  res <- bed_random(genome, length = len, n = 1e6, seed = seed) %>%
    mutate(chrom = as.character(chrom)) %>%
    left_join(genome, by = 'chrom')
  expect_true(all(res$end <= res$size))
})

test_that("chrom sizes less than length throws an error", {
  genome <- tibble::tribble(
    ~chrom, ~size,
    'chr1',      125
  )
  expect_error(bed_random(genome, seed = seed))
})

test_that("intervals are sorted by default", {
  x <- bed_random(genome, seed = seed)
  y <- bed_random(genome, sort_by = NULL, seed = seed)
  expect_false(all(x == y))
})
