context('bed_random')

genome <- read_genome(valr_example('hg19.chrom.sizes.gz'))

# Seed for reproducible bed_random tests
seed <- 1010486

test_that("returns correct number of intervals", {
  res <- bed_random(genome, n = 1e2, seed = seed)
  expect_equal(nrow(res), 1e2)
})

test_that("returns correctly sized intervals", {
  len <- 1000
  res <- bed_random(genome, length = len, n = 1e2, seed = seed)
  expect_true(all(res$end - res$start == len))
})

test_that("all ends are less or equal to than chrom size", {
  len <- 1000
  res <- bed_random(genome, length = len, n = 1e4, seed = seed) %>%
    mutate(chrom = as.character(chrom)) %>%
    left_join(genome, by = 'chrom')
  expect_true(all(res$end <= res$size))
})

test_that("chrom sizes less than length throws an error", {
  genome <- tibble::tribble(
    ~chrom, ~size,
    'chr1', 125
  )
  expect_error(bed_random(genome, seed = seed))
})

test_that("intervals are sorted by default", {
  x <- bed_random(genome, n = 1e4, seed = seed)
  y <- bed_random(genome, n = 1e4, sort_by = NULL, seed = seed)
  expect_false(all(x == y))

  # default sort
  x_sort <- x %>% arrange(chrom, start)
  expect_true(all(x == x_sort))
})
