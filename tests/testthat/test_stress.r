# Stress tests for large-scale interval operations
# These tests verify correctness at scale and catch performance regressions

# Skip on CRAN to avoid long test times
skip_on_cran()

# Use a realistic genome for testing
genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 248956422,
  "chr2", 242193529,
  "chr3", 198295559,
  "chr4", 190214555,
  "chr5", 181538259
)

test_that("bed_intersect handles 100K x 100K intervals",
{
  x <- bed_random(genome, n = 1e5, seed = 42, length = 1000)
  y <- bed_random(genome, n = 1e5, seed = 43, length = 1000)

  res <- bed_intersect(x, y, min_overlap = 1L)

  # Verify we got results

  expect_gt(nrow(res), 0)

  # Verify all overlaps are valid (overlap > 0)
  expect_true(all(res$.overlap > 0))

  # Verify coordinates are sensible
  expect_true(all(res$start.x < res$end.x))
  expect_true(all(res$start.y < res$end.y))
})

test_that("bed_intersect handles dense overlaps (many overlaps per interval)",
{
  # Small genome = high density = many overlaps per interval
  small_genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 10000000
  )

  x <- bed_random(small_genome, n = 1e4, seed = 42, length = 5000)
  y <- bed_random(small_genome, n = 1e4, seed = 43, length = 5000)

  res <- bed_intersect(x, y, min_overlap = 1L)

  # With 10K intervals of length 5K on a 10M genome, expect high overlap density
  # Average overlaps per x interval should be > 10
  avg_overlaps <- nrow(res) / nrow(x)
  expect_gt(avg_overlaps, 10)

  # Verify all results are valid
  expect_true(all(res$.overlap > 0))
})

test_that("bed_subtract handles 100K x 100K intervals",
{
  x <- bed_random(genome, n = 1e5, seed = 42, length = 1000)
  y <- bed_random(genome, n = 1e5, seed = 43, length = 500)

  res <- bed_subtract(x, y, min_overlap = 1L)

  # Result should have fewer total bases than input
  x_bases <- sum(x$end - x$start)
  res_bases <- sum(res$end - res$start)
  expect_lt(res_bases, x_bases)

  # All result intervals should be valid
  expect_true(all(res$start < res$end))
})

test_that("bed_coverage handles 100K x 100K intervals",
{
  x <- bed_random(genome, n = 1e5, seed = 42, length = 1000)
  y <- bed_random(genome, n = 1e5, seed = 43, length = 500)

  res <- bed_coverage(x, y, min_overlap = 1L)

  # Should have same number of rows as x

  expect_equal(nrow(res), nrow(x))

  # Coverage fraction should be between 0 and 1

  expect_true(all(res$.frac >= 0 & res$.frac <= 1))

  # .len should match interval lengths
  expect_equal(res$.len, res$end - res$start)
})

test_that("bed_closest handles 100K x 100K intervals",
{
  x <- bed_random(genome, n = 1e5, seed = 42, length = 100)
  y <- bed_random(genome, n = 1e5, seed = 43, length = 100)

  res <- bed_closest(x, y)

  # Should have at least as many rows as x (could have ties)
  expect_gte(nrow(res), nrow(x))

  # Distance should be defined for all results
  expect_true(all(!is.na(res$.dist)))
})

test_that("bed_intersect with invert handles large datasets",
{
  x <- bed_random(genome, n = 1e5, seed = 42, length = 1000)
  y <- bed_random(genome, n = 1e4, seed = 43, length = 1000)

  res_overlap <- bed_intersect(x, y, min_overlap = 1L)
  res_invert <- bed_intersect(x, y, invert = TRUE, min_overlap = 1L)

  # x intervals should partition into overlapping and non-overlapping
  x_with_overlap <- length(unique(res_overlap$start.x))
  x_without_overlap <- nrow(res_invert)

  # These should account for all x intervals (approximately, due to duplicates)
  expect_lte(x_with_overlap + x_without_overlap, nrow(x) * 1.1)
})

test_that("results are consistent between grouped and ungrouped operations",
{
  x <- bed_random(genome, n = 1e4, seed = 42, length = 1000)
  y <- bed_random(genome, n = 1e4, seed = 43, length = 1000)

  # Ungrouped (internally groups by chrom)
  res_ungrouped <- bed_intersect(x, y, min_overlap = 1L)

  # Explicitly grouped by chrom
  res_grouped <- bed_intersect(
    group_by(x, chrom),
    group_by(y, chrom),
    min_overlap = 1L
  )

  # Results should be identical
  expect_equal(nrow(res_ungrouped), nrow(res_grouped))
  expect_equal(
    sort(res_ungrouped$.overlap),
    sort(res_grouped$.overlap)
  )
})
