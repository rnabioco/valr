x <- tibble::tribble(
  ~chrom, ~start, ~end, ~strand,
  "chr1", 20, 70, "+",
  "chr1", 50, 100, "-",
  "chr1", 200, 250, "+",
  "chr1", 220, 250, "+"
)

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 500,
  "chr2", 1000
)
test_that("bed_genomecov works", {
  ex <- tibble::tribble(
    ~chrom, ~start, ~end, ~.depth,
    "chr1", 20L, 50L, 1L,
    "chr1", 50L, 70L, 2L,
    "chr1", 70L, 100L, 1L,
    "chr1", 200L, 220L, 1L,
    "chr1", 220L, 250L, 2L
  )
  obs <- bed_genomecov(x, genome)
  expect_identical(obs, ex)

  res <- bed_genomecov(x, genome, zero_depth = TRUE)
  expect_true(sum(res$.depth == 0) > 0)
})

test_that("groups are respected", {
  ex <- tibble::tribble(
    ~chrom, ~start, ~end, ~strand, ~.depth,
    "chr1", 20L, 70L, "+", 1L,
    "chr1", 200L, 220L, "+", 1L,
    "chr1", 220L, 250L, "+", 2L,
    "chr1", 50L, 100L, "-", 1L
  )
  obs <- bed_genomecov(group_by(x, strand), genome)
  expect_identical(obs, ex)

  res <- bed_genomecov(group_by(x, strand), genome, zero_depth = TRUE)
  expect_equal(sum(res$chrom == "chr2"), 2L)
})

test_that("grouping is retained for zero depth intervals", {
  xx <- tibble::tribble(
    ~chrom, ~start, ~end, ~strand, ~grp,
    "chr1", 20L, 70L, "+", 1,
    "chr1", 200L, 220L, "-", 1,
    "chr1", 20L, 70L, "+", 2,
    "chr1", 200L, 220L, "-", 2
  ) |>
    group_by(strand, grp)

  many_chroms_genome <- tibble(
    chrom = c("chr1", LETTERS),
    size  = 500
  )

  res <- bed_genomecov(xx, many_chroms_genome, zero_depth = TRUE)
  expect_equal(length(setdiff(many_chroms_genome$chrom, res$chrom)), 0L)

  ex <- tibble::tribble(
    ~strand, ~grp, ~n,
    "+", 1, 26L,
    "+", 2, 26L,
    "-", 1, 26L,
    "-", 2, 26L
  )

  lr <- res[res$chrom %in% LETTERS, ]
  nlr <- lr |>
    group_by(strand, grp) |>
    summarize(n = n()) |>
    ungroup()
  expect_identical(nlr, ex)
})

test_that("chroms in bed, not in genome, are dropped", {
  xx <- tibble::tribble(
    ~chrom, ~start, ~end, ~strand, ~.depth,
    "hello", 20L, 70L, "+", 1L,
    "world", 200L, 220L, "+", 1L,
    "chr1", 220L, 250L, "+", 2L
  )
  expect_warning(res <- bed_genomecov(xx, genome))
  expect_true(all(res$chrom == "chr1"))
})

test_that("zero length input is handled", {
  xx <- tibble::tribble(
    ~chrom, ~start, ~end, ~strand, ~.depth,
    "hello", 20L, 70L, "+", 1L,
  )

  expect_warning(res <- bed_genomecov(xx, genome))
  expect_true(nrow(res) == 0)

  xx <- xx[xx$chrom != "hello", ]

  res <- bed_genomecov(xx, genome)
  expect_true(nrow(res) == 0)
})

test_that("check edge cases with 1 bp intervals", {
  # base-level coverage equals number of basepairs in input intervals
  genome <- tribble(
    ~chrom, ~size,
    "chr1", 1e5
  )
  seed <- 1010486
  ivls <- bed_random(genome, length = 1, n = 1e3, seed = seed)

  res <- bed_genomecov(ivls, genome)
  expect_true(sum(res$.depth) == 1e3)

  ivls <- tibble(
    chrom = "chr1",
    start = seq(0, 999),
    end = start + 2
  )
  res <- bed_genomecov(ivls, genome)
  expect_true(sum(res$.depth) == (1e3 * 2))

  set.seed(seed)
  ivls <- tibble(
    chrom = "chr1",
    start = seq(0, 999),
    end = start + sample(1:100, length(start), replace = TRUE)
  )

  res <- bed_genomecov(ivls, genome)
  n_bp <- sum(ivls$end - ivls$start)
  n_cov <- sum(res$.depth * (res$end - res$start))
  expect_equal(n_bp, n_cov)
})

test_that("check edge cases at beginning and end", {
  genome <- tribble(
    ~chrom, ~size,
    "chr1", 1000
  )
  ex <- tibble::tribble(
    ~chrom, ~start, ~end, ~.depth,
    "chr1", 0L, 1L, 3L,
    "chr1", 1L, 2L, 1L,
    "chr1", 999L, 1000L, 1L
  )

  # oob intervals are ignored with a warning
  ivls <- tibble(
    chrom = "chr1",
    start = c(rep(0, 3), 1, 999, 1000),
    end = start + 1
  )

  expect_warning(res <- bed_genomecov(ivls, genome))
  expect_true(all(res$start < 1000))

  expect_identical(res, ex)
})

# bed related tests from #https://github.com/arq5x/bedtools2/blob/master/test/genomecov/test-genomecov.sh

y <- tibble::tribble(
  ~chrom, ~start, ~end, ~group, ~score, ~strand,
  "1", 15L, 20L, "y1", 1L, "+",
  "1", 17L, 22L, "y2", 2L, "+"
)

genome <- tibble::tribble(
  ~chrom, ~size,
  "1", 100L,
  "2", 100L,
  "3", 100L
)

test_that("Test with chroms that have no coverage", {
  ex <- tibble::tribble(
    ~chrom, ~start, ~end, ~.depth,
    "1", 15L, 17L, 1L,
    "1", 17L, 20L, 2L,
    "1", 20L, 22L, 1L
  )
  obs <- bed_genomecov(y, genome)
  expect_equal(ex, obs)
})

test_that("Test with chroms that have no coverage", {
  ex <- tibble::tribble(
    ~chrom, ~start, ~end, ~.depth,
    "1", 0L, 15L, 0L,
    "1", 15L, 17L, 1L,
    "1", 17L, 20L, 2L,
    "1", 20L, 22L, 1L,
    "1", 22L, 100L, 0L,
    "2", 0L, 100L, 0L,
    "3", 0L, 100L, 0L
  )

  obs <- bed_genomecov(y, genome, zero_depth = TRUE)
  expect_equal(ex, obs)
})
