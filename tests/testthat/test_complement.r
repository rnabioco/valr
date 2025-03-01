test_that("complement with covering interval", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 500
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 500
  )
  res <- bed_complement(bed, genome)
  expect_equal(nrow(res), 0)
})

test_that("complement with middle interval", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 500
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 200
  )
  # fmt: skip
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 100,
    "chr1", 200, 500
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("complement adds final interval", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 500
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 50,
    "chr1", 100, 200
  )
  # fmt: skip
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 50, 100,
    "chr1", 200, 500
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("multiple chroms", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 500,
    "chr2", 500
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 500,
    "chr2", 100, 500
  )
  # fmt: skip
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 100,
    "chr2", 0, 100
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("multiple chroms, chr1 is covered", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 500,
    "chr2", 500
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 500,
    "chr2", 100, 500
  )
  # fmt: skip
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr2", 0, 100
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("record exceeds chrom length", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 100
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 90, 110
  )
  # fmt: skip
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 90
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("test complement only on left side", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 100
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 90, 100
  )
  # fmt: skip
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 90
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("overlapping intervals", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 100
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 10, 100,
    "chr1", 20, 80
  )
  # fmt: skip
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 10
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("starting intervals are complemented (#78)", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 100
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 20
  )
  res <- bed_complement(bed, genome)
  expect_equal(nrow(res), 1)
  expect_equal(res$start, 20)
  expect_equal(res$end, 100)
})

test_that("non-overlapping genome chroms are in output (#78)", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 100,
    "chr2", 100
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 10, 100
  )
  # fmt: skip
  expected <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 10,
    "chr2", 0, 100
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expected))
})

test_that("`x` intervals not in genome are ignored (#78)", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 500
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr2", 10, 100
  )
  # fmt: skip
  expected <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 500
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expected))
})

# from https://github.com/arq5x/bedtools2/blob/master/test/complement/test-complement.sh
test_that("ends are covered", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 20
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 1,
    "chr1", 19, 20
  )
  # fmt: skip
  expected <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1, 19
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expected))
})

test_that("entirety is covered", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 20
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 20
  )
  res <- bed_complement(bed, genome)
  expect_true(nrow(res) == 0)
})

test_that("nothing is covered", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 20,
    "chr2", 20
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr2", 0, 20
  )
  # fmt: skip
  expected <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 20
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expected))
})

test_that("Issue #356", {
  # fmt: skip
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 249250621
  )
  # fmt: skip
  bed <- tibble::tribble(
    ~chrom, ~start, ~end, ~name,
    "chr1", 0, 10000, "telomere"
  )
  # fmt: skip
  expected <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 10000, 249250621
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expected))
})
