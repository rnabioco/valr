context("bed_shuffle")

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 1e6,
  "chr2", 1e7,
  "chr3", 1e8
)

x <- bed_random(genome) %>% bed_sort

test_that('within = TRUE maintains chroms', {
  res <- bed_shuffle(x, genome, within = TRUE)
  expect_true(all(x$chrom == res$chrom))
})

test_that('within = FALSE shuffles chroms', {
  res <- bed_shuffle(x, genome, within = FALSE)
  expect_false(all(x$chrom == res$chrom))
})

test_that('`incl` includes intervals',{
  incl <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1e4,    1e6
  )
  res <- bed_shuffle(x, genome, incl = incl)
  expect_true(all(res$chrom == 'chr1'))
  expect_true(all(res$start >= 1e4))
  expect_true(all(res$end <= 1e6))
})

test_that('`excl` excludes intervals',{
  excl <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1e4,    1e6,
    "chr2", 1,      1e7,
    "chr3", 1,      1e8
  )
  res <- bed_shuffle(x, genome, excl = excl)
  expect_true(all(res$chrom == 'chr1'))
  expect_false(any(res$chrom == 'chr2'))
  expect_false(any(res$chrom == 'chr3'))
  expect_true(all(res$start < 1e4))
})

test_that('completely excluded intervals throw an error',{
  # all intervals completely excluded 
  excl <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1,      1e6,
    "chr2", 1,      1e7,
    "chr3", 1,      1e8
  )
  expect_error(bed_shuffle(x, genome, excl = excl))
})

test_that('exceeding `max_tries` yields an error',{
  # intervals partialy excluded 
  excl <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100,    1e6,
    "chr2", 1,      1e7,
    "chr3", 1,      1e8
  )
  expect_error(bed_shuffle(x, genome, excl = excl))
})

test_that('`seed` generates reproducible intervals',{
   seed <- 1010486
   res1 <- bed_shuffle(x, genome, seed = seed)
   res2 <- bed_shuffle(x, genome, seed = seed)
   expect_identical(res1, res2) 
})
