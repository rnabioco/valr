context("bed_shuffle")

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 1e6,
  "chr2", 1e7,
  "chr3", 1e8
)

x <- bed_random(genome, n = 100) %>% bed_sort

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
    "chr1", 10000,  1000000
  )
  res <- bed_shuffle(x, genome, incl = incl)
  expect_true(all(res$chrom == 'chr1'))
  expect_true(all(res$start >= 1e4))
  expect_true(all(res$end <= 1e6))
})

test_that('`excl` excludes intervals',{
  excl <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 10000,  1000000,
    "chr2", 1,      10000000,
    "chr3", 1,      100000000
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
    "chr1", 1,      1000000,
    "chr2", 1,      10000000,
    "chr3", 1,      100000000
  )
  expect_error(bed_shuffle(x, genome, excl = excl))
})

test_that('`incl` and `excl` are handled', {
  excl <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1,      500000,
    "chr2", 1,      10000000
  )
  incl <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1,  1000000
  )
  res <- bed_shuffle(x, genome, incl, excl)
  expect_true(all(res$chrom == 'chr1'))
  expect_true(all(res$start > 500000))
})

test_that('empty intervals derived from `incl` and `excl` is handled', {
  excl <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1,  1000000
  )
  incl <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1,  1000000
  )
  expect_error(bed_shuffle(x, genome, incl, excl))
})

test_that('exceeding `max_tries` yields an error', {
  # 100 bp interval is left but x intervals are 1kb
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

test_that('all supplied x interval columns are passed to the result', {
  
    x <- tibble::tribble(
      ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
      "chr1", 80, 100,    "q1",   1,  "+"
    )
    
  res <- bed_shuffle(x, genome)
  expect_true(all(c("strand", "score", "name", "start") %in% colnames(res)))
})



