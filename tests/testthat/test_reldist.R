context('bed_reldist')

x <- tibble::tribble(
  ~chrom,   ~start,    ~end,
  "chr1",    75,       125
)

y <- tibble::tribble(
  ~chrom,   ~start,    ~end,
  "chr1",    50,       100,
  "chr1",    100,       150
)

test_that("reldist calculation is correct", {
  res <- bed_reldist(x, y)
  expect_true(res$reldist == 0.5)   
})

test_that("self reldist is 0", {
  res <- bed_reldist(y, y)
  expect_true(res$reldist == 0)   
})

test_that("detail argument works", {
  res <- bed_reldist(x, y, detail = T)
  expect_true(all(names(res) %in% c("chrom", "start", "end", "reldist")))
})


