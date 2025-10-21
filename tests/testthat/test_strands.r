x <- tibble::tribble(
  ~chrom , ~start , ~end , ~strand ,
  "chr1" ,      1 ,  100 , "+"     ,
  "chr2" ,      1 ,  100 , "-"
)

test_that("strands are flipped", {
  res <- flip_strands(x)
  expect_equal(res$strand, c("-", "+"))
})

y <- tibble::tribble(
  ~chrom , ~start , ~end ,
  "chr1" ,      1 ,  100 ,
  "chr2" ,      1 ,  100
)

test_that("unstranded tbls throw an error", {
  expect_error(flip_strands(y))
})
