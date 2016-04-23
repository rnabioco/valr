context('bed_slop')

genome <- tibble::frame_data(
 ~chrom, ~size,
 "chr1", 5000
)

bed_df <- tibble::frame_data(
 ~chrom, ~start, ~end, ~name, ~score, ~strand,
 "chr1", 500,    1000, '.',   '.',     '+',
 "chr1", 1000,   1500, '.',   '.',     '-'
)

test_that("left arg works", {
  dist <- 100
  out <- bed_df %>% bed_slop(genome, left = dist)
  expect_true(all(bed_df$start - out$start == dist))
})

test_that("right arg works", {
  dist <- 100
  out <- bed_df %>% bed_slop(genome, right = dist)
  expect_true(all(out$end - bed_df$end == dist))
})

test_that("both arg works", {
  dist <- 100
  out <- bed_df %>% bed_slop(genome, both = dist)
  expect_true(all(bed_df$start - out$start == dist))
  expect_true(all(out$end - bed_df$end == dist))
})

# test fraction
# test off chrom
