context('bed_window')

genome <- tibble::frame_data(
  ~chrom, ~size,
  "chr1", 5000,
  "chr2", 10000
)

bed_df_x <- tibble::frame_data(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 500,    1000, '.',   '.',     '+',
  "chr1", 1000,   1500, '.',   '.',     '-',
  "chr2", 1000,   1200, '.',   '.',     '-'
)

bed_df_y <- tibble::frame_data(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 400,    450, '.',   '.',     '+',
  "chr1", 1000,   1200, '.',   '.',     '-',
  "chr1", 1100,    1500, '.',   '.',     '+',
  "chr2", 1300,   1500, '.',   '.',     '-'
)

test_that("both arg works", {
  out <- bed_window(bed_df_x, bed_df_y, genome, both = 110)
  expect_equal(nrow(out), 8)
})

test_that("left arg works", {
  out <- bed_window(bed_df_x, bed_df_y, genome, left = 110)
  expect_equal(nrow(out), 5)
})

test_that("right arg works", {
  out <- bed_window(bed_df_x, bed_df_y, genome, right = 110)
  expect_equal(nrow(out), 7)
})

test_that("strand position arg works", {
  out <- bed_window(bed_df_x, bed_df_y, genome, right = 110, 
                    strand_pos = TRUE)
  expect_equal(nrow(out), 4)
})

test_that("strand intersect arg works", {
  dist <- 100
  out <- bed_window(bed_df_x, bed_df_y, genome, right = 110, 
                    strand_pos = TRUE, strand = TRUE)
  expect_equal(nrow(out), 1)
})

# test fraction
