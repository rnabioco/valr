context("bed_makewindows")

genome <- tibble::frame_data(
  ~chrom, ~size,
  "chr1", 5000,
  "chr2", 400
)

x <- tibble::frame_data(
  ~chrom, ~start, ~end, ~name,
  "chr1", 100, 200, 'A',
  "chr2", 300, 350, 'B'
) 

# Fixed number of windows 
res <- bed_makewindows(x, genome, num_win = 10)
#test_that('num_win are correct', {
#  expect_that(nrow(res), 20)
#})

# Fixed window size
res <- bed_makewindows(x, genome, win_size = 10)

# Fixed window size with overlaps
res <- bed_makewindows(x, genome, win_size = 10, step_size = 5)

# Reversed window numbering
res <- bed_makewindows(x, genome, win_size = 10, reverse = TRUE)

test_that('window IDs are generated', {
  res <- bed_makewindows(x, genome, win_size = 10)
  expect_true('win_id' %in% colnames(res))
})
