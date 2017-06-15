context("bed_makewindows")

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 5000,
  "chr2", 400
)

x <- tibble::tribble(
  ~chrom, ~start, ~end, ~name,
  "chr1", 100, 200, 'A',
  "chr2", 300, 350, 'B'
)

# Window IDs are generated
test_that('window IDs are generated', {
  res <- bed_makewindows(x, genome, win_size = 10)
  expect_true('.win_id' %in% colnames(res))
})

# Fixed win_size with foward numbering
test_that('win_size fwd', {
  res <- bed_makewindows(x, genome, win_size = 10)
  # test number of windows
  expect_equal(nrow(res), 15)
  # test forward window numbering
  expect_true(all(res$.win_id == c(1:10, 1:5)))
  # test interval size
  expect_true(all(res$end - res$start == 10))
})

# Fixed win_size with reverse numbering
test_that('win_size rev', {
  res <- bed_makewindows(x, genome, reverse = TRUE, win_size = 10)
  # test number of windows
  expect_equal(nrow(res), 15)
  # test forward window numbering
  expect_true(all(res$.win_id == c(10:1, 5:1)))
  # test interval size
  expect_true(all(res$end - res$start == 10))
})

# Fixed win_size +step_size with forward numbering
test_that('win_size +step_size fwd', {
  res <- bed_makewindows(x, genome, win_size = 10, step_size = 5)
  # test number of windows
  expect_equal(nrow(res), 30)
  # test forward window numbering
  expect_true(all(res$.win_id == c(1:20, 1:10)))
  # test interval size
  expect_true(all(res[1:19, 'end'] - res[1:19, 'start'] == 10 & res$end[20] - res$start[20] == 5))
  expect_true(all(res[21:29, 'end'] - res[21:29, 'start'] == 10 & res$end[30] - res$start[30] == 5))
})

# Fixed win_size +step_size with reverse numbering
test_that('win_size +step_size rev', {
  res <- bed_makewindows(x, genome, reverse = TRUE, win_size = 10, step_size = 5)
  # test number of windows
  expect_equal(nrow(res), 30)
  # test forward window numbering
  expect_true(all(res$.win_id == c(20:1, 10:1)))
  # test interval size
  expect_true(all(res[1:19, 'end'] - res[1:19, 'start'] == 10 & res$end[20] - res$start[20] == 5))
  expect_true(all(res[21:29, 'end'] - res[21:29, 'start'] == 10 & res$end[30] - res$start[30] == 5))
})

# Fixed number of windows with forward numbering
test_that('num_win fwd', {
  res <- bed_makewindows(x, genome, num_win = 10)
  # test number of windows
  expect_equal(nrow(res), 20)
  # test forward window numbering
  expect_true(all(res$.win_id == c(1:10, 1:10)))
  # test interval size
  expect_true(all(res[1:10, 'end'] - res[1:10, 'start'] == 10))
  expect_true(all(res[11:20, 'end'] - res[11:20, 'start'] == 5))
})

# Fixed number of windows with reverse numbering
test_that('num_win rev', {
  res <- bed_makewindows(x, genome, reverse = TRUE, num_win = 10)
  # test number of windows
  expect_equal(nrow(res), 20)
  # test forward window numbering
  expect_true(all(res$.win_id == c(10:1, 10:1)))
  # test interval size
  expect_true(all(res[1:10, 'end'] - res[1:10, 'start'] == 10))
  expect_true(all(res[11:20, 'end'] - res[11:20, 'start'] == 5))
})



