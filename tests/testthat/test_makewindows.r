context("bed_makewindows")

x <- tibble::tribble(
  ~chrom, ~start, ~end, ~name,
  "chr1", 100, 200, "A",
  "chr2", 300, 350, "B"
)

# Window IDs are generated
test_that("window IDs are generated", {
  res <- bed_makewindows(x, win_size = 10)
  expect_true(".win_id" %in% colnames(res))
})

# Fixed win_size with foward numbering
test_that("win_size fwd", {
  res <- bed_makewindows(x, win_size = 10)
  # test number of windows
  expect_equal(nrow(res), 15)
  # test forward window numbering
  expect_true(all(res$.win_id == c(1:10, 1:5)))
  # test interval size
  expect_true(all(res$end - res$start == 10))
})

# Fixed win_size with reverse numbering
test_that("win_size rev", {
  res <- bed_makewindows(x, reverse = TRUE, win_size = 10)
  # test number of windows
  expect_equal(nrow(res), 15)
  # test forward window numbering
  expect_true(all(res$.win_id == c(10:1, 5:1)))
  # test interval size
  expect_true(all(res$end - res$start == 10))
})

# Fixed win_size +step_size with forward numbering
test_that("win_size +step_size fwd", {
  res <- bed_makewindows(x, win_size = 10, step_size = 5)
  # test number of windows
  expect_equal(nrow(res), 30)
  # test forward window numbering
  expect_true(all(res$.win_id == c(1:20, 1:10)))
  # test interval size
  expect_true(all(res[1:19, "end"] - res[1:19, "start"] == 10 & res$end[20] - res$start[20] == 5))
  expect_true(all(res[21:29, "end"] - res[21:29, "start"] == 10 & res$end[30] - res$start[30] == 5))
})

# Fixed win_size +step_size with reverse numbering
test_that("win_size +step_size rev", {
  res <- bed_makewindows(x, reverse = TRUE, win_size = 10, step_size = 5)
  # test number of windows
  expect_equal(nrow(res), 30)
  # test forward window numbering
  expect_true(all(res$.win_id == c(20:1, 10:1)))
  # test interval size
  expect_true(all(res[1:19, "end"] - res[1:19, "start"] == 10 & res$end[20] - res$start[20] == 5))
  expect_true(all(res[21:29, "end"] - res[21:29, "start"] == 10 & res$end[30] - res$start[30] == 5))
})

# Fixed number of windows with forward numbering
test_that("num_win fwd", {
  res <- bed_makewindows(x, num_win = 10)
  # test number of windows
  expect_equal(nrow(res), 20)
  # test forward window numbering
  expect_true(all(res$.win_id == c(1:10, 1:10)))
  # test interval size
  expect_true(all(res[1:10, "end"] - res[1:10, "start"] == 10))
  expect_true(all(res[11:20, "end"] - res[11:20, "start"] == 5))
})

# Fixed number of windows with reverse numbering
test_that("num_win rev", {
  res <- bed_makewindows(x, reverse = TRUE, num_win = 10)
  # test number of windows
  expect_equal(nrow(res), 20)
  # test forward window numbering
  expect_true(all(res$.win_id == c(10:1, 10:1)))
  # test interval size
  expect_true(all(res[1:10, "end"] - res[1:10, "start"] == 10))
  expect_true(all(res[11:20, "end"] - res[11:20, "start"] == 5))
})

test_that("interval is smaller than n windows", {
  # test warning
  expect_warning(
    bed_makewindows(x, num_win = 150),
    "Interval [^:]+:\\d+-\\d+, smaller than requested number of windows. skipping"
  )
  # test that intervals are dropped if num_win > than interval size
  res <- suppressWarnings(bed_makewindows(x, num_win = 150))
  expect_equal(nrow(res), 0)
})

# from https://github.com/arq5x/bedtools2/blob/master/test/makewindows/test-makewindows.sh
test_that("always get the number of requested windows. issue #322", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~name,
    "chr1", 11, 44, "A"
  )
  res <- bed_makewindows(x, num_win = 10)
  expect_equal(nrow(res), 10)
  expect_equal(res$start[10], 38)
  expect_equal(res$end[10], 44)
})

test_that("report error if negative value win_size or num_win arguments supplied", {
  expect_error(bed_makewindows(x, num_win = -1))
  expect_error(bed_makewindows(x, win_size = -1))
})

test_that("check num_win reported correctly for additional intervals (related to #322)", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 9437053, 9438070,
    "chr1", 75360291, 75368579,
    "chr1", 86351980, 86352127
  )
  res <- bed_makewindows(x, num_win = 100)
  expect_equal(max(res$.win_id), 100)
})
