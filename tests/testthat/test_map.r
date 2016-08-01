context('bed_map')

test_that("`chrom` grouping throws an error", {
  x <- tibble::frame_data(
    ~chrom, ~start, ~end,
    "chr1", 100, 250,
    "chr2", 250, 500
  ) %>% group_by(chrom)
  
  y <- tibble::frame_data(
    ~chrom, ~start, ~end, ~value,
    "chr1", 100, 250, 10,
    "chr1", 150, 250, 20,
    "chr2", 250, 500, 500
  ) %>% group_by(chrom)
  
  expect_error(bed_map(x, y))
})

test_that("x/y groupings are respected", {
  x <- tibble::frame_data(
    ~chrom, ~start, ~end, ~id,
    "chr1", 100, 250, 1,
    "chr2", 250, 500, 2,
    "chr2", 250, 500, 3
  ) %>% group_by(id)
  
  y <- tibble::frame_data(
    ~chrom, ~start, ~end, ~value, ~id,
    "chr1", 100, 250, 10,  1,
    "chr1", 150, 250, 20,  2,
    "chr2", 250, 500, 500, 3
  ) %>% group_by(id)
  
  res <- bed_map(x, y, vals = sum(value.y))
  expect_equal(res$vals, c(10,20,500,500))
})
