context('utils')

x <- trbl_interval(
  ~chrom, ~start, ~end,
  "chr1", 1     , 100,
  "chr1", 200,    500
)
y <- trbl_interval(
  ~chrom, ~start, ~end,
  "chr1", 1     , 100,
  "chr1", 200,    500
) %>% group_by(chrom)

test_that("NULL is return when there are no shared groups, shared_groups()", {
  res <- shared_groups(x, y)
  expect_null(res)
})

test_that("only shared groups are return, shared_groups()", {
  x <- trbl_interval(
    ~chrom, ~start, ~end,
    "chr1", 1     , 100,
    "chr1", 200,    500
  ) %>% group_by(chrom, start, end)

  res <- shared_groups(x, y)
  expect_true(length(res) == 1 && res == "chrom")
})


x <- trbl_interval(
  ~end,  ~chrom,   ~start, ~value,
  75,  "chr1",    125,    10
)

y <- trbl_interval(
  ~chrom,   ~start,    ~end,  ~scores,
  "chr1",    50,       100,  1.2,
  "chr1",    100,       150,  2.4
)

test_that("x columns are reordered based on y, reorder_names()", {
  res <- reorder_names(x, y)

  pred <- intersect(colnames(y), colnames(res))
  n <- length(pred)

  expect_true(all(colnames(res[1:n]) == pred[1:n]))
})
