context("spacing")

x <- tibble::tribble(
~chrom, ~start, ~end,
"chr1", 1,      100,
"chr1", 150,    200,
"chr2", 200,    300
)

test_that("start intervals are NA", {
  res <- interval_spacing(x)
  first <- res %>% group_by(chrom) %>% slice(1) %>% select(chrom:end)
  nas <- filter(res, is.na(.spacing)) %>% select(chrom:end)
  expect_true(all(first == nas))
})
