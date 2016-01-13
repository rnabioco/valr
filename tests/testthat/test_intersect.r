library(Rbedtools)
library(dplyr)
context('intersect_cpp')

test_that("simple overlap", {
   bed1_df <- tibble(
    ~start,    ~end,
    100,       200,
    150,       250
  )
  
  bed2_df <- tibble(
    ~start,    ~end,
    175,       200,
    175,       225
  )
  
  res <- intersect_cpp(bed1_df, bed2_df)
  expect_equal(nrow(res), 4)   
})

test_that("multple a's", {
   bed1_df <- tibble(
    ~start,    ~end,
    100,       200,
    100,       200,
    100,       200,
    100,       200,
    100,       200
  )
  
  bed2_df <- tibble(
    ~start,    ~end,
    175,       200
  )
  
  res <- intersect_cpp(bed1_df, bed2_df)
  expect_equal(nrow(res), 5)   
})

  
test_that("multple b's", {
   bed1_df <- tibble(
    ~start,    ~end,
    100,       200
  )
  
  bed2_df <- tibble(
    ~start,    ~end,
    175,       200,
    175,       200,
    175,       200,
    175,       200,
    175,       200
  )
  
  res <- intersect_cpp(bed1_df, bed2_df)
  expect_equal(nrow(res), 5)   
})

  