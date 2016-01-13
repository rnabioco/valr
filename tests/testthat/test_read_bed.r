library(Rbedtools)
context("reading BED files")

bed3_tbl <- read_bed('3fields.bed.gz')

test_that('read BED3', {
  expect_equal(ncol(bed3_tbl), 3)
  expect_equal(nrow(bed3_tbl), 10)
})

bed6_tbl <- read_bed('6fields.bed.gz', n_fields = 6)

test_that('read BED6', {
  expect_equal(ncol(bed6_tbl), 6)
  expect_equal(nrow(bed6_tbl), 10)
})

bedgraph_tbl <- read_bedgraph('test.bg.gz')

test_that('read bedGraph', {
  expect_equal(ncol(bedgraph_tbl), 4)
  expect_equal(nrow(bedgraph_tbl), 4)
})

# sorted_tbl <- read_bed('3fields.bed.gz')
# unsorted_tbl <- read_bed('3fields.bed.gz', sort = FALSE)
# 
# test_that('sort attribute is set correctly', {
#   expect_equal(attr(sorted_tbl, 'sorted'), TRUE) 
#   expect_equal(attr(unsorted_tbl, 'sorted'), FALSE) 
# })