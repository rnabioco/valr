context("reading BED files")

bed3_path <- system.file('extdata', '3fields.bed.gz', package = 'Rbedtools')
bed6_path <- system.file('extdata', '6fields.bed.gz', package = 'Rbedtools')
bed12_path <- system.file('extdata', 'mm9.bed12.gz', package = 'Rbedtools')
bedgraph_path <- system.file('extdata', 'test.bg.gz', package = 'Rbedtools')

bed3_tbl <- read_bed(bed3_path)

test_that('read BED3', {
  expect_equal(ncol(bed3_tbl), 3)
  expect_equal(nrow(bed3_tbl), 10)
})

bed6_tbl <- read_bed(bed6_path, n_fields = 6)

test_that('read BED6', {
  expect_equal(ncol(bed6_tbl), 6)
  expect_equal(nrow(bed6_tbl), 10)
})

bed12_tbl <- read_bed12(bed12_path)
test_that('read BED12', {
  expect_equal(ncol(bed12_tbl), 12)
})

bedgraph_tbl <- read_bedgraph(bedgraph_path)

test_that('read bedGraph', {
  expect_equal(ncol(bedgraph_tbl), 4)
  expect_equal(nrow(bedgraph_tbl), 4)
})
