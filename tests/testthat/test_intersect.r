context('bed_intersect')

x <- tibble::tribble(
  ~chrom,   ~start,    ~end,
  "chr1",    100,       200,
  "chr1",    150,       250,
  "chr1",    400,       500
)

y <- tibble::tribble(
  ~chrom,   ~start,    ~end,
  "chr1",    175,       200,
  "chr1",    175,       225
)

test_that("simple overlap works", {
  res <- bed_intersect(x, y)
  expect_equal(nrow(res), 4)   
})

test_that("invert param works", {
  res <- bed_intersect(x, y, invert = TRUE)
  expect_equal(nrow(res), 1)   
})

test_that("multple a's", {
   x <- tibble::tribble(
    ~chrom,    ~start,    ~end,
    "chr1",    100,       200,
    "chr1",    100,       200,
    "chr1",    100,       200,
    "chr1",    100,       200,
    "chr1",    100,       200
  )
  
  y <- tibble::tribble(
    ~chrom,    ~start,    ~end,
    "chr1",    175,       200
  )
  
  res <- bed_intersect(x, y)
  expect_equal(nrow(res), 5)   
})

test_that("multple b's", {
   x <- tibble::tribble(
    ~chrom,    ~start,    ~end,
    "chr1",    100,       200
  )
  
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",    175,       200,
    "chr1",    175,       200,
    "chr1",    175,       200,
    "chr1",    175,       200,
    "chr1",    175,       200
  )
  
  res <- bed_intersect(x, y)
  expect_equal(nrow(res), 5)   
})

test_that("no overlaps returns empty df", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100,    200
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 300,    400
  )
  res <- bed_intersect(x, y)
  expect_is(res, "data.frame")
  expect_equal(nrow(res), 0)
}) 

test_that("duplicate intervals are removed (#23)", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100,    500,
    "chr1", 175,    200
  )
  
  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 150,    400,
    "chr1", 151,    401
  )
  
  res <- bed_intersect(x, y)
  expect_equal(nrow(res), 4)
})

test_that("suffixes disambiguate x/y columns (#28)", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 1000,   1500, '.',   '.',     '-'
  )
  
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 1000,   1200, '.',   '.',     '-'
  )
  
  res <- bed_intersect(x, y)
  test_that("start.y" %in% colnames(res), TRUE)
})

test_that("incorrect `suffix` args throw errors", {
   x <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score,
    "chr1", 1000,   1500, '.',   '.'
  )
  
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score,
    "chr1", 1000,   1200, '.',   '.'
  )
  
  expect_error(bed_intersect(x, y, suffix = 'TESTING')) 
})

test_that("intersections from x bed_tbl with more chroms than y are captured", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",    100,       200,
    "chr3",    400,       500
  )
  
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr3",    425,       475) 
  
  res <- bed_intersect(x, y)
  expect_true("chr3" %in% res$chrom)
})

test_that("intersections from y bed_tbl with more chroms are captured", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr3",    400,       500
  )
  
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",    100,       200,
    "chr3",    425,       475) 
  
  res <- bed_intersect(x, y)
  expect_true("chr3" %in% res$chrom)
})

test_that("input x groups are used for comparing intervals issue #108",{
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    'chr1', 100,    200,  'A',
    'chr1', 200,    400,  'A',
    'chr1', 300,    500,  'A',
    'chr1', 125,    175,  'B',
    'chr1', 150,    200,  'B'
  )
  x <- bed_sort(x)
  x <- group_by(x, group)
  res <- bed_intersect(x, x)
  expect_true(all(res$group.x == res$group.y))
  
})

test_that("tbls grouped by strand are processed", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 1000,   1500, '.',   '.',     '+'
  )
  
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 1000,   1200, '.',   '.',     '-'
  )
 
  res <- bed_intersect(group_by(x, strand), group_by(y, strand))
  expect_equal(nrow(res), 0)
 
  # flip strands 
  res <- bed_intersect(group_by(x, strand), group_by(flip_strands(y), strand))
  expect_equal(nrow(res), 1)
})
