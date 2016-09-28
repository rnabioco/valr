context('bed_complement')

test_that("complement with covering interval", {
 
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 500
  ) 
  bed <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",    1,       500
  )
  res <- bed_complement(bed, genome)
  expect_equal(nrow(res), 0)   
})

test_that("complement with middle interval", {
 
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 500
  ) 
  bed <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",    100,       200
  )
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1",      1,   100,
    "chr1",    200,   500
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("complement adds final interval", {
  
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 500
  ) 
  bed <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",    1,        50,
    "chr1",    100,      200
  )
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1",     50,   100,
    "chr1",    200,   500
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("multiple chroms", {
  
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 500,
    "chr2", 500
  ) 
  bed <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",      100,     500,
    "chr2",      100,     500
  )
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1",      1,   100,
    "chr2",      1,   100
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("multiple chroms, chr1 is covered", {
  
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 500,
    "chr2", 500
  ) 
  bed <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",        1,     500,
    "chr2",      100,     500
  )
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr2",      1,   100
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("record exceeds chrom length", {
  
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 100
  ) 
  bed <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",       90,     110
  )
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1",      1,   90
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("test complement only on right side", {  ## This test fails ##
  
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 100
  ) 
  bed <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",        1,      20
  )
  res <- bed_complement(bed, genome)
  expect_equal(nrow(res), 1)
})

test_that("test complement only on left side", {
  
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 100
  ) 
  bed <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",        90,      100
  )
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1",      1,   90
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("nothing is covered", {  ## This test fails ##
  
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 100,
    "chr2", 100
  ) 
  bed <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",       10,     100
  )
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1",      1,   10,
    "chr2",      1,  100
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("chrom is not in genome file", {  ## This test fails ##
  
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 100
  ) 
  bed <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr2",       10,     100   
  )
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1",      1,  100
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("overlapping intervals", {  
  
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1",   100
  ) 
  bed <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",       10,     100,
    "chr1",       20,      80
  )
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1",      1,  10
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})

test_that("bedtools issue #356", {  
  
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 249250621
  ) 
  bed <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",       1,    10000
  )
  expect <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1",  10000, 249250621
  )
  res <- bed_complement(bed, genome)
  expect_true(all(res == expect))
})



