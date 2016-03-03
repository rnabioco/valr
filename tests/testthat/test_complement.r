context('bed_complement')

test_that("complement with covering interval", {
 
  genome <- tibble(
    ~chrom, ~size,
    "chr1", 500
  ) 
  bed <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    1,       500
  )
  
  res <- bed_complement(bed, genome)
  expect_equal(nrow(res), 0)   
})

test_that("complement with middle interval", {
 
  genome <- tibble(
    ~chrom, ~size,
    "chr1", 500
  ) 
  bed <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    100,       200
  )
  
  res <- bed_complement(bed, genome)
  expect_equal(nrow(res), 2)   
})

test_that("complement adds final interval", {
  
  genome <- tibble(
    ~chrom, ~size,
    "chr1", 500
  ) 
  bed <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    1,        50,
    "chr1",    100,      200
  )
  
  res <- bed_complement(bed, genome)
  expect_equal(nrow(res), 2)   
})

