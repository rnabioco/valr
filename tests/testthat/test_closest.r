#https://github.com/arq5x/bedtools2/blob/master/test/closest/test-closest.sh
context("bed_closest")


test_that("1bp closer, check for off-by-one errors", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  ) 
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10,
    "chr1", 19, 20,
    "chr1",	20, 21
  )
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 3)   
})

test_that("reciprocal test of 1bp closer, check for off-by-one errors", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10,
    "chr1", 19, 20,
    "chr1",	20, 21
  ) 
  res <- bed_closest(y, x)
  expect_equal(nrow(res), 3)   
})

test_that("0bp apart closer, check for off-by-one errors", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10,
    "chr1", 19, 21,
    "chr1",	20, 21
  ) 
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 3)   
})

test_that("reciprocal of 0bp apart closer, check for off-by-one errors", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10,
    "chr1", 19, 21,
    "chr1",	20, 21
  ) 
  res <- bed_closest(y, x)
  expect_equal(nrow(res), 3)   
})

test_that("check that first left interval at index 0 is not lost", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10
  ) 
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 1) 
}
)

test_that("check that first right interval at index 0 is not lost", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1", 20, 21
  ) 
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 1) 
}
)

test_that("check that strand closest works (strand = TRUE)", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 100, 200, "a", 10,	"+")
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 90, 120, "b", 1,	"-")
  
  res <- bed_closest(x, y, strand = TRUE)
  expect_equal(nrow(res), 0)
}
)

test_that("check that reciprocal strand closest works (strand_opp = TRUE) ", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 100, 200, "a", 10,	"+")
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 80, 90, "b", 1,	"-")
  
  res <- bed_closest(x, y, strand_opp = TRUE)
  expect_equal(nrow(res), 1)
}
)

test_that("check that stranded distance reporting works ( distance_type = 'strand') ", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 100, 200, "a", 10,	"+",
    "chr1", 100, 200, "a", 10,	"-"
    )
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 310, 320, "b", 1,	"-")
  
  res <- bed_closest(x, y, distance_type = "strand")
  expect_true(res$.distance[1] > 0 &&  res$.distance[2] < 0)
}
)

test_that("check that abs distance reporting works (distance_type = 'abs')", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 100, 200, "a", 10,	"+",
    "chr1", 350, 400, "a", 10,	"+"
  )
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 310, 320, "b", 1,	"+")
  
  res <- bed_closest(x, y, distance_type = "abs")
  expect_true(res$.distance[1] > 0 &&  res$.distance[2] > 0)
}
)

test_that("overlapping intervals are removed (overlap = F)", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10,
    "chr1", 19, 21,
    "chr1",	20, 21
  ) 
  
  res <- bed_closest(x, y, overlap = F)
  expect_true(res[2, "start.y"] != 19)
}
)

test_that("duplicate intervals are not reported", {
  x <- tibble::frame_data(
    ~chrom, ~start, ~end,
    "chr1", 100,    200
    )
  
  y <- tibble::frame_data(
    ~chrom, ~start, ~end,
    "chr1", 100,    200,
    "chr1", 150,    200,
    "chr1", 550,    580,
    "chr2", 7000,   8500
    )
  res <- bed_closest(x, y)
  expect_false(any(duplicated(res)))
}
)

test_that("all overlapping features are reported", {
  x <- tibble::frame_data(
    ~chrom, ~start, ~end,
    "chr1", 100,    200
  )
  
  y <- tibble::frame_data(
    ~chrom, ~start, ~end,
    "chr1", 100,    200,
    "chr1", 150,    200,
    "chr1", 50,    100,
    "chr1", 200,   300
  )
  exp <- tibble::frame_data(
    ~chrom, ~start.x, ~start.y,
    "chr1", 100,    200
  )
  res <- bed_closest(x, y)
  expect_true(nrow(res) == 4)
}
)


