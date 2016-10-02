#https://github.com/arq5x/bedtools2/blob/master/test/closest/test-closest.sh
context("bed_closest")


test_that("1bp closer, check for off-by-one errors", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  ) 
  
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10,
    "chr1", 19, 20,
    "chr1",	20, 21
  )
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 3)   
})

test_that("reciprocal test of 1bp closer, check for off-by-one errors", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10,
    "chr1", 19, 20,
    "chr1",	20, 21
  ) 
  res <- bed_closest(y, x)
  expect_equal(nrow(res), 3)   
})

test_that("0bp apart closer, check for off-by-one errors", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10,
    "chr1", 19, 21,
    "chr1",	20, 21
  ) 
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 3)   
})

test_that("reciprocal of 0bp apart closer, check for off-by-one errors", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10,
    "chr1", 19, 21,
    "chr1",	20, 21
  ) 
  res <- bed_closest(y, x)
  expect_equal(nrow(res), 3)   
})

test_that("check that first left interval at index 0 is not lost", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10
  ) 
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 1) 
}
)

test_that("check that first right interval at index 0 is not lost", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1", 20, 21
  ) 
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 1) 
}
)

test_that("check that strand closest works (strand = TRUE)", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 100, 200, "a", 10,	"+")
  
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 90, 120, "b", 1,	"-")
  
  res <- bed_closest(x, y, strand = TRUE)
  expect_equal(nrow(res), 0)
}
)

test_that("check that same strand is reported (strand = TRUE", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1",	80,	100,	"q1",	1,	"+"
  )
    
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1",	5,	  15,	"d1.1",	1,	"+",
    "chr1",	20,	  60,	"d1.2",	2,	"-",
    "chr1",	200,	220,	"d1.3",	3,	"-"
  ) 
  
  pred <- tibble::tribble(
    ~chrom,   ~start.x,    ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y,    ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",	 80,	100,	"q1",	1,	"+",	5, 	15, 	"d1.1",	1,	"+", 0, -65 
  )
  
  bed_closest(x, y, strand = T)
  res <- bed_closest(x, y, strand = T)
  expect_true(all(pred == res))
}
)

test_that("check that different strand is reported (strand_opp = TRUE", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1",	80,	100,	"q1",	1,	"+"
  )
  
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1",	5,	  15,	"d1.1",	1,	"+",
    "chr1",	20,	  60,	"d1.2",	2,	"-",
    "chr1",	200,	220,	"d1.3",	3,	"-"
  ) 
  
  pred <- tibble::tribble(
    ~chrom,   ~start.x,    ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y,    ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",	 80,	100,	"q1",	1,	"+",	20, 	60, 	"d1.2",	2,	"-", 0, -20 
  )
  
  res <- bed_closest(x, y, strand_opp = T)
  expect_true(all(pred == res))
}
)

test_that("check that reciprocal strand closest works (strand_opp = TRUE) ", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 100, 200, "a", 10,	"+")
  
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 80, 90, "b", 1,	"-")
  
  res <- bed_closest(x, y, strand_opp = TRUE)
  expect_equal(nrow(res), 1)
}
)

test_that("check that stranded distance reporting works ( distance_type = 'strand') ", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 100, 200, "a", 10,	"+",
    "chr1", 100, 200, "a", 10,	"-"
    )
  
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 310, 320, "b", 1,	"-")
  
  res <- bed_closest(x, y, distance_type = "strand")
  expect_true(res$.distance[1] > 0 &&  res$.distance[2] < 0)
}
)

test_that("check that abs distance reporting works (distance_type = 'abs')", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 100, 200, "a", 10,	"+",
    "chr1", 350, 400, "a", 10,	"+"
  )
  
  y <- tibble::tribble(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 310, 320, "b", 1,	"+")
  
  res <- bed_closest(x, y, distance_type = "abs")
  expect_true(res$.distance[1] > 0 &&  res$.distance[2] > 0)
}
)

test_that("overlapping intervals are removed (overlap = F)", {
  x <- tibble::tribble(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  
  y <- tibble::tribble(
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
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100,    200
    )
  
  y <- tibble::tribble(
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
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100,    200
  )
  
  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100,    200,
    "chr1", 150,    200,
    "chr1", 50,    100,
    "chr1", 200,   300
  )
  exp <- tibble::tribble(
    ~chrom, ~start.x, ~start.y,
    "chr1", 100,    200
  )
  res <- bed_closest(x, y)
  expect_true(nrow(res) == 4)
}
)

test_that("test reporting of first overlapping feature and 
           overlap = F excludes overlapping intervals", {
  
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1",	100,	101,
    "chr1",	200,	201,
    "chr1",	300,	301,
    "chr1",	100000,	100010,
    "chr1",	100020,	100040,
    "chr2",	1,	10,
    "chr2",	20,	30)

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1",	100,	101,
    "chr1",	150,	201,
    "chr1",	175,	375
)
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~start.y, ~end.y, ~.distance,
    "chr1",	100,	101,  150,  201,  49,
    "chr1",	200,	201,  100,  101, -99,
    "chr1",	300,	301,  150,  201, -99,
    "chr1",	100000,  100010,	175,  375, -99625,
    "chr1",	100020,  100040,  175,  375, -99645
)
  res <- bed_closest(x, y, overlap = F)
  expect_true(all(pred == res))
}
)

### test all distance reporting conditions ###

### tbls to test
d_q1 <- tibble::tribble(
  ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
  "chr1",	80,	100,	"d_q1.1",	5,	"+"
)

d_q2 <- tibble::tribble(
  ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
  "chr1",	80,	100,	"d_q2.1",	5,	"-"
)

d_d1F <- tibble::tribble(
  ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
  "chr1",	40,	60,	"d1F.1",	10,	"+"
)

d_d1R <- tibble::tribble(
  ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
  "chr1",	40,	60,	"d1R.1",	10,	"-"
)

d_d2F <- tibble::tribble(
  ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
  "chr1",	140,	160,	"d2F.1",	10,	"+"
)

d_d2R <- tibble::tribble(
  ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
  "chr1",	140,	160,	"d2R.1",	10,	"-"
)

test_that("default distance reporting works for forward hit on left, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q1.1",       5,        "+",      40,    60,  "d1F.1",      10,        "+",        0,       -20
  ) 
  res <- bed_closest(d_q1, d_d1F)
  expect_true(all(pred == res))
}
)

test_that("strand distance reporting works for forward hit on left, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q1.1",       5,        "+",      40,    60,  "d1F.1",      10,        "+",        0,       -20
  ) 
  res <- bed_closest(d_q1, d_d1F, distance_type = "strand")
  expect_true(all(pred == res))
}
)

test_that("abs distance reporting works for forward hit on left, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q1.1",       5,        "+",      40,    60,  "d1F.1",      10,        "+",        0,       20
  ) 
  res <- bed_closest(d_q1, d_d1F, distance_type = "abs")
  expect_true(all(pred == res))
}
)

test_that("default distance reporting works for reverse hit on left, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q1.1",       5,        "+",      40,    60,  "d1R.1",      10,        "-",        0,       -20
  ) 
  res <- bed_closest(d_q1, d_d1R)
  expect_true(all(pred == res))
}
)

test_that("strand distance reporting works for reverse hit on left, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q1.1",       5,        "+",      40,    60,  "d1R.1",      10,        "-",        0,       -20
  ) 
  res <- bed_closest(d_q1, d_d1R, distance_type = "strand")
  expect_true(all(pred == res))
}
)

test_that("default distance reporting works for forward hit on left, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q2.1",       5,        "-",      40,    60,  "d1F.1",      10,        "+",        0,       -20
  ) 
  res <- bed_closest(d_q2, d_d1F)
  expect_true(all(pred == res))
}
)

test_that("strand distance reporting works for forward hit on left, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q2.1",       5,        "-",      40,    60,  "d1F.1",      10,        "+",        0,       20
  ) 
  res <- bed_closest(d_q2, d_d1F, distance_type = "strand")
  expect_true(all(pred == res))
}
)

test_that("abs distance reporting works for forward hit on left, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q2.1",       5,        "-",      40,    60,  "d1F.1",      10,        "+",        0,       20
  ) 
  res <- bed_closest(d_q2, d_d1F, distance_type = "abs")
  expect_true(all(pred == res))
}
)

test_that("default distance reporting works for reverse hit on left, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q2.1",       5,        "-",      40,    60,  "d1R.1",      10,        "-",        0,       -20
  ) 
  res <- bed_closest(d_q2, d_d1R)
  expect_true(all(pred == res))
}
)

test_that("strand distance reporting works for reverse hit on left, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q2.1",       5,        "-",      40,    60,  "d1R.1",      10,        "-",        0,       20
  ) 
  res <- bed_closest(d_q2, d_d1R, distance_type = "strand")
  expect_true(all(pred == res))
}
)

test_that("abs distance reporting works for reverse hit on left, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q2.1",       5,        "-",      40,    60,  "d1R.1",      10,        "-",        0,       20
  ) 
  res <- bed_closest(d_q2, d_d1R, distance_type = "abs")
  expect_true(all(pred == res))
}
)


test_that("default distance reporting works for forward hit on right, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q1.1",       5,        "+",      140,    160,  "d2F.1",      10,        "+",        0,      40 
  ) 
  res <- bed_closest(d_q1, d_d2F)
  expect_true(all(pred == res))
}
)

test_that("strand distance reporting works for forward hit on right, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q1.1",       5,        "+",      140,    160,  "d2F.1",      10,        "+",        0,       40
  ) 
  res <- bed_closest(d_q1, d_d2F, distance_type = "strand")
  expect_true(all(pred == res))
}
)

test_that("abs distance reporting works for forward hit on right, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q1.1",       5,        "+",      140,    160,  "d2F.1",      10,        "+",        0,       40
  ) 
  res <- bed_closest(d_q1, d_d2F, distance_type = "abs")
  expect_true(all(pred == res))
}
)

test_that("default distance reporting works for reverse hit on right, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q1.1",       5,        "+",      140,    160,  "d2R.1",      10,        "-",        0,      40 
  ) 
  res <- bed_closest(d_q1, d_d2R)
  expect_true(all(pred == res))
}
)

test_that("strand distance reporting works for rverse hit on right, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q1.1",       5,        "+",      140,    160,  "d2R.1",      10,        "-",        0,       40
  ) 
  res <- bed_closest(d_q1, d_d2R, distance_type = "strand")
  expect_true(all(pred == res))
}
)

test_that("abs distance reporting works for reverse hit on right, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q1.1",       5,        "+",      140,    160,  "d2R.1",      10,        "-",        0,       40
  ) 
  res <- bed_closest(d_q1, d_d2R, distance_type = "abs")
  expect_true(all(pred == res))
}
)

test_that("default distance reporting works for forward hit on right, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q2.1",       5,        "-",      140,    160,  "d2F.1",      10,        "+",        0,      40 
  ) 
  res <- bed_closest(d_q2, d_d2F)
  expect_true(all(pred == res))
}
)

test_that("strand distance reporting works for forward hit on right, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q2.1",       5,        "-",      140,    160,  "d2F.1",      10,        "+",        0,       -40
  ) 
  res <- bed_closest(d_q2, d_d2F, distance_type = "strand")
  expect_true(all(pred == res))
}
)

test_that("abs distance reporting works for forward hit on right, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q2.1",       5,        "-",      140,    160,  "d2F.1",      10,        "+",        0,       40
  ) 
  res <- bed_closest(d_q2, d_d2F, distance_type = "abs")
  expect_true(all(pred == res))
}
)

test_that("default distance reporting works for reverse hit on right, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q2.1",       5,        "-",      140,    160,  "d2R.1",      10,        "-",        0,      40 
  ) 
  res <- bed_closest(d_q2, d_d2R)
  expect_true(all(pred == res))
}
)

test_that("strand distance reporting works for reverse hit on right, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q2.1",       5,        "-",      140,    160,  "d2R.1",      10,        "-",        0,       -40
  ) 
  res <- bed_closest(d_q2, d_d2R, distance_type = "strand")
  expect_true(all(pred == res))
}
)

test_that("abs distance reporting works for reverse hit on right, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      80,   100, "d_q2.1",       5,        "-",      140,    160,  "d2R.1",      10,        "-",        0,       40
  ) 
  res <- bed_closest(d_q2, d_d2R, distance_type = "abs")
  expect_true(all(pred == res))
}
)

### additional tbls for tests ###
a2 <- tibble::tribble(
  ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
  "chr1",	10,	20,	"a1",	1,	"-"
)

b2 <- tibble::tribble(
  ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
  "chr1",	8,	9,	"b1",	1,	"+",
  "chr1",	21,	22,	"b2",	1, "-"
)

test_that("Make sure non-overlapping ties are reported ", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      10,   20, "a1",       1,        "-",      21,    22,  "b2",      1,        "-",        0,       1,
    "chr1",      10,   20, "a1",       1,        "-",      8,    9,  "b1",      1,        "+",        0,       -1
  ) 
  res <- bed_closest(a2, b2)
  expect_true(all(pred == res))
}
)

test_that("Make sure non-overlapping ties are reported with strand = T ", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      10,   20, "a1",       1,        "-",      21,    22,  "b2",      1,        "-",        0,       1
  ) 
  res <- bed_closest(a2, b2, strand = T)
  expect_true(all(pred == res))
}
)


test_that("Make sure non-overlapping ties are reported with strand_opp = T ", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1",      10,   20, "a1",       1,        "-",      8,    9,  "b1",      1,        "+",        0,       -1
  ) 
  res <- bed_closest(a2, b2, strand_opp = T)
  expect_true(all(pred == res))
}
)










