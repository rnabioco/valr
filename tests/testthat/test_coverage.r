context('bed_coverage')

x <- tibble::tribble(
  ~chrom,   ~start,    ~end,
  "chr1",	0,	10,
  "chr1",	15,	20,
  "chr1",	21,	25
)

y <- tibble::tribble(
  ~chrom,   ~start,    ~end,
  "chr1", 3,	15
)

a <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1",	20,	70,	6,	25,	"+",
  "chr1",	50,	100,	1,	25,	"-",
  "chr1",	200,	250,	3,	25,	"+",
  "chr2",	80,	130,	5,	25,	"-",
  "chr2",	150,	200,	4,	25,	"+",
  "chr2",	180,	230,	2,	25,	"-"
)

b <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1",	0,	75,	19,	75,	"+",
  "chr1",	3,	78,	15,	75,	"+",
  "chr1",	74,	149,	7,	75,	"+",
  "chr1",	77,	152,	16,	75,	"-",
  "chr1",	92,	167,	6,	75,	"+",
  "chr1",	116,	191,	9,	75,	"+",
  "chr1",	128,	203,	12,	75,	"-",
  "chr1",	132,	207,	10,	75,	"-",
  "chr1",	143,	218,	20,	75,	"-",
  "chr1",	163,	238,	8,	75,	"-",
  "chr2",	6,	81,	17,	75,	"-",
  "chr2",	39,	114,	4,	75,	"+",
  "chr2",	74,	149,	11,	75,	"-",
  "chr2",	77,	152,	1,	75,	"+",
  "chr2",	114,	189,	3,	75,	"-",
  "chr2",	127,	202,	5,	75,	"-",
  "chr2",	131,	206,	13,	75,	"+",
  "chr2",	137,	212,	2,	75,	"-",
  "chr2",	139,	214,	18,	75,	"-",
  "chr2",	163,	238,	14,	75,	"+"
  )

c <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1",	0, 50,
  "chr1",	12,	20
  )

test_that("default coverage works", {
  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,~.intersections,  ~.basecoverage,  ~.x_length,  ~.fraction,
    "chr1",	20,	70,	6,	25,	"+",	2,	50,	50,	1.0000000,
     "chr1",	50,	100,	1,	25,	"-",	5,	50,	50,	1.0000000,
     "chr1",	200,	250,	3,	25,	"+",	4,	38,	50,	0.7600000,
     "chr2",	80,	130,	5,	25,	"-",	6,	50,	50,	1.0000000,
     "chr2",	150,	200,	4,	25,	"+",	7,	50,	50,	1.0000000,
     "chr2",	180,	230,	2,	25,	"-",	6,	50,	50,	1.0000000
  )
  res <- bed_coverage(a, b)
  expect_true(all(res == pred))   
})

test_that(" strand coverage works (strand = TRUE)", {
  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,~.intersections,  ~.basecoverage,  ~.x_length,  ~.fraction,
    "chr1",	20,	70,	6,	25,	"+",	2,	50,	50,	1.0000000,
    "chr1",	50,	100,	1,	25,	"-",	1,	23,	50,	0.4600000,
    "chr1",	200,	250,	3,	25,	"+",	0,	0,	50,	0.000000,
    "chr2",	80,	130,	5,	25,	"-",	4,	50,	50,	1.0000000,
    "chr2",	150,	200,	4,	25,	"+",	3,	50,	50,	1.0000000,
    "chr2",	180,	230,	2,	25,	"-",	4,	34,	50,	0.6800000
  )
  res <- bed_coverage(a, b, strand = TRUE) %>% bed_sort()
  expect_true(all(res == pred))   
})

test_that(" strand_opp coverage works (strand_opp = TRUE)", {
  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,~.intersections,  ~.basecoverage,  ~.x_length,  ~.fraction,
    "chr1",	20,	70,	6,	25,	"+",	0,	0,	50,	0.0000000,
    "chr1",	50,	100,	1,	25,	"-",	4,	50,	50,	1.0000000,
    "chr1",	200,	250,	3,	25,	"+",	4,	38,	50,	0.760000,
    "chr2",	80,	130,	5,	25,	"-",	2,	50,	50,	1.0000000,
    "chr2",	150,	200,	4,	25,	"+",	4,	50,	50,	1.0000000,
    "chr2",	180,	230,	2,	25,	"-",	2,	50,	50,	1.0000000
  )
  res <- bed_coverage(a, b, strand_opp = TRUE) %>% bed_sort()
  expect_true(all(res == pred))   
})
