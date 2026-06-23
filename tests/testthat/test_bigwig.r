test_that("is_bigwig_path() recognizes bigWig/bigBed paths and URLs", {
  expect_true(is_bigwig_path("signal.bw"))
  expect_true(is_bigwig_path("signal.bigWig"))
  expect_true(is_bigwig_path("anno.bb"))
  expect_true(is_bigwig_path("anno.bigBed"))
  expect_true(is_bigwig_path("https://host.org/data/signal.bw?token=abc"))

  expect_false(is_bigwig_path("intervals.bed"))
  expect_false(is_bigwig_path("signal.bw.gz"))
  expect_false(is_bigwig_path(c("a.bw", "b.bw")))
  expect_false(is_bigwig_path(NA_character_))
  expect_false(is_bigwig_path(42))
})

test_that("read_bigwig_regions() errors on an unsupported extension", {
  x <- tibble::tibble(chrom = "chr1", start = 1L, end = 10L)
  expect_error(read_bigwig_regions(x, "data.txt"), "bigWig or bigBed")
})

test_that("bed_map() accepts a bigWig file path for `y`", {
  skip_if_not_installed("cpp11bigwig")

  bw <- valr_example("hg19.dnase1.bw")
  whole <- read_bigwig(bw)
  x <- head(dplyr::select(whole, "chrom", "start", "end"), 10)

  from_file <- bed_map(x, bw, .sum = sum(value), .n = length(value))
  from_tbl <- bed_map(x, whole, .sum = sum(value), .n = length(value))

  expect_identical(from_file, from_tbl)
})

test_that("bed_intersect() accepts a bigBed file path for `y`", {
  skip_if_not_installed("cpp11bigwig")

  bb <- valr_example("test.bb")
  whole <- read_bigbed(bb)
  x <- dplyr::mutate(
    dplyr::select(whole, "chrom", "start", "end"),
    start = .data[["start"]] - 10L
  )

  from_file <- bed_intersect(x, bb, min_overlap = 1L)
  from_tbl <- bed_intersect(x, whole, min_overlap = 1L)

  expect_identical(from_file, from_tbl)
})

test_that("bed_subtract() accepts a bigBed file path for `y`", {
  skip_if_not_installed("cpp11bigwig")

  bb <- valr_example("test.bb")
  whole <- read_bigbed(bb)
  x <- dplyr::mutate(
    dplyr::select(whole, "chrom", "start", "end"),
    start = .data[["start"]] - 10L
  )

  from_file <- bed_subtract(x, bb, min_overlap = 1L)
  from_tbl <- bed_subtract(x, whole, min_overlap = 1L)

  expect_identical(from_file, from_tbl)
})

test_that("bed_coverage() accepts a bigBed file path for `y`", {
  skip_if_not_installed("cpp11bigwig")

  bb <- valr_example("test.bb")
  whole <- read_bigbed(bb)
  x <- dplyr::mutate(
    dplyr::select(whole, "chrom", "start", "end"),
    start = .data[["start"]] - 10L
  )

  from_file <- bed_coverage(x, bb, min_overlap = 1L)
  from_tbl <- bed_coverage(x, whole, min_overlap = 1L)

  expect_identical(from_file, from_tbl)
})

test_that("bed_window() accepts a bigBed file path for `y` (scoped to the window)", {
  skip_if_not_installed("cpp11bigwig")

  bb <- valr_example("test.bb")
  whole <- read_bigbed(bb)
  genome <- tibble::tibble(
    chrom = c("chr1", "chr10", "chr20"),
    size = 6000000L
  )
  # x sits just outside each feature; only a window picks the feature up
  x <- dplyr::mutate(
    dplyr::select(whole, "chrom", "start", "end"),
    start = .data[["end"]] + 100L,
    end = .data[["end"]] + 200L
  )

  from_file <- bed_window(x, bb, genome, both = 1000)
  from_tbl <- bed_window(x, whole, genome, both = 1000)

  expect_identical(from_file, from_tbl)
})

test_that("overlapping query intervals do not double-count file entries", {
  skip_if_not_installed("cpp11bigwig")

  bw <- valr_example("hg19.dnase1.bw")
  whole <- read_bigwig(bw)

  # two identical (fully overlapping) query intervals spanning real signal
  rng <- range(whole$start, whole$end)
  x <- tibble::tibble(
    chrom = "chr22",
    start = rep(rng[1], 2),
    end = rep(rng[2], 2)
  )

  y <- read_bigwig_regions(x, bw)
  expect_identical(y, dplyr::distinct(y))
})
