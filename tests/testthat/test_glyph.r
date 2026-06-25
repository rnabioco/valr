x <- tibble::tribble(
  ~chrom , ~start , ~end , ~name , ~score , ~strand ,
  "chr1" ,    500 , 1000 , "."   , "."    , "+"     ,
  "chr1" ,   1000 , 1500 , "."   , "."    , "-"
)

test_that("glyphs are rendered", {
  expect_doppelganger("merge glyph is ok", bed_glyph(bed_merge(x)))
})

test_that("glyph labels are applied", {
  res <- bed_glyph(bed_merge(x, id = dplyr::n()), label = "id")
  if ("get_labs" %in% getNamespaceExports("ggplot2")) {
    expect_equal(ggplot2::get_labs(res)$label, "id")
  } else {
    expect_equal(res$labels$label, "id")
  }
})

test_that("an absent label column is an error", {
  expect_snapshot(
    bed_glyph(bed_merge(x), label = "does_not_exist"),
    error = TRUE
  )
})

test_that("namespace-qualified calls are supported", {
  expect_s3_class(bed_glyph(valr::bed_merge(x)), "ggplot")
})

test_that("a genome argument is not faceted as an input", {
  ivl <- tibble::tribble(
    ~chrom , ~start , ~end ,
    "chr1" ,     25 ,   50 ,
    "chr1" ,    100 ,  125
  )
  other <- tibble::tribble(
    ~chrom , ~start , ~end ,
    "chr1" ,     60 ,   75
  )
  g <- tibble::tribble(
    ~chrom , ~size ,
    "chr1" ,   125
  )
  # facets should be ivl, other, result (3) -- genome excluded
  built <- ggplot2::ggplot_build(
    bed_glyph(bed_window(ivl, other, g, both = 15))
  )
  expect_length(built$layout$panel_params, 3)
})

a <- tibble::tribble(
  ~chrom , ~start , ~end ,
  "chr1" ,     25 ,   50 ,
  "chr1" ,    100 ,  125
)

b <- tibble::tribble(
  ~chrom , ~start , ~end , ~value ,
  "chr1" ,     30 ,   75 ,     50
)

test_that("expr arguments do not need to be x and/or y", {
  expect_doppelganger(
    "intersect glyph is ok",
    bed_glyph(bed_intersect(a, b, min_overlap = 0L))
  )
})

# deterministic, non-overlapping intervals so `bed_intersect(z, z)` yields
# exactly nrow(z) rows (each interval overlaps only its own copy)
z <- tibble::tibble(
  chrom = "chr1",
  start = seq(0L, by = 10L, length.out = 101L),
  end = seq(0L, by = 10L, length.out = 101L) + 5L
)

test_that("exceeding max_rows throws an error", {
  expect_snapshot(
    bed_glyph(bed_intersect(z, z, min_overlap = 0L)),
    error = TRUE
  )
})

test_that("max_rows is configurable", {
  small <- head(z, 6)
  expect_error(
    bed_glyph(bed_intersect(small, small, min_overlap = 0L), max_rows = 3)
  )
  expect_s3_class(
    bed_glyph(bed_intersect(small, small, min_overlap = 0L), max_rows = 1000),
    "ggplot"
  )
})
