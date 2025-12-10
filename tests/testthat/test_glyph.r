x <- tibble::tribble(
  ~chrom , ~start , ~end , ~name , ~score , ~strand ,
  "chr1" ,    500 , 1000 , "."   , "."    , "+"     ,
  "chr1" ,   1000 , 1500 , "."   , "."    , "-"
)

test_that("glyphs are rendered", {
  skip_on_ci()
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
  skip_on_ci()
  expect_doppelganger(
    "intersect glyph is ok",
    bed_glyph(bed_intersect(a, b, min_overlap = 0L))
  )
})

genome <- tibble::tribble(
  ~chrom , ~size ,
  "chr1" ,   1e6
)

x <- bed_random(genome, n = 101)
test_that("exceeding max intervals throws an error", {
  expect_snapshot(bed_glyph(bed_intersect(x, x, min_overlap = 0L)), error = TRUE)
})
