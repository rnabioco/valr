test_that("start intervals are NA", {
  x <- tibble::tribble(
    ~chrom , ~start , ~end ,
    "chr1" ,      1 ,  100 ,
    "chr1" ,    150 ,  200 ,
    "chr2" ,    200 ,  300
  )

  res <- interval_spacing(x)

  first <- res |>
    group_by(chrom) |>
    slice(1) |>
    select(chrom:end)

  nas <- filter(res, is.na(.spacing)) |>
    select(chrom:end)

  expect_true(all(first == nas))
})

# from bedtools2
test_that("bt test succeeds", {
  x <- tibble::tribble(
    ~chrom , ~start , ~end ,
    "chr1" ,     20 ,   30 ,
    "chr1" ,     25 ,   40 ,
    "chr1" ,     40 ,   50 ,
    "chr1" ,     60 ,   80 ,
    "chr1" ,     75 ,  100 ,
    "chr1" ,    105 ,  110 ,
    "chr2" ,    115 ,  130 ,
    "chr2" ,    120 ,  160 ,
    "chr2" ,    170 ,  180
  )

  res <- interval_spacing(x)

  exp <- c(NA, -5, 0, 10, -5, 5, NA, -10, 10)
  expect_equal(res$.spacing, exp)
})
