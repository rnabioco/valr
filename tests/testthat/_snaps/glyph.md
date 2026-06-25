# an absent label column is an error

    Code
      bed_glyph(bed_merge(x), label = "does_not_exist")
    Condition
      Error in `glyph_plot()`:
      ! `label` ("does_not_exist") is not a column in the result.

# exceeding max_rows throws an error

    Code
      bed_glyph(bed_intersect(z, z, min_overlap = 0L))
    Condition
      Error in `bed_glyph()`:
      ! Result has 101 rows, exceeding `max_rows` (100). Use a smaller example or raise `max_rows`.

