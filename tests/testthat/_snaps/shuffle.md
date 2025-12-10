# completely excluded intervals throw an error

    Code
      bed_shuffle(x, genome, excl = excl, seed = seed)
    Condition
      Error:
      ! maximum iterations exceeded in bed_shuffle

# empty intervals derived from `incl` and `excl` is handled

    Code
      bed_shuffle(x, genome, incl, excl, seed = seed)
    Condition
      Error in `bed_shuffle()`:
      ! no intervals to sample from.

# exceeding `max_tries` yields an error

    Code
      bed_shuffle(x, genome, excl = excl, seed = seed)
    Condition
      Error:
      ! maximum iterations exceeded in bed_shuffle

