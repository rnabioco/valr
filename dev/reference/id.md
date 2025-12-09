# Plyr function id packaged due to plyr being retired Compute a unique numeric id for each unique row in a data frame.

Properties:

- `order(id)` is equivalent to `do.call(order, df)`

- rows containing the same data have the same value

- if `drop = FALSE` then room for all possibilites

## Usage

``` r
id(.variables, drop = FALSE)
```

## Arguments

- .variables:

  list of variables

- drop:

  drop unusued factor levels?

## Value

a numeric vector with attribute n, giving total number of possibilities

## See also

[`id_var`](https://rnabioco.github.io/valr/dev/reference/id_var.md)
