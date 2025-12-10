# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Package Overview

**valr** is an R package for genome interval arithmetic, providing
BEDtools-like functionality natively in R with full dplyr/tidyverse
integration.

## Build and Test Commands

``` r
# Run tests
devtools::test()

# Run a single test file
devtools::test(filter = "intersect")

# Full R CMD check
devtools::check()

# Generate documentation after modifying roxygen comments
devtools::document()

# Install dev dependencies
devtools::install_dev_deps()
```

## Architecture

### Core Data Structures

- **Interval dataframe (`ivl_df`)**: tibble with required columns
  `chrom`, `start`, `end`
- **Genome dataframe (`genome_df`)**: tibble with required columns
  `chrom`, `size`

### Code Organization

- `R/` - R functions (~34 files)
  - `bed_*.r` - Main interval operations (intersect, merge, closest,
    etc.)
  - `read_*.r` - File I/O (BED, VCF, genome files)
  - `tbls.r` - Data validation (`check_interval`, `check_genome`)
  - `utils.r` - Helpers including
    [`valr_example()`](https://rnabioco.github.io/valr/reference/valr_example.md)
    for example data
- `src/` - C++ implementations via Rcpp (~14 .cpp files)
  - Performance-critical operations use IntervalTree data structure
  - Headers in `inst/include/`
- `tests/testthat/` - Test suite (~33 test files)
  - Uses testthat edition 3 with parallel execution
  - Visual regression tests via vdiffr

### C++ Integration

The package uses Rcpp with an IntervalTree library for efficient
interval lookups. Key headers in `inst/include/`: - `IntervalTree.h` -
Core interval tree implementation - `grouped_dataframe.h` - Handles
grouped operations - `DataFrameBuilder.h` - Efficient dataframe
construction

## Code Style

- Follow tidyverse style guide (<https://style.tidyverse.org>)
- Use `styler::style_pkg()` but don’t restyle unrelated code
- Documentation uses roxygen2 with Markdown syntax
- Code formatter: air (configured in `air.toml` for tribble formatting)

## Key Conventions

- Use `.data[["column"]]` syntax for column references in NSE operations
- Use `all_of()` for column selection in dplyr verbs
- Output columns from operations use `.` prefix (e.g., `.overlap`,
  `.dist`)
- All `bed_*` functions support grouped dataframes via
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)

## Contributing

- Add tests for new functionality
- Update `NEWS.md` for user-facing changes
- Run
  [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  before submitting PRs
