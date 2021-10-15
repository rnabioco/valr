#' Identify intersecting intervals.
#'
#' Report intersecting intervals from `x` and `y` tbls. Book-ended intervals
#' have `.overlap` values of `0` in the output.
#'
#' @param x [ivl_df]
#' @param ... one or more (e.g. a list of) `y` [ivl_df()]s
#' @param invert report `x` intervals not in `y`
#' @param suffix colname suffixes in output
#'
#' @return
#'
#' [ivl_df] with original columns from `x` and `y` suffixed with `.x`
#' and `.y`, and a new `.overlap` column with the extent of overlap for the
#' intersecting intervals.
#'
#' If  multiple `y` tbls are supplied, the `.source` contains variable names
#' associated with each interval. All original columns from the `y` are suffixed
#' with `.y` in the output.
#'
#' If `...` contains named inputs (i.e `a = y, b = z` or `list(a = y, b =  z)`),
#' then `.source` will contain supplied names (see examples).
#'
#' @template groups
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 25,      50,
#'   'chr1', 100,     125
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 30,     75
#' )
#'
#' bed_glyph(bed_intersect(x, y))
#'
#' bed_glyph(bed_intersect(x, y, invert = TRUE))
#'
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 100,    500,
#'   'chr2', 200,    400,
#'   'chr2', 300,    500,
#'   'chr2', 800,    900
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~value,
#'   'chr1', 150,    400,  100,
#'   'chr1', 500,    550,  100,
#'   'chr2', 230,    430,  200,
#'   'chr2', 350,    430,  300
#' )
#'
#' bed_intersect(x, y)
#'
#' bed_intersect(x, y, invert = TRUE)
#'
#' # start and end of each overlapping interval
#' res <- bed_intersect(x, y)
#' dplyr::mutate(res, start = pmax(start.x, start.y),
#'                    end = pmin(end.x, end.y))
#'
#' z <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~value,
#'   'chr1', 150,    400,  100,
#'   'chr1', 500,    550,  100,
#'   'chr2', 230,    430,  200,
#'   'chr2', 750,    900,  400
#' )
#'
#' bed_intersect(x, y, z)
#'
#' bed_intersect(x, exons = y, introns = z)
#'
#' # a list of tbl_intervals can also be passed
#' bed_intersect(x, list(exons = y, introns = z))
#'
#' @family multiple set operations
#'
#' @seealso
#' \url{https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html}
#'
#' @export
bed_intersect <- function(x, ..., invert = FALSE, suffix = c(".x", ".y")) {

  y_tbl <- parse_dots(...)
  multiple_tbls <- FALSE

  if (length(y_tbl) > 1) {
    multiple_tbls <- TRUE
    # bind_rows preserves grouping
    y <- bind_rows(y_tbl, .id = ".source")
    select_vars <- list(quo(-one_of(".source")), quo(everything()), ".source")
    y <- select(y, !!! select_vars)
  } else {
    # only one tbl supplied, so extract out single tbl from list
    y <- y_tbl[[1]]
  }

  x <- check_interval(x)
  y <- check_interval(y)

  check_suffix(suffix)

  # establish grouping with shared groups (and chrom)
  groups_xy <- shared_groups(x, y)
  groups_xy <- unique(as.character(c("chrom", groups_xy)))
  groups_vars <- rlang::syms(groups_xy)

  # type convert grouping factors to characters if necessary and ungroup
  x <- convert_factors(x, groups_xy)
  y <- convert_factors(y, groups_xy)

  x <- group_by(x, !!! groups_vars)
  y <- group_by(y, !!! groups_vars)

  suffix <- list(x = suffix[1], y = suffix[2])

  grp_indexes <- shared_group_indexes(x, y)

  res <- intersect_impl(x, y,
                        grp_indexes$x,
                        grp_indexes$y,
                        invert,
                        suffix$x,
                        suffix$y)

  if (invert) {
    res <- filter(res, is.na(.overlap))
    res <- select(
      res, chrom, ends_with(".x")
    )
    names(res) <- str_replace(names(res), fixed(".x"), "")
    return(res)
  }

  if (multiple_tbls) {
    # rename .source.y to .source
    source_col <- paste0(".source", suffix$y)
    replace_col <- str_replace(
      source_col,
      fixed(suffix$y), ""
    )
    cols <- colnames(res)
    colnames(res) <- ifelse(cols == source_col, replace_col, cols)
  }

  res
}

# handle objects passed to ... in bed_intersect
parse_dots <- function(...){
  # determine if list supplied to ... or series of variables
  n_inputs <- ...length()
  res <- list(...)
  if( typeof(substitute(...)) == "symbol"){
    if(length(res) == 1 & inherits(res[[1]], 'list')){
      # list was passed to ...
      res <- res[[1]]
    } else if (n_inputs > 1 & is.null(names(res))) {
      # multiple objests passed e.g. ... == x, y or a = x, b = y
      # name each tbl based on supplied variable
      dots <- eval(substitute(alist(...)))
      names(res) <- dots
    }
  } else if (n_inputs == 1) {
    # handles list initialized in ...
    # e.g. ...  = list(y, z) or lst[1]
    res <- res[[1]]
    if (is.null(names(res))) {
      # name each tbl based on supplied variable if list initialized in dots
      dots <- eval(substitute(alist(...)))[[1]]
      # extract out variables from language object list(a, b, c)
      dots <- as.character(dots)
      if(dots[1] == "list"){
        dots <- dots[2:length(dots)]
        names(res) <- as.character(dots)
      }
    }
  }
  res
}
