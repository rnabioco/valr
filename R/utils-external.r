#' Plyr function id packaged due to plyr being retired
#' Compute a unique numeric id for each unique row in a data frame.
#'
#' Properties:
#' \itemize{
#'   \item \code{order(id)} is equivalent to \code{do.call(order, df)}
#'   \item rows containing the same data have the same value
#'   \item if \code{drop = FALSE} then room for all possibilites
#' }
#'
#' @param .variables list of variables
#' @param drop drop unusued factor levels?
#' @return a numeric vector with attribute n, giving total number of
#'   possibilities
#' @keywords internal
#' @seealso \code{\link{id_var}}
id <- function(.variables, drop = FALSE) {
  # Drop all zero length inputs
  lengths <- vapply(.variables, length, integer(1))
  .variables <- .variables[lengths != 0]

  if (length(.variables) == 0) {
    # inlined "%||%" function to avoid packaging
    nvars <- nrow(.variables)
    n <- if (is.null(nvars)) {
      0L
    } else {
      nvars
    }
    return(structure(seq_len(n), n = n))
  }

  # Special case for single variable
  if (length(.variables) == 1) {
    return(id_var(.variables[[1]], drop = drop))
  }

  # Calculate individual ids
  ids <- rev(lapply(.variables, id_var, drop = drop))
  p <- length(ids)

  # Calculate dimensions
  ndistinct <- vapply(ids, attr, "n", FUN.VALUE = numeric(1), USE.NAMES = FALSE)

  n <- prod(ndistinct)
  if (n > 2^31) {
    # Too big for integers, have to use strings, which will be much slower :(

    char_id <- do.call("paste", c(ids, sep = "\r"))
    res <- match(char_id, unique(char_id))
  } else {
    combs <- c(1, cumprod(ndistinct[-p]))

    mat <- do.call("cbind", ids)
    res <- c((mat - 1L) %*% combs + 1L) # nolint
  }
  attr(res, "n") <- n

  if (drop) {
    id_var(res, drop = TRUE)
  } else {
    structure(as.integer(res), n = attr(res, "n"))
  }
}

#' Plyr function id_var packaged due to plyr being retired
#' Numeric id for a vector.
#' @keywords internal
id_var <- function(x, drop = FALSE) {
  if (length(x) == 0) {
    return(structure(integer(), n = 0L))
  }
  if (!is.null(attr(x, "n")) && !drop) {
    return(x)
  }

  if (is.factor(x) && !drop) {
    x <- addNA(x, ifany = TRUE)
    id <- as.integer(x)
    n <- length(levels(x))
  } else {
    levels <- sort(unique(x), na.last = TRUE)
    id <- match(x, levels)
    n <- max(id)
  }
  structure(id, n = n)
}
