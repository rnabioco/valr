#' @export
run_example <- function() {
  appDir <- system.file("shiny", package = "valr")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "showcase")
}
