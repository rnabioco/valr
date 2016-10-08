#' launch valr shiny demo
#' 
#' @export
launch_shiny <- function() {
  #http://deanattali.com/2015/04/21/r-package-shiny-app/
  appDir <- system.file("shiny", package = "valr")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `valr`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "showcase")
}
