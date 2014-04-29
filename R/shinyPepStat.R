##' Launch the pepStat Shiny Application
##'
##' Launches the \code{pepStat} Shiny application, providing an interactive
##' interface for constructing peptide sets, normalizing intensities, generating
##' calls. Quality control is also facilitated through interactive plotting
##' features.
##' @export
##' @examples
##' if (interactive()) {
##'   shinyPepStat()
##' }
shinyPepStat <- function() {
  shinyDir <- system.file(package="pepStat", "shiny")
  shiny::runApp(shinyDir)
}
