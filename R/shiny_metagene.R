#' Launch interactive metagene session
#'
#' @examples
#' if (interactive()) {
#' shiny_metagene()
#' }
shiny_metagene <- function() {
  app <- system.file("shiny/Imetagene/", package = "Imetagene")
  runApp(app)
}
