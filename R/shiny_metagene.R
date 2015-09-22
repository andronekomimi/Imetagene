#' Launch interactive metagene session
#'
#' @examples
#' shiny_metagene
shiny_metagene <- function() {
    app <- system.file("shiny/Imetagene/", package = "Imetagene")
    runApp(app)
}
