# IMETAGENE UI.R
library(shiny)
library(shinyBS)
library(shinyFiles)

app_name = "Imetagene"

footer<-function(){
  tags$div(
    class = "footer",
    style = "text-align:center",
    tags$div(class = "foot-inner",
             list(
               hr(),
               "Imetage is an ", 
               tags$a(href="http://bioinformatique.genome.ulaval.ca/recherche/","Arnaud Droit Lab"),
               " project (2015).",
               br(),
               "Development of metagene library : ",
               tags$a(href="https://ca.linkedin.com/pub/charles-joly-beauparlant/89/491/3b3", "Charles Joly-Beauparlant "),
               br(),
               "Development of Imetagene web interface ",
               tags$a(href="https://ca.linkedin.com/in/audreylemacon", "Audrey Lemacon ")
             )
             
    )
  )
}

shinyUI(fluidPage(
  #   tags$head(
  #     tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css"),
  #     tags$link(rel = "icon", type = "image/x-icon", href = "logo.png")
  #   ),
  
  headerPanel(app_name),
  
  sidebarLayout(
    NULL,
    mainPanel(
      fluidRow(
        column(width = 12,
               navbarPage("",
                          tabPanel("INPUTS",
                                   bsCollapse(id = "load", 
                                              open = c("Load existing metagene object"), 
                                              multiple = TRUE,
                                              bsCollapsePanel("Load existing metagene object", 
                                                              list(
                                                                shinyFilesButton('files', 'Choose file', 
                                                                                 'Please select a file', FALSE),
                                                                textInput(inputId = "path", label = "")
                                                              ), 
                                                              style = "primary"),
                                              bsCollapsePanel("Create a new metagene object", 
                                                              list(h3("test1.2")), 
                                                              style = "primary")
                                   ),
                                   bsButton(inputId = "runMetagene", label = "Run metagene", style = "btn btn-primary", disabled = TRUE)
                          ),
                          tabPanel("DESIGN",
                                   h3("test2")
                          ),
                          tabPanel("MATRIX",
                                   h3("test3")
                          ),tabPanel("PLOT",
                                     h3("test4")
                          )
               )
        )
      ), width="100%",footer()
    )
  )
))
