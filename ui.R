# IMETAGENE UI.R
library(shiny)
library(shinyBS)
library(shinyFiles)
library(shinythemes)

footer<-function(){
  tags$div(
    class = "footer",
    style = "text-align:center",
    tags$div(class = "foot-inner",
             list(
               hr(),
               "Imetagene is an ", 
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


customNumInputs<-function(id1, id2, before = "", after = "", value1 = "", value2 = ""){
  div(style="display:inline-block;font-weight:normal",
      tags$label(before, `for` = "none"), 
      tags$input(id = id1, type = "number", value = value1, class="input-small"),
      tags$label(after, `for` = "none"),
      tags$input(id = id2, type = "number", value = value2, class="input-small")
  )
}


app_name = "Imetagene"

shinyUI(fluidPage(theme = shinytheme("flatly"),
  #   tags$head(
  #     tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css"),
  #     tags$link(rel = "icon", type = "image/x-icon", href = "logo.png")
  #   ),
  
  headerPanel(app_name),
  
  sidebarLayout(
    NULL,
    mainPanel(
      p("This application enables to use interactively the metagene package.
        This R package produces Metagene-like plots to compare the behavior of 
        DNA-interacting proteins at selected groups of features."),
      fluidRow(
        column(width = 12,
               navbarPage("I-metagene",
                          tabPanel("INPUTS",
                                   helpText("There is no hard limit in the number of BAM files that can be included in an analysis 
                                             (but with too many BAM files, memory may become an issue). BAM files must be indexed.
                                              For instance, if you use a file names file.bam, a file named file.bam.bai must be present in the same directory."),
                                   helpText("To compare custom regions of interest, it is possible to use a list of one or more BED files."),
                                   bsCollapse(id = "inputs", 
                                              open = c("Load existing metagene object"), 
                                              multiple = FALSE,
                                              bsCollapsePanel("Load existing metagene object", 
                                                              list(
                                                                shinyFilesButton('loadMetagene', 'Choose file', 
                                                                                 'Please select a file', FALSE),
                                                                textInput(inputId = "path_load_metagene", label = "")
                                                              ), 
                                                              style = "primary"),
                                              bsCollapsePanel("Create a new metagene object", 
                                                              list(
                                                                fluidRow(
                                                                  column(width = 6,
                                                                         list(
                                                                           shinyFilesButton('bams', 'Select BAM files', 
                                                                                            'Please select a file', FALSE),
                                                                           br(),br(),
                                                                           verbatimTextOutput("bam_list")
                                                                         )
                                                                  ),
                                                                  column(width = 6,
                                                                         list(
                                                                           shinyFilesButton('regions', 'Select BED files', 
                                                                                            'Please select a file', FALSE),
                                                                           br(),br(),
                                                                           verbatimTextOutput("bed_list")
                                                                         )
                                                                  )
                                                                ),
                                                                br(),br(),
                                                                fluidRow(
                                                                  column(width = 4,
                                                                         bsButton(inputId = "runMetagene", label = "Run metagene", style = "btn btn-primary", disabled = TRUE)),
                                                                  column(width = 8,
                                                                         bsButton(inputId = "saveMetagene", label = "Save metagene", style = "btn btn-primary", disabled = TRUE))
                                                                )
                                                              ), 
                                                              style = "primary")
                                   ),
                                   bsButton(inputId = "go2design", label = "Next step", style = "btn btn-primary", disabled = TRUE, icon = icon("arrow-right"))
                          ),
                          tabPanel("DESIGN",
                                   helpText("A design group contains a set of BAM files that, when pull together, 
                                            represent a logical analysis. Furthermore, a design group contains 
                                            the relationship between every BAM files present."),
                                   helpText("Samples (with or without replicates) and controls can be assigned to a same design group."),
                                   helpText("There can be as many groups as necessary. A BAM file can be assigned to more than one group."),
                                   bsCollapse(id = "design", 
                                              open = c("Load existing design file"), 
                                              multiple = FALSE,
                                              bsCollapsePanel("Load existing design file", 
                                                              list(
                                                                helpText("ADD HELP BUBBLE"),
                                                                shinyFilesButton('loadDesign', 'Choose file', 
                                                                                 'Please select a file', FALSE),
                                                                textInput(inputId = "path_load_design", label = "")
                                                              ), 
                                                              style = "primary"),
                                              bsCollapsePanel("Create a new design", 
                                                              list(
                                                                textInput(inputId = "exp_name", label = "Experience name", value = "Enter text..."),
                                                                fluidRow(
                                                                  column(width = 6,
                                                                         list(
                                                                           h3("CHIP")
                                                                         )
                                                                  ),
                                                                  column(width = 6,
                                                                         list(
                                                                           h3("CONTROLS")
                                                                         )
                                                                  )
                                                                ),
                                                                br(),br(),
                                                                bsButton(inputId = "saveExp", label = "Save Experiment", style = "btn btn-primary", disabled = TRUE),
                                                                br(),br(),
                                                                dataTableOutput('design'),
                                                                bsButton(inputId = "saveDesign", label = "Save Design", style = "btn btn-primary", disabled = TRUE)
                                                              ), 
                                                              style = "primary")
                                   ),
                                   bsButton(inputId = "go2matrix", label = "Next step", style = "btn btn-primary", disabled = TRUE, icon = icon("arrow-right"))
                          ),
                          tabPanel("MATRIX",
                                   helpText("PUT SOME DESCRIPTION HERE"),
                                   bsCollapse(id = "MATRIX", 
                                              open = c("Matrix parameters"), 
                                              multiple = FALSE,
                                              bsCollapsePanel("Matrix parameters", 
                                                              list(
                                                                customNumInputs(id1="bin_count", id2 = "bin_size",
                                                                                before = "bin count ", 
                                                                                after = "or bin size", 
                                                                                value1 = 10,
                                                                                value2 = 100),
                                                                br(),br(),
                                                                selectInput(inputId = "noise", choices = c("NOISE1","NOISE2","NOISE3"), label = "noise removal"),
                                                                selectInput(inputId = "norm", choices = c("NORM1","NORM2","NORM3"), label = "normalization"),
                                                                checkboxInput(inputId = "flip", value = FALSE, label = "flip regions"),
                                                                checkboxInput(inputId = "design", value = TRUE, label = "use design")
                                                              )
                                                              ,style = "primary"
                                              )
                                   )
                          ),tabPanel("PLOT",
                                     h3("test4")
                          )
               )
        )
      ), width="100%",footer()
    )
  )
))
