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
                  tags$head(
                    tags$style(
                      HTML("     
                      h1 {
                          font-family: 'Courier New', Courier, monospace;
                      }
                      .navbar-brand {
                          font-family: 'Courier New', Courier, monospace;
                          float:left;
                      }
                    ")
                    ),
                    tags$script('Shiny.addCustomMessageHandler("myCallbackHandler",
                       function(typeMessage) {console.log(typeMessage)
                           if(typeMessage == 1){
                              $("a:contains(DESIGN)").click();
                           } else {
                              if(typeMessage == 2) {
                                  $("a:contains(MATRIX)").click();
                              } else {
                                if(typeMessage == 3) {
                                    $("a:contains(PLOT)").click();
                                }
                              }
                           }
                        });')  
                    
                  ),
                  
                  headerPanel(app_name),
                  sidebarLayout(
                    NULL,
                    mainPanel(
                      p("This application enables to use interactively the metagene package.
        This R package produces Metagene-like plots to compare the behavior of 
        DNA-interacting proteins at selected groups of features."),
                      fluidRow(
                        column(width = 12,
                               navbarPage("Imetagene",
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
                                                                                textInput(inputId = "path_load_metagene", label = ""),
                                                                                bsAlert("load_alert")
                                                                              ), 
                                                                              style = "primary"),
                                                              bsCollapsePanel("Create a new metagene object", 
                                                                              list(
                                                                                fluidRow(
                                                                                  column(width = 6,
                                                                                         list(
                                                                                           shinyFilesButton('bams', 'Select BAM files', 
                                                                                                            'Please select a file', TRUE),
                                                                                           br(),br(),
                                                                                           verbatimTextOutput("bam_list")
                                                                                           #  ,
                                                                                           #  bsButton('clear_bams', 'Clear BAM files', 
                                                                                           #           , style = "btn btn-primary", disabled = FALSE)
                                                                                         )
                                                                                  ),
                                                                                  column(width = 6,
                                                                                         list(
                                                                                           shinyFilesButton('beds', 'Select BED files', 
                                                                                                            'Please select a file', TRUE),
                                                                                           br(),br(),
                                                                                           verbatimTextOutput("bed_list")
                                                                                           #  ,
                                                                                           #  bsButton('clear_beds', 'Clear BED files', 
                                                                                           #           , style = "btn btn-primary", disabled = FALSE)
                                                                                         )
                                                                                  )
                                                                                ),
                                                                                br(),hr(),br(),
                                                                                fluidRow(
                                                                                  column(width = 3,
                                                                                         bsButton(inputId = "runMetagene", label = "Run metagene", style = "btn btn-primary", disabled = FALSE)),
                                                                                  column(width = 8,
                                                                                         downloadButton(outputId = "saveMetagene", label = "Save metagene", class = "btn btn-primary"))
                                                                                ),
                                                                                br(),
                                                                                bsAlert("run_alert"),
                                                                                uiOutput(outputId = "resizingUI")
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
                                                                                fluidRow(
                                                                                  column(width = 4,
                                                                                         list(
                                                                                           shinyFilesButton('loadDesign', 'Choose file', 
                                                                                                            'Please select a file', FALSE),
                                                                                           textInput(inputId = "path_load_design", label = "")
                                                                                         )
                                                                                  ),
                                                                                  column(width = 2,
                                                                                         list(
                                                                                           checkboxInput(inputId = "header", label = "File with header", value = TRUE)
                                                                                         )
                                                                                  )
                                                                                ),
                                                                                dataTableOutput('loaded_design'),
                                                                                bsAlert("load_alert_d")
                                                                              ), 
                                                                              style = "primary"),
                                                              bsCollapsePanel("Create a new design", 
                                                                              list(
                                                                                textInput(inputId = "exp_name", label = "Experience name", value = "Exp1"),
                                                                                fluidRow(
                                                                                  column(width = 6,
                                                                                         list(
                                                                                           selectInput(inputId = "chips", 
                                                                                                       label = "Select the CHiP for your experiment",
                                                                                                       selected = NULL,
                                                                                                       multiple = TRUE,
                                                                                                       width = '100%',
                                                                                                       selectize = TRUE,
                                                                                                       choices = c(''))
                                                                                         )
                                                                                  ),
                                                                                  column(width = 6,
                                                                                         list(
                                                                                           selectInput(inputId = "ctrls", 
                                                                                                       label = "Select the associated controls for your experiment",
                                                                                                       selected = NULL,
                                                                                                       multiple = TRUE,
                                                                                                       width = '100%',
                                                                                                       selectize = TRUE,
                                                                                                       choices = c(''))
                                                                                         )
                                                                                  )
                                                                                ),
                                                                                br(),br(),
                                                                                bsAlert("save_exp_alert"),
                                                                                bsButton(inputId = "saveExp", label = "Save Experiment", style = "btn btn-primary", disabled = FALSE),
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
                                                   bsCollapse(id = "matrix", 
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
                                                                                checkboxInput(inputId = "design", value = TRUE, label = "use design"),
                                                                                br(),br(),
                                                                                fluidRow(
                                                                                  column(width = 3,
                                                                                         bsButton(inputId = "runMatrix", label = "Produce matrix", style = "btn btn-primary", disabled = TRUE)),
                                                                                  column(width = 8,
                                                                                         bsButton(inputId = "updateMetagene", label = "Update metagene", style = "btn btn-primary", disabled = TRUE))
                                                                                )
                                                                              )
                                                                              ,style = "primary"
                                                              )
                                                   )
                                          ),tabPanel("PLOT",
                                                     helpText("PUT SOME DESCRIPTION HERE"),
                                                     bsCollapse(id = "plot", 
                                                                open = c("Plot parameters"), 
                                                                multiple = FALSE,
                                                                bsCollapsePanel("Plot parameters", 
                                                                                list(
                                                                                  fluidRow(
                                                                                    column(width = 6,
                                                                                           list(
                                                                                             selectInput(inputId = "plot_regions", 
                                                                                                         label = "Select the regions you want to draw",
                                                                                                         selected = NULL,
                                                                                                         multiple = TRUE,
                                                                                                         width = '100%',
                                                                                                         selectize = TRUE,
                                                                                                         choices = c("REGION1","REGION2","REGION3"))
                                                                                           )
                                                                                    ),
                                                                                    column(width = 6,
                                                                                           list(
                                                                                             textInput(inputId = "plot_title", label = "plot title", value = "Enter text..."),
                                                                                             sliderInput(inputId = "alpha", label = "alpha", min = 0, max = 1, value = 0),
                                                                                             numericInput(inputId = "sample_count", label = "sample count",value = 0)
                                                                                           )
                                                                                    )
                                                                                  ),
                                                                                  br(),br(),
                                                                                  bsButton(inputId = "runPlot", label = "Plot", style = "btn btn-primary", disabled = FALSE),
                                                                                  br(),br(),
                                                                                  plotOutput(outputId = "plot"),
                                                                                  fluidRow(
                                                                                    column(width = 6,
                                                                                           bsButton(inputId = "savePlotPNG", label = "Export PNG", style = "btn btn-primary", disabled = TRUE)),
                                                                                    column(width = 6,
                                                                                           bsButton(inputId = "savePlotPDF", label = "Export PDF", style = "btn btn-primary", disabled = TRUE))
                                                                                  )
                                                                                )
                                                                                ,style = "primary"
                                                                )
                                                     )
                                          )
                               )
                        )
                      ), width="100%",footer()
                    )
                  )
))
