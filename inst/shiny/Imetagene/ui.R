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

app_name = "Imetagene"

shinyUI(fluidPage(theme = shinythemes::shinytheme("flatly"),
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
                                                   shinyBS::bsCollapse(id = "inputs", 
                                                                       open = c("Load existing metagene object"), 
                                                                       multiple = FALSE,
                                                                       shinyBS::bsCollapsePanel("Load existing metagene object", 
                                                                                                list(
                                                                                                  helpText("Load a metagene object produced during a previous session."),
                                                                                                  shinyFiles::shinyFilesButton('loadMetagene', 'Choose file',
                                                                                                                               'Please select a file', FALSE),
                                                                                                  shiny::actionLink(inputId = "hloadMetagene", label = "", icon = icon("question-circle")),
                                                                                                  textInput(inputId = "path_load_metagene", label = "", width = "100%"),
                                                                                                  shinyBS::bsAlert("load_alert")
                                                                                                ), 
                                                                                                style = "primary"),
                                                                       shinyBS::bsCollapsePanel("Create a new metagene object", 
                                                                                                list(
                                                                                                  helpText("Create a metagene object from scratch with BAM and BED files."),
                                                                                                  fluidRow(
                                                                                                    column(width = 6,
                                                                                                           list(
                                                                                                             shinyFiles::shinyFilesButton('bams', 'Select BAM files', 
                                                                                                                                          'Please select a file', TRUE),
                                                                                                             shiny::actionLink(inputId = "hbams", label = "", icon = icon("question-circle")),
                                                                                                             br(),
                                                                                                             uiOutput("bam_list")
                                                                                                           )
                                                                                                    ),
                                                                                                    column(width = 6,
                                                                                                           list(
                                                                                                             shinyFiles::shinyFilesButton('beds', 'Select BED files', 
                                                                                                                                          'Please select a file', TRUE),
                                                                                                             shiny::actionLink(inputId = "hbeds", label = "", icon = icon("question-circle")),
                                                                                                             br(),
                                                                                                             uiOutput("bed_list")
                                                                                                           )
                                                                                                    )
                                                                                                  ),
                                                                                                  br(),hr(),br(),
                                                                                                  fluidRow(
                                                                                                    column(width = 3,
                                                                                                           shinyBS::bsButton(inputId = "runMetagene", label = "Run metagene", style = "btn btn-primary", disabled = FALSE)),
                                                                                                    column(width = 8,
                                                                                                           downloadButton(outputId = "saveMetagene", label = "Save metagene", class = "btn btn-primary"))
                                                                                                  ),
                                                                                                  br(),
                                                                                                  shinyBS::bsAlert("run_alert"),
                                                                                                  uiOutput(outputId = "resizingUI")
                                                                                                ), 
                                                                                                style = "primary")
                                                   ),
                                                   shinyBS::bsButton(inputId = "go2design", label = "Next step", style = "btn btn-primary", disabled = TRUE, icon = icon("arrow-right"))
                                          ),
                                          tabPanel("DESIGN",
                                                   shinyBS::bsCollapse(id = "design", 
                                                                       open = c("Current design"), 
                                                                       multiple = TRUE,
                                                                       shinyBS::bsCollapsePanel("Current design", 
                                                                                                list(
                                                                                                  fluidRow(
                                                                                                    column(width = 12,
                                                                                                           list(
                                                                                                             helpText("Here is the design contained in the current metagene if any."),
                                                                                                             dataTableOutput('current_mg_design')
                                                                                                           )
                                                                                                    )
                                                                                                  )
                                                                                                ), 
                                                                                                style = "primary"),
                                                                       shinyBS::bsCollapsePanel("Load existing design file", 
                                                                                                list(
                                                                                                  fluidRow(
                                                                                                    column(width = 8,
                                                                                                           list(
                                                                                                             shinyFiles::shinyFilesButton('loadDesign', 'Choose file', 
                                                                                                                                          'Please select a file', FALSE),
                                                                                                             shiny::actionLink(inputId = "hloadDesign", label = "", icon = icon("question-circle")),
                                                                                                             textInput(inputId = "path_load_design", label = "", width = "100%")
                                                                                                           )
                                                                                                    ),
                                                                                                    column(width = 3,
                                                                                                           list(
                                                                                                             checkboxInput(inputId = "header", label = "File with header", value = TRUE)
                                                                                                           )
                                                                                                    )
                                                                                                  ),
                                                                                                  shinyBS::bsAlert("load_alert_d")
                                                                                                ), 
                                                                                                style = "primary"),
                                                                       shinyBS::bsCollapsePanel("Create a new design", 
                                                                                                list(
                                                                                                  textInput(inputId = "exp_name", label = "Experience name", value = "Exp1"),
                                                                                                  fluidRow(
                                                                                                    column(width = 6,
                                                                                                           list(
                                                                                                             selectInput(inputId = "chips", 
                                                                                                                         label = "Select the ChIP for your experiment",
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
                                                                                                  shinyBS::bsAlert("save_exp_alert"),
                                                                                                  shinyBS::bsButton(inputId = "saveExp", label = "Save Experiment", style = "btn btn-primary", disabled = FALSE),
                                                                                                  br(),br(),
                                                                                                  downloadButton(outputId = "saveDesign", label = "Save Current Design")
                                                                                                ), 
                                                                                                style = "primary")
                                                   ),
                                                   shinyBS::bsButton(inputId = "go2matrix", label = "Next step", style = "btn btn-primary", disabled = FALSE, icon = icon("arrow-right"))
                                          ),
                                          tabPanel("MATRIX",
                                                   helpText("In order to represent the genomic coverages in a metagene plot, the coverages must be converted in matrix format. Regions are binned to reduce computation time for the following steps. During this step, the noise can be removed using one or more controls and the coverages can be normalized to allow comparison of multiple experiments."),
                                                   shinyBS::bsCollapse(id = "matrix", 
                                                                       open = c("Current matrix subset"), 
                                                                       multiple = TRUE,
                                                                       shinyBS::bsCollapsePanel("Current matrix subset",
                                                                                                list(
                                                                                                  fluidRow(
                                                                                                    column(width = 12,
                                                                                                           list(
                                                                                                             helpText("Here is a subset of the matrix contained in the current metagene if any."),
                                                                                                             verbatimTextOutput('current_mg_matrix_content'),
                                                                                                             uiOutput('current_mg_matrix_heatmap')
                                                                                                           )
                                                                                                    )
                                                                                                  )
                                                                                                ), 
                                                                                                style = "primary"),
                                                                       shinyBS::bsCollapsePanel("Matrix parameters", 
                                                                                                list(
                                                                                                  fluidRow(
                                                                                                    column(width = 3,
                                                                                                           numericInput(inputId = "bin_size", value = 100, label = "bin size")),
                                                                                                    column(width = 8,
                                                                                                           shiny::actionLink(inputId = "hbin_size", label = "", icon = icon("question-circle")))
                                                                                                  ),
                                                                                                  fluidRow(
                                                                                                    column(width = 3,
                                                                                                           selectInput(inputId = "noise", choices = c("NONE","NCIS"), label = "noise removal")),
                                                                                                    column(width = 8,
                                                                                                           shiny::actionLink(inputId = "hnoise", label = "", icon = icon("question-circle")))
                                                                                                  ),
                                                                                                  fluidRow(
                                                                                                    column(width = 3,
                                                                                                           selectInput(inputId = "norm", choices = c("NONE","RPM"), label = "normalization")),
                                                                                                    column(width = 8,
                                                                                                           shiny::actionLink(inputId = "hnorm", label = "", icon = icon("question-circle")))
                                                                                                  ),
                                                                                                  fluidRow(
                                                                                                    column(width = 3,
                                                                                                           checkboxInput(inputId = "flip", value = FALSE, label = "flip regions")),
                                                                                                    column(width = 8,
                                                                                                           shiny::actionLink(inputId = "hflip", label = "", icon = icon("question-circle")))
                                                                                                  ),
                                                                                                  fluidRow(
                                                                                                    column(width = 3,
                                                                                                           checkboxInput(inputId = "use_design", value = TRUE, label = "use design")),
                                                                                                    column(width = 8,
                                                                                                           shiny::actionLink(inputId = "huse_design", label = "", icon = icon("question-circle")))
                                                                                                  ),
                                                                                                  br(),
                                                                                                  fluidRow(
                                                                                                    column(width = 3,
                                                                                                           shinyBS::bsButton(inputId = "runMatrix", label = "Produce matrix", style = "btn btn-primary", disabled = FALSE)),
                                                                                                    column(width = 8,
                                                                                                           downloadButton(outputId = "updateMetagene", label = "Save metagene", class = "btn btn-primary"))
                                                                                                  ),
                                                                                                  br(),
                                                                                                  shinyBS::bsAlert("run_matrix_alert")
                                                                                                )
                                                                                                ,style = "primary"
                                                                       )
                                                   ),
                                                   shinyBS::bsButton(inputId = "go2plot", label = "Next step", style = "btn btn-primary", disabled = FALSE, icon = icon("arrow-right"))
                                          ),tabPanel("PLOT",
                                                     helpText("Produce a metagene plot."),
                                                     shinyBS::bsCollapse(id = "plot", 
                                                                         open = c("Plot parameters"), 
                                                                         multiple = FALSE,
                                                                         shinyBS::bsCollapsePanel("Plot parameters", 
                                                                                                  list(
                                                                                                    fluidRow(
                                                                                                      column(width = 4,
                                                                                                             list(
                                                                                                               selectInput(inputId = "plot_regions", 
                                                                                                                           label = "Regions to display",
                                                                                                                           selected = NULL,
                                                                                                                           multiple = TRUE,
                                                                                                                           width = '75%',
                                                                                                                           selectize = TRUE,
                                                                                                                           choices = c("")),
                                                                                                               selectInput(inputId = "plot_exps", 
                                                                                                                           label = "Experiments to display",
                                                                                                                           selected = NULL,
                                                                                                                           multiple = TRUE,
                                                                                                                           width = '75%',
                                                                                                                           selectize = TRUE,
                                                                                                                           choices = c("")),
                                                                                                               textInput(inputId = "plot_title", label = "plot title", value = "Enter text..."),
                                                                                                               fluidRow(
                                                                                                                 column(width = 9,
                                                                                                                        sliderInput(inputId = "alpha", label = "alpha", min = 0, max = 1, value = 0.05)),
                                                                                                                 column(width = 3,
                                                                                                                        shiny::actionLink(inputId = "halpha", label = "", icon = icon("question-circle")))
                                                                                                               ),
                                                                                                               fluidRow(
                                                                                                                 column(width = 9,
                                                                                                                        numericInput(inputId = "sample_count", label = "sample count",value = 1000)),
                                                                                                                 column(width = 3,
                                                                                                                        shiny::actionLink(inputId = "hsample_count", label = "", icon = icon("question-circle")))
                                                                                                               )
                                                                                                             )
                                                                                                      ),
                                                                                                      column(width = 8,
                                                                                                             list(
                                                                                                               plotOutput(outputId = "mg_plot")
                                                                                                             )
                                                                                                      )
                                                                                                    ),
                                                                                                    br(),
                                                                                                    shinyBS::bsButton(inputId = "runPlot", label = "Plot", style = "btn btn-primary", disabled = FALSE),
                                                                                                    br(),br(),
                                                                                                    fluidRow(
                                                                                                      column(width = 6,
                                                                                                             downloadButton(outputId = "savePlotPNG", label = "Export PNG")),
                                                                                                      column(width = 6,
                                                                                                             downloadButton(outputId = "savePlotPDF", label = "Export PDF"))
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
