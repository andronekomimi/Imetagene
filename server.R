# IMETAGENE SERVER.R
library(metagene)
source("helper.R")
roots <- c(wd='.')
bams <- NULL
beds <- NULL
sizes <- 0
metagene_object <- NULL

shinyServer(function(input, output, session) {
  
  #### INPUTS ####
  ### LOAD EXISTING METAGENE
  observe({
    shinyFileChoose(input, 'loadMetagene', session=session, roots=roots,
                    filetypes=c('Rda'))
    pfile = parseFilePaths(roots, input$loadMetagene)
    updateTextInput(session, "path_load_metagene",  value = pfile$datapath)
    
    ## SOME CHECK : 
    if(nrow(pfile) > 0) {
      f <- load(as.character(pfile$datapath))
      metagene_object <<- eval(as.symbol(f))
      if(ismetagene(metagene_object)) {
        
        createAlert(session, "load_alert", "load_ok", title = "Loading status",
                    content = "Metagene object successfully loaded", append = FALSE, dismiss = TRUE,
                    style = "info")
        updateButton(session, inputId = "go2design", disabled = FALSE)
        
      } else {
        
        createAlert(session, "load_alert", "load_fail", title = "Loading status",
                    content = "The loaded R object seems not be a metagene object", append = FALSE, dismiss = TRUE,
                    style="warning")
      }
      
    }
  })
  
  
  ### LOAD FILE TO CREATE METAGENE OBJECT
  udpateBamList <- reactive({
    shinyFileChoose(input, 'bams', session=session, roots=roots, filetypes=c('bam'))
    pfile = parseFilePaths(roots, input$bams)
    bams <<- c(bams, as.character(pfile$datapath))
    bams <<- unique(bams)
    paste(bams, collapse = '\n')
    
  })
  
  output$bam_list <- renderText({
    udpateBamList()
  })
  
  udpateBedList <- reactive({
    shinyFileChoose(input, 'beds', session=session, roots=roots, filetypes=c('bed'))
    pfile = parseFilePaths(roots, input$beds)
    beds <<- c(beds, as.character(pfile$datapath))
    beds <<- unique(beds)
    paste(beds, collapse = '\n')
    
  })
  
  output$bed_list <- renderText({
    udpateBedList()
  })
  
  #   clear_bam <- eventReactive(input$clear_bams, {
  #     print(bams)
  #     bams <<- NULL
  #     print(bams)
  #   })
  # 
  #   observe({
  #     clear_bam()
  #   })
  
  
  
  ### RUN METAGENE
  runMetagene <- eventReactive(input$runMetagene,{
    ## make some check
    can_run <- TRUE
    print(bams)
    print(beds)
    msg <- "Please wait. This process may take a few minutes..."
    
    ## CHECK 1 : BAM & BED provided ?
    if(length(bams) < 1 || length(beds) < 1) {
      msg <- "You should provide at least one BAM file and one BED file"
      can_run <- FALSE
      createAlert(session, "run_alert", "run_fail", title = "Metagene status",
                  content = msg, append = FALSE, dismiss = TRUE,
                  style="danger")
      return(NULL)
    }
    
    ## CHECK 2 : ALL BAM indexed ?
    unindexed_bams <- c()
    for(bam in bams) {
      if(!file.exists(paste0(bam,".bai"))) {
        unindexed_bams <- c(unindexed_bams, bam)
      }
    }
    
    if(length(unindexed_bams) > 0) {
      msg <- paste(unindexed_bams, collapse = ",")
      createAlert(session, "run_alert", "run_fail2", title = "Metagene status",
                  content = paste0("Can't find indexes for the following BAM files : ", msg), append = FALSE, dismiss = TRUE,
                  style="dander")
      return(NULL)
    }
    
    ## CHECK 3 : ALL REGIONS present in BAM (at least the chrom)
    bed_levels <- getBEDlevels(beds)
    unicomplete_bams <- c()
    
    for(bam in bams) {
      bam_levels <- getBAMlevels(bam)
      if(sum(bed_levels %in% bam_levels) != length(bed_levels)) {
        unicomplete_bams <- c(unicomplete_bams, bam)
      }
    }
    
    if(length(unicomplete_bams) > 0) {
      msg <- paste(unicomplete_bams, collapse = ",")
      createAlert(session, "run_alert", "run_fail3", title = "Metagene status",
                  content = paste0("Some of the wanted regions can't be found in the following BAM files : ", msg), append = FALSE, dismiss = TRUE,
                  style="danger")
      return(NULL)
    }
    
    
    ## CHECK 4 : REGIONS SIZE
    sizes <<- getBEDregionSize(beds)
    
    if(sizes$min != sizes$max){ # WE NEED TO RESIZE 
      can_run = FALSE
      msg <- paste0("All your regions must have the same size. Do you want to shrink all the regions to ", sizes$min, 
                    "bp or extand to ", sizes$max, "bp ?")
      createAlert(session, anchorId = "run_alert", alertId = "run_fail4", title = "Metagene status",
                  content = paste0("Some of the wanted regions can't be found in the following BAM files : ", msg), append = FALSE, dismiss = TRUE,
                  style="warning")
      
      output$resizingUI <- renderUI({
        list(
          fluidRow(
            column(width = 3,
                   bsButton(inputId = "shrink", label = paste0("Shrink to ",sizes$min), icon = icon("cut"))
            ),
            column(width = 8,
                   bsButton(inputId = "extand", label = paste0("Extand to ",sizes$max), icon = icon("arrows-h"))
            )
          )
        )
      })
      
      #updateButton(session,inputId = "runMetagene", disabled = TRUE)
      return(NULL)
    }
    
    if(can_run) {
      closeAlert(session, alertId = "run_alert")
      createAlert(session, "run_alert", "run_succeed", title = "Metagene status",
                  content = msg, append = FALSE, dismiss = TRUE,
                  style="info")
      
      updateButton(session,inputId = "runMetagene", "Metagene running...", disabled = TRUE)
      
      tryCatch ({
        metagene_object <<- metagene$new(beds, bams)
        createAlert(session, "run_alert", "run_succeed1", title = "Metagene status",
                    content = "Done running metagene !", append = FALSE, dismiss = TRUE,
                    style = "success")
        updateButton(session, inputId = "go2design", disabled = FALSE)
        #updateButton(session, inputId = "saveMetagene", disabled = FALSE)
        
      },error = function(e) {
        print(e)
      },warning = function(w) {
        print(w)
      })
      
      
    }
  })
  
  observe({
    runMetagene()
  })
  
  
  extandRegionsAndRunMetagene <- eventReactive(input$extand,{
    updateButton(session,inputId = "extand", "Extanding...", disabled = TRUE)
    updateButton(session,inputId = "shrink", disabled = TRUE)
    updateButton(session,inputId = "runMetagene", "Metagene running...", disabled = TRUE)
    regions <- resizeRegions(beds,sizes$max)
    
    tryCatch ({
      metagene_object <<- metagene$new(regions, bams)
      createAlert(session, "run_alert", "run_succeed2", title = "Metagene status",
                  content = "Done running metagene !", append = FALSE, dismiss = TRUE,
                  style = "success")
      updateButton(session, inputId = "go2design", disabled = FALSE)
      #updateButton(session, inputId = "saveMetagene", disabled = FALSE)
      
    },error = function(e) {
      print(e)
    },warning = function(w) {
      print(w)
    })
    
  })
  
  
  shrinkRegionsAndRunMetagene <- eventReactive(input$shrink,{
    updateButton(session,inputId = "shrink", "Shrinking...", disabled = TRUE)
    updateButton(session,inputId = "extand", disabled = TRUE)
    updateButton(session,inputId = "runMetagene", "Metagene running...", disabled = TRUE)
    
    regions <- resizeRegions(bed_files = beds,new_size = sizes$min)
    
    tryCatch ({
      metagene_object <<- metagene$new(regions, bams)
      createAlert(session, "run_alert", "run_succeed3", title = "Metagene status",
                  content = "Done running metagene !", append = FALSE, dismiss = TRUE,
                  style = "success")
      updateButton(session, inputId = "go2design", disabled = FALSE)
      #updateButton(session, inputId = "saveMetagene", disabled = FALSE)
      
    },error = function(e) {
      print(e)
    },warning = function(w) {
      print(w)
    })
  })
  
  observe({
    extandRegionsAndRunMetagene()
  })
  
  observe({
    shrinkRegionsAndRunMetagene()
  })
  
  
  goToDesign <- eventReactive(input$go2design,{
    print(names(metagene_object$coverages))
  })
  
  observe({
    goToDesign()
  })
  
  output$saveMetagene <- downloadHandler(
    filename = function() { 
      paste("metagene_",format(Sys.time(), "%m_%d_%y_%H_%M_%S"),'.Rda', sep='') 
    },
    content = function(file) {
      save(metagene_object, file = file)
    }
  )
  
  observeEvent(input$go2design, {
    # renvoyer sur le panel principal
    session$sendCustomMessage("myCallbackHandler", "1")
  })
  
  
  
  
})
