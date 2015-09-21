# IMETAGENE SERVER.R
#library(metagene)
#library(d3heatmap)
#library(ggplot2)
#source("helper.R")



shinyServer(function(input, output, session) {
  
  # GLOBALS
  roots <- c(wd='./demo')
  bams <- NULL
  beds <- NULL
  sizes <- 0
  metagene_object <- NULL
  design <- NULL
  
  
  # create empty plot
  waiting_plot <- function(msg) {
    df = data.frame(x=c(1), 
                    y=c(1), 
                    name = c(msg))
    g = ggplot(data=df, mapping=aes(x=x, y=y)) +
      geom_blank() + ylab("") + xlab("") + 
      geom_text(aes(x = x, y = y, label=name), size = 7, color = "midnightblue") +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
    g
  }
  
  ###############################################################
  # Copie de MHmakeRandomString 
  # https://ryouready.wordpress.com/2008/12/18/generate-random-string-name/
  # makeRandomID(n, length)
  # function generates a random string random string of the
  # length (length), made up of numbers, small and capital letters
  
  makeRandomID <- function(n=1, lenght=12)
  {
    randomString <- c(1:n)                  # initialize vector
    for (i in 1:n)
    {
      randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                      lenght, replace=TRUE),
                               collapse="")
    }
    return(randomString)
  }
  
  #  > MHmakeRandomString()
  #  [1] "XM2xjggXX19r"
  ###############################################################
  randomID <- makeRandomID()
  
  #### INPUTS ####
  ### LOAD EXISTING METAGENE
  shiny::observe({
     shinyFiles::shinyFileChoose(input, 'loadMetagene', session=session, roots=roots,
                    filetypes=c('Rda'))
    pfile = shinyFiles::parseFilePaths(roots, input$loadMetagene)
    shiny::updateTextInput(session, "path_load_metagene",  value = pfile$datapath)
    
    ## SOME CHECK : 
    if(nrow(pfile) > 0) {
      f <- load(as.character(pfile$datapath))
      metagene_object <<- eval(as.symbol(f))
      ret <- ismetagene(metagene_object)
      if(ret == 0) {
        
        shinyBS::createAlert(session, "load_alert", "load_ok", title = "Loading status",
                    content = "Metagene object successfully loaded", append = FALSE, dismiss = TRUE,
                    style = "info")
        shinyBS::updateButton(session, inputId = "go2design", disabled = FALSE)
        
        bin_def_size = 100
        regions_size = as.numeric(unlist(width(metagene_object$regions))[1])
        bin_num =  regions_size %/% bin_def_size
        shiny::updateNumericInput(session, inputId = "bin_count", value = bin_num, min = 1, max = regions_size)
        shiny::updateSelectInput(session, inputId = "plot_regions", choices = names(metagene_object$regions), 
                          selected = names(metagene_object$regions))
        
      } else {
        if(ret == 1) {
	  shinyBS::createAlert(session, "load_alert", "load_fail", title = "Loading status",
                      content = "The loaded R object seems not be a metagene object", append = FALSE, dismiss = TRUE,
                      style="warning")
        } else {
	  shinyBS::createAlert(session, "load_alert", "load_fail_old", title = "Loading status",
                      content = "The loaded metagene object has been generated with an old version of metagene (needed version of metagene 2.2.0).", append = FALSE, dismiss = TRUE,
                      style="warning")
        }
      }
      
    }
  })
  
  ### LOAD FILE TO CREATE METAGENE OBJECT
  udpateBamList <- shiny::reactive({
      shinyFiles::shinyFileChoose(input, 'bams', session=session, roots=roots, filetypes=c('bam'))
    pfile = shinyFiles::parseFilePaths(roots, input$bams)
    bams <<- c(bams, as.character(pfile$datapath))
    bams <<- unique(bams)
    
    list(
      shiny::selectInput(inputId = "edit_bams_list", 
                  label = "",
                  selected = bams, 
                  choices = bams,
                  multiple = TRUE,
                  width = '100%',
                  selectize = TRUE)
    )  
  })
  
  output$bam_list <- shiny::renderUI({
    udpateBamList()
  })
  
  udpateBedList <- shiny::reactive({
      shinyFiles::shinyFileChoose(input, 'beds', session=session, roots=roots, filetypes=c('bed'))
    pfile = shinyFiles::parseFilePaths(roots, input$beds)
    beds <<- c(beds, as.character(pfile$datapath))
    beds <<- unique(beds)
    
    list(
      shiny::selectInput(inputId = "edit_beds_list", 
                  label = "",
                  selected = beds, 
                  choices = beds,
                  multiple = TRUE,
                  width = '100%',
                  selectize = TRUE)
    )
    
  })
  
  output$bed_list <- shiny::renderUI({
    udpateBedList()
  })
  
  shiny::observe({
    if(is.null(input$edit_bams_list))
      return(NULL)
    bams <<- input$edit_bams_list
    
  })
  
  shiny::observe({
    if(is.null(input$edit_beds_list))
      return(NULL)
    beds <<- input$edit_beds_list
    
  })
  
  ### RUN METAGENE
  runMetagene <- shiny::eventReactive(input$runMetagene,{
    ## make some check
    cat("mes bams : ", bams, "\n")
    cat("mes beds : ", beds, "\n")
    
    can_run <- TRUE
    msg <- "Please wait. This process may take a few minutes..."
    
    ## CHECK 1 : BAM & BED provided ?
    if(length(bams) < 1 || length(beds) < 1) {
      msg <- "You should provide at least one BAM file and one BED file"
      can_run <- FALSE
      shinyBS::createAlert(session, "run_alert", "run_fail", title = "Metagene status",
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
      shinyBS::createAlert(session, "run_alert", "run_fail2", title = "Metagene status",
                  content = paste0("Can't find indexes for the following BAM files : ", msg), append = FALSE, dismiss = TRUE,
                  style="danger")
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
      shinyBS::createAlert(session, "run_alert", "run_fail3", title = "Metagene status",
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
      shinyBS::createAlert(session, anchorId = "run_alert", alertId = "run_fail4", title = "Metagene status",
                  content = paste0("Some of the wanted regions can't be found in the following BAM files : ", msg), append = FALSE, dismiss = TRUE,
                  style="warning")
      
      output$resizingUI <- shiny::renderUI({
        list(
          shiny::fluidRow(
            column(width = 3,
                   shinyBS::bsButton(inputId = "shrink", label = paste0("Shrink to ",sizes$min), icon = icon("cut"))
            ),
            column(width = 8,
                   shinyBS::bsButton(inputId = "extand", label = paste0("Extand to ",sizes$max), icon = icon("arrows-h"))
            )
          )
        )
      })
      
      return(NULL)
    }
    
    if(can_run) {
      shinyBS::closeAlert(session, alertId = "run_alert")
      shinyBS::createAlert(session, "run_alert", "run_succeed", title = "Metagene status",
                  content = msg, append = FALSE, dismiss = TRUE,
                  style="info")
      
      shinyBS::updateButton(session,inputId = "runMetagene", "Metagene running...", disabled = TRUE)
      
      tryCatch ({
        metagene_object <<- metagene$new(beds, bams)
        shinyBS::createAlert(session, "run_alert", "run_succeed1", title = "Metagene status",
                    content = "Done running metagene !", append = FALSE, dismiss = TRUE,
                    style = "success")
        shinyBS::updateButton(session, inputId = "go2design", disabled = FALSE)
        
        bin_def_size = 100
        regions_size = as.numeric(unlist(width(metagene_object$regions))[1])
        bin_num =  regions_size %/% bin_def_size
        shiny::updateNumericInput(session, inputId = "bin_count", value = bin_num, min = 1, max = regions_size)
        shiny::updateSelectInput(session, inputId = "plot_regions", choices = names(metagene_object$regions), 
                          selected = names(metagene_object$regions))
        
      },error = function(e) {
        print(e)
      },warning = function(w) {
        print(w)
      })
      
      
    }
  })
  
  shiny::observe({
    runMetagene()
  })
  
  
  extandRegionsAndRunMetagene <- shiny::eventReactive(input$extand,{
    shinyBS::updateButton(session,inputId = "extand", "Extanding...", disabled = TRUE)
    shinyBS::updateButton(session,inputId = "shrink", disabled = TRUE)
    shinyBS::updateButton(session,inputId = "runMetagene", "Metagene running...", disabled = TRUE)
    regions <- resizeRegions(beds,sizes$max)
    
    tryCatch ({
      metagene_object <<- metagene$new(regions, bams)
      shinyBS::createAlert(session, "run_alert", "run_succeed2", title = "Metagene status",
                  content = "Done running metagene !", append = FALSE, dismiss = TRUE,
                  style = "success")
      shinyBS::updateButton(session, inputId = "go2design", disabled = FALSE)
      
      bin_def_size = 100
      regions_size = as.numeric(unlist(width(metagene_object$regions))[1])
      bin_num =  regions_size %/% bin_def_size
      shiny::updateNumericInput(session, inputId = "bin_count", value = bin_num, min = 1, max = regions_size)
      shiny::updateSelectInput(session, inputId = "plot_regions", choices = names(metagene_object$regions), 
                        selected = names(metagene_object$regions))
      
    },error = function(e) {
      print(e)
    },warning = function(w) {
      print(w)
    })
    
  })
  
  
  shrinkRegionsAndRunMetagene <- shiny::eventReactive(input$shrink,{
    shinyBS::updateButton(session,inputId = "shrink", "Shrinking...", disabled = TRUE)
    shinyBS::updateButton(session,inputId = "extand", disabled = TRUE)
    shinyBS::updateButton(session,inputId = "runMetagene", "Metagene running...", disabled = TRUE)
    
    regions <- resizeRegions(bed_files = beds,new_size = sizes$min)
    
    tryCatch ({
      metagene_object <<- metagene$new(regions, bams)
      shinyBS::createAlert(session, "run_alert", "run_succeed3", title = "Metagene status",
                  content = "Done running metagene !", append = FALSE, dismiss = TRUE,
                  style = "success")
      shinyBS::updateButton(session, inputId = "go2design", disabled = FALSE)
      
      bin_def_size = 100
      regions_size = as.numeric(unlist(width(metagene_object$regions))[1])
      bin_num =  regions_size %/% bin_def_size
      shiny::updateNumericInput(session, inputId = "bin_count", value = bin_num, min = 1, max = regions_size)
      shiny::updateSelectInput(session, inputId = "plot_regions", choices = names(metagene_object$regions),
                        selected = names(metagene_object$regions))
      
    },error = function(e) {
      print(e)
    },warning = function(w) {
      print(w)
    })
  })
  
  shiny::observe({
    extandRegionsAndRunMetagene()
  })
  
  shiny::observe({
    shrinkRegionsAndRunMetagene()
  })
  
  output$saveMetagene <- shiny::downloadHandler(
    filename = function() { 
      paste("metagene_",format(Sys.time(), "%m_%d_%y_%H_%M_%S"),'.Rda', sep='') 
    },
    content = function(file) {
      save(metagene_object, file = file)
    }
  )
  
  shiny::observeEvent(input$go2design, {
    # preparation des donnees pour le design
    
    if(is.null(bams)){
      bams <<- names(metagene_object$coverages)
    }
    
    #bam_choice <- sapply(X=bams, FUN=extract_file_name)
    #names(bam_choice) = NULL
    bam_choice <-  bams
    
    shiny::updateSelectInput(session, inputId = "chips", choices = bam_choice)
    shiny::updateSelectInput(session, inputId = "ctrls", choices = bam_choice)
    
    #update design
    output$current_mg_design = shiny::renderDataTable({
      metagene_object$design
    })
    
    # envoyer sur le panel DESIGN
    session$sendCustomMessage("myCallbackHandler", "1")
  })
  
  
  #### DESIGNS ####
  ### LOAD EXISTING DESIGN
  shiny::observe({
      shinyFiles::shinyFileChoose(input, 'loadDesign', session=session, roots=roots,
                    filetypes=c('','txt','csv','tsv'))
    pfile = shinyFiles::parseFilePaths(roots, input$loadDesign)
    shiny::updateTextInput(session, "path_load_design",  value = pfile$datapath)
    
    ## SOME CHECK : 
    if(nrow(pfile) > 0) {
      
      tryCatch ({
        design <<- read.table(file = as.character(pfile$datapath), header = input$header)
        metagene_object$add_design(design)
        
        msg <- "Design successfully loaded"
        shinyBS::createAlert(session, "load_alert_d", "load_d_ok", title = "Loading status",
                    content = msg, append = FALSE, dismiss = TRUE,
                    style = "info")
        
        output$current_mg_design = shiny::renderDataTable({
          design
        })
	shinyBS::updateButton(session, inputId = "go2matrix", disabled = FALSE)
        
      },error = function(e) {
        shinyBS::closeAlert(session,alertId = "load_d_fail_error")
        print(e)
        shinyBS::createAlert(session, "load_alert_d", "load_d_fail_error", title = "Loading status",
                    content = as.character(e), append = FALSE, dismiss = TRUE,
                    style="danger")
      },warning = function(w) {
        print(w)
        shinyBS::createAlert(session, "load_alert_d", "load_d_fail_warn", title = "Loading status",
                    content = as.character(w), append = FALSE, dismiss = TRUE,
                    style="warning")
      })
    }
  })
  
  ### DISPLAY OF DESIGN IN LOADED METAGENE OBJECT
  shiny::observeEvent(input$path_load_metagene, {
    output$current_mg_design = shiny::renderDataTable({
      if(!is.null(metagene_object)) {
        if(nrow(metagene_object$design) > 0) {
	  shinyBS::updateButton(session, inputId = "go2matrix", disabled = FALSE)
          metagene_object$design
        }
      }
    })
  })
  
  ### CONSTRUCT NEW DESIGN
  
  shiny::observeEvent(input$chips, {
    if(is.null(bams)){
      bams <<- names(metagene_object$coverages)
    }
    
    #bam_choice <- sapply(X=bams, FUN=extract_file_name)
    #names(bam_choice) = NULL
    bam_choice <-  bams
    
    current_chips_select <- input$chips
    current_ctrls_select <- input$ctrls
    
    control_choice <- bam_choice[! (bam_choice %in% current_chips_select)]
    
    shiny::updateSelectInput(session, inputId = "ctrls", choices = control_choice, selected = current_ctrls_select)
  })
  
  shiny::observeEvent(input$ctrls, {
    if(is.null(bams)){
      bams <<- names(metagene_object$coverages)
    }
    
    #bam_choice <- sapply(X=bams, FUN=extract_file_name)
    #names(bam_choice) = NULL
    bam_choice <-  bams
    
    current_chips_select <- input$chips
    current_ctrls_select <- input$ctrls
    
    chips_choice <- bam_choice[! (bam_choice %in% current_ctrls_select)]
    
    shiny::updateSelectInput(session, inputId = "chips", choices = chips_choice, selected = current_chips_select)
  })
  
  addExperiments <- shiny::eventReactive(input$saveExp, {
    
    ### TO DO : TEST AVANT AJOUT EXPERIENCE !!!
    
    can_add = TRUE
    
    if(can_add) {
      if(is.null(design)) {
        shinyBS::updateButton(session, inputId = "go2matrix", disabled = FALSE)
        design <<- data.frame(Samples = bams)
      }
      
      n <- nrow(design)
      v = rep(0,n)
      names(v) = bams
      
      chips <- input$chips
      ctrls <- input$ctrls
      
      for(chip in chips) {
        v[chip] <- 1         
      }
      
      for(ctrl in ctrls) {
        v[ctrl] <- 2         
      }
      
      design[input$exp_name] <<- v
      
      
      output$current_mg_design = shiny::renderDataTable({
        metagene_object$add_design(design)
        design
      })
      
      shinyBS::updateCollapse(session, "design", open = "Current design")
    }
    
  })
  
  shiny::observe({
    addExperiments()
  })
  
  ### SAVE DESIGN
  output$saveDesign <- shiny::downloadHandler(
    filename = function() {
      paste("metagene_design_",format(Sys.time(), "%m_%d_%y_%H_%M_%S"),'.csv', sep='')
    },
    content = function(file) {
      write.csv(metagene_object$design, file, row.names=FALSE, quote=FALSE )
    }
  )
  
  
  shiny::observeEvent(input$go2matrix, {
    # preparation des donnees pour la matrix
    
    if(is.null(bams)){
      bams <<- names(metagene_object$coverages)
    }
    
    # afficher les données matrice si presente
    if(!is.null(metagene_object) && length(metagene_object$matrices) > 0) {
      ## nom des regions
      m <- metagene_object$matrices
      regions_names <- paste0("Regions : ",paste(names(m), collapse = ", "))
      
      experiences_names_for_all_regions <- c()
      experiences_list <- c()
      for(region_name in names(m)) {
        experiences_list <- c(experiences_list, paste(region_name, names(m[[region_name]])))
        experiences_names <- names(m[[region_name]])
        experiences_names_for_all_regions = c(experiences_names_for_all_regions, 
                                              paste("Experiments for",region_name,":", paste(experiences_names, collapse = ", ")))
      }
      
      output$current_mg_matrix_content <- renderText({
        paste0(regions_names,"\n",paste(experiences_names_for_all_regions, collapse = "\n"))
      })
      
      output$current_mg_matrix_heatmap <- shiny::renderUI({
        list(
          shiny::selectInput(inputId = "select_exp_4_heatmap", 
                      label = "Select an experiment to see its matrix",
                      selected = NULL, 
                      choices = experiences_list,
                      multiple = FALSE,
                      width = '100%',
                      selectize = TRUE),
          d3heatmap::d3heatmapOutput("heatmap")
        )  
      })
    }
    
    
    # envoyer sur le panel MATRIX
    session$sendCustomMessage("myCallbackHandler", "2")
  })
  
  shiny::observeEvent(input$select_exp_4_heatmap, {
    
    selected_exp <- strsplit(x = input$select_exp_4_heatmap, split = " ")[[1]]
    m <- metagene_object$matrices[[selected_exp[1]]][[selected_exp[2]]]$input
    
    color = rev(heat.colors(256))
    
    output$heatmap <- d3heatmap::renderD3heatmap({
      d3heatmap::d3heatmap(m, scale = "column", dendrogram = "none", color = "Blues")
    })
  })
  
  #### MATRIX PARAMETERS
  
  shiny::observeEvent(input$runMatrix, {
    if(!is.null(metagene_object)) {
      
      if(! input$use_design){
        used_design <- NA
      } else {
        used_design <- metagene_object$adesign
      }
      
      if(input$noise == "NONE"){
        noise <- NA
      } else {
        noise <- input$noise
      }
      
      if(input$norm == "NONE"){
        norm <- NA
      } else {
        norm <- input$norm
      }
      
      tryCatch ({
        metagene_object$produce_matrices(select_regions = names(metagene_object$regions), 
                                         design = used_design,
                                         bin_size = input$bin_size, 
                                         noise_removal = noise,
                                         normalization = norm, 
                                         flip_regions = input$flip)
        
        shinyBS::updateButton(session,inputId = "runMatrix", "Producing matrix...", disabled = TRUE)
        shinyBS::updateButton(session,inputId = "go2plot", disabled = FALSE)
        shinyBS::updateButton(session,inputId = "updateMetagene", disabled = FALSE)
        
        shinyBS::createAlert(session, "run_matrix_alert", "run_matrix_success", title = "Produce matrix status",
                    content = "Matrix successfully produced", append = FALSE, dismiss = TRUE,
                    style="success")
        
        ### Construction visualisation de la matrice
        if(!is.null(metagene_object) && length(metagene_object$matrices) > 0) {
          ## nom des regions
          m <- metagene_object$matrices
          regions_names <- paste0("Regions : ",paste(names(m), collapse = ", "))
          
          experiences_names_for_all_regions <- c()
          experiences_list <- c()
          for(region_name in names(m)) {
            experiences_list <- c(experiences_list, paste(region_name, names(m[[region_name]])))
            experiences_names <- names(m[[region_name]])
            experiences_names_for_all_regions = c(experiences_names_for_all_regions, 
                                                  paste("Experiments for",region_name,":", paste(experiences_names, collapse = ", ")))
          }
          
          output$current_mg_matrix_content <- renderText({
            paste0(regions_names,"\n",paste(experiences_names_for_all_regions, collapse = "\n"))
          })
          
          output$current_mg_matrix_heatmap <- shiny::renderUI({
            list(
              shiny::selectInput(inputId = "select_exp_4_heatmap", 
                          label = "Select an experiment to see its matrix",
                          selected = NULL, 
                          choices = experiences_list,
                          multiple = FALSE,
                          width = '100%',
                          selectize = TRUE),
              d3heatmap::d3heatmapOutput("heatmap")
            )  
          })
        }
        
	shinyBS::updateCollapse(session, "matrix", open = "Current matrix")
        
      },error = function(e) {
        shinyBS::closeAlert(session,alertId = "load_d_fail_error")
        print(e)
        shinyBS::createAlert(session, "run_matrix_alert", "run_matrix_error", title = "Produce matrix status",
                    content = as.character(e), append = FALSE, dismiss = TRUE,
                    style="danger")
      },warning = function(w) {
        print(w)
        shinyBS::createAlert(session, "run_matrix_alert", "run_matrix_warn", title = "Produce matrix status",
                    content = as.character(w), append = FALSE, dismiss = TRUE,
                    style="warning")
      })
    }
  })
  
  
  output$updateMetagene <- shiny::downloadHandler(
    filename = function() { 
      paste("metagene_",format(Sys.time(), "%m_%d_%y_%H_%M_%S"),'.Rda', sep='') 
    },
    content = function(file) {
      save(metagene_object, file = file)
    }
  )
  
  
  shiny::observeEvent(input$go2plot, {
    session$sendCustomMessage("myCallbackHandler", "3")
  })
  
  output$mg_plot <- renderPlot({
    return(waiting_plot("Waiting for your request..."))
  })
  
  shiny::observeEvent(input$runPlot, {
    if(!is.null(metagene_object) && length(metagene_object$matrices) > 0) {
      
      
      if(! input$use_design){
        used_design <- NA
      } else {
        used_design <- metagene_object$design
      }
      
      shinyBS::updateButton(session,inputId = "runPlot", "Plotting in progress...", disabled = TRUE)
      
      metagene_plot <- metagene_object$plot(design = used_design, regions_group = input$plot_regions,
                                     bin_size = input$bin_size, alpha = input$alpha, sample_count = input$sample_count,
                                     range = c(-1,1), title = input$plot_title, flip_regions = input$flip)
      
      pdf(paste0("plots/metagene_plot_",randomID,".pdf"), onefile=T, paper="USr")
      print(metagene_plot)
      dev.off()
      
      png(paste0("plots/metagene_plot_",randomID,".png"))
      print(metagene_plot)
      dev.off()
      
      shinyBS::updateButton(session,inputId = "runPlot", "Plotting", disabled = FALSE)
      
      output$mg_plot <- renderPlot({
        metagene_plot
      })
    }
  })
  
  
  output$savePlotPDF <- shiny::downloadHandler(
    filename = function() {
      "metagene_plot.pdf"
    },
    content = function(file) {
      file.copy(paste0("plots/metagene_plot_",randomID,".pdf"), file, overwrite = TRUE)
    }
  )
  
  output$savePlotPNG <- shiny::downloadHandler(
    filename = function() {
      "metagene_plot.png"
    },
    content = function(file) {
      file.copy(paste0("plots/metagene_plot_",randomID,".png"), file, overwrite = TRUE)
    }
  )
  
})
