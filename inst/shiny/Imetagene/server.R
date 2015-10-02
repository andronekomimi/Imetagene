# IMETAGENE SERVER.R
library(metagene)
library(d3heatmap)
library(ggplot2)
source("helper.R")



shinyServer(function(input, output, session) {
  
  # GLOBALS
  roots <- c(wd= path.expand("~"))
  bams <- NULL
  beds <- NULL
  sizes <- 0
  metagene_object <- NULL
  design <- NULL
  tempfiles <- tempfile(pattern = "imetagene_", fileext = c(".pdf",".png"))
  
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
  
  #### INPUTS ####
  ### LOAD EXISTING METAGENE
  shiny::observe({
    shinyFiles::shinyFileChoose(input, 'loadMetagene', session=session, roots=roots,
                                filetypes=c('Rda', 'RData', 'RDA', 'RDATA'))
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
        
        regions_size = as.numeric(unlist(width(metagene_object$get_regions()))[1])
        bin_num =  regions_size %/% bin_def_size
        
        shiny::updateNumericInput(session, inputId = "bin_count", value = bin_num, min = 1, max = regions_size)
        shiny::updateSelectInput(session, inputId = "plot_regions", choices = names(metagene_object$get_regions()), 
                                 selected = names(metagene_object$get_regions()))
        #ICI
        
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
    shinyFiles::shinyFileChoose(
      input, 'beds', session=session, roots=roots, 
      filetypes=c('bed', 'BED', 'narrowPeak', 'broadPeak'))
    
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
    
    uncomplete_bams <- c()
    
    for(bam in bams) {
      bam_levels <- getBAMlevels(bam)
      
      if(sum(bed_levels %in% bam_levels) != length(bed_levels)) {
        uncomplete_bams <- c(uncomplete_bams, bam)
      }
    }
    
    
    if(length(uncomplete_bams) > 0) {
      msg <- paste(uncomplete_bams, collapse = ",")
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
        metagene_object <<- metagene$new(regions = beds, bam_files = bams)
        shinyBS::createAlert(session, "run_alert", "run_succeed1", title = "Metagene status",
                             content = "Done running metagene !", append = FALSE, dismiss = TRUE,
                             style = "success")
        shinyBS::updateButton(session, inputId = "go2design", disabled = FALSE)
        
        bin_def_size = 100
        regions_size = as.numeric(unlist(width(metagene_object$get_regions()))[1])
        bin_num =  regions_size %/% bin_def_size
        
        shiny::updateNumericInput(session, inputId = "bin_count", value = bin_num, min = 1, max = regions_size)
        shiny::updateSelectInput(session, inputId = "plot_regions", choices = names(metagene_object$get_regions()),
                                 selected = names(metagene_object$get_regions()))
        
        #ICI
        
        
      },error = function(e) {
        shinyBS::createAlert(session, anchorId = "run_alert", alertId = "run_fail5", title = "Metagene status",
                             content = as.character(e), append = FALSE, dismiss = TRUE,
                             style="danger")
        print(e)
      },warning = function(w) {
        shinyBS::createAlert(session, anchorId = "run_alert", alertId = "run_fail6", title = "Metagene status",
                             content = as.character(w), append = FALSE, dismiss = TRUE,
                             style="warning")
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
      regions_size = as.numeric(unlist(width(metagene_object$get_regions()))[1])
      bin_num =  regions_size %/% bin_def_size
      
      shiny::updateNumericInput(session, inputId = "bin_count", value = bin_num, min = 1, max = regions_size)
      shiny::updateSelectInput(session, inputId = "plot_regions", choices = names(metagene_object$get_regions()), 
                               selected = names(metagene_object$get_regions()))
      
      #ICI
      
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
      regions_size = as.numeric(unlist(width(metagene_object$get_regions()))[1])
      bin_num =  regions_size %/% bin_def_size
      
      shiny::updateNumericInput(session, inputId = "bin_count", value = bin_num, min = 1, max = regions_size)
      shiny::updateSelectInput(session, inputId = "plot_regions", choices = names(metagene_object$get_regions()),
                               selected = names(metagene_object$get_regions()))
      
      #ICI 
      
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
      bams <<- names(metagene_object$get_raw_coverages())
    }
    
    #bam_choice <- sapply(X=bams, FUN=extract_file_name)
    #names(bam_choice) = NULL
    
    bam_choice <- names(metagene_object$get_normalized_coverages())
    
    shiny::updateSelectInput(session, inputId = "chips", choices = bam_choice)
    shiny::updateSelectInput(session, inputId = "ctrls", choices = bam_choice)
    
    #update design
    output$current_mg_design = shiny::renderDataTable({
      metagene_object$get_design()
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
        design <<- read.table(file = as.character(pfile$datapath), header = input$header, stringsAsFactors = FALSE)
        #         print(names(metagene_object$get_raw_coverages()))
        metagene_object$add_design(design, check_bam_file = TRUE)
        
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
        if(nrow(metagene_object$get_design()) > 0) {
          shinyBS::updateButton(session, inputId = "go2matrix", disabled = FALSE)
          metagene_object$get_design()
        }
      }
    })
  })
  
  ### CONSTRUCT NEW DESIGN
  
  shiny::observeEvent(input$chips, {
    if(is.null(bams)){
      bams <<- names(metagene_object$get_raw_coverages())
    }
    
    #bam_choice <- sapply(X=bams, FUN=extract_file_name)
    #names(bam_choice) = NULL
    bam_choice <-  names(metagene_object$get_raw_coverages())
    
    current_chips_select <- input$chips
    current_ctrls_select <- input$ctrls
    
    control_choice <- bam_choice[! (bam_choice %in% current_chips_select)]
    
    shiny::updateSelectInput(session, inputId = "ctrls", choices = control_choice, selected = current_ctrls_select)
  })
  
  shiny::observeEvent(input$ctrls, {
    if(is.null(bams)){
      bams <<- names(metagene_object$get_raw_coverages())
    }
    
    #bam_choice <- sapply(X=bams, FUN=extract_file_name)
    #names(bam_choice) = NULL
    bam_choice <-  names(metagene_object$get_raw_coverages())
    
    current_chips_select <- input$chips
    current_ctrls_select <- input$ctrls
    
    chips_choice <- bam_choice[! (bam_choice %in% current_ctrls_select)]
    
    shiny::updateSelectInput(session, inputId = "chips", choices = chips_choice, selected = current_chips_select)
  })
  
  addExperiments <- shiny::eventReactive(input$saveExp, {
    
    ### TO DO : TEST AVANT AJOUT EXPERIENCE !!!
    can_add = TRUE
    
    if(is.null(input$chips)) {
      can_add = FALSE
    }
    
    
    if(can_add) {
      if(nrow(metagene_object$get_design()) == 0) {
        
        shinyBS::updateButton(session, inputId = "go2matrix", disabled = FALSE)
        design <<- data.frame(Samples = names(metagene_object$get_raw_coverages()), stringsAsFactors = FALSE)
      }
      
      n <- nrow(design)
      v = rep(0,n)
      names(v) = names(metagene_object$get_raw_coverages())
      
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
        metagene_object$add_design(design, check_bam_file = TRUE)
        design
      })
      
      shinyBS::createAlert(session, "save_exp_alert", "run_matrix_success", title = "Add experiment status",
                           content = "Experiment successfully added to design", append = FALSE, dismiss = TRUE,
                           style="success")
      
      
      shinyBS::updateCollapse(session, "design", open = "Current design")
    } else {
      
      shinyBS::createAlert(session, "save_exp_alert", "run_matrix_fail", title = "Add experiment status",
                           content = "Fail while adding the experiment. Please choose at least one chip", append = FALSE, dismiss = TRUE,
                           style="danger")
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
      write.csv(metagene_object$get_design(), file, row.names=FALSE, quote=FALSE )
    }
  )
  
  
  shiny::observeEvent(input$go2matrix, {
    # preparation des donnees pour la matrix
    
    if(is.null(bams)){
      bams <<- names(metagene_object$get_raw_coverages())
    }
    
    # afficher les données matrice si presente
    if(!is.null(metagene_object) && !is.null(metagene_object$get_matrices())) {
      ## nom des regions
      m <- metagene_object$get_matrices()
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
    m <- metagene_object$get_matrices()[[selected_exp[1]]][[selected_exp[2]]]$input
    
    m_dim = dim(m)
    rownames(m) <-paste0("region_",seq(1,m_dim[1]))
    colnames(m) <- paste0("position_",seq(1,m_dim[2]))
    
    max_dim_1 <- 50
    max_dim_2 <- 100
    max_dim_prod <- max_dim_1 * max_dim_2
    
    if( (m_dim[1] * m_dim[2]) > max_dim_prod) { # matrice trop volumineuse
      idx_row <- sort(sample(1:m_dim[1], max_dim_1, replace = FALSE))
      idx_col <- sort(sample(1:m_dim[2], max_dim_2, replace = FALSE))
      m <- m[idx_row,idx_col]
    }
    
    color = rev(heat.colors(256))
    
    output$heatmap <- d3heatmap::renderD3heatmap({
      d3heatmap::d3heatmap(m, scale = "column", dendrogram = "none", color = "Blues", 
                           xaxis_font_size = "0px", yaxis_font_size = "0px")
      
    })
  })
  
  #### MATRIX PARAMETERS
  
  shiny::observe({
    if(!input$use_design) {
      shiny::updateSelectInput(session, inputId = "noise", choices = "NONE", selected = "NONE")
    } else {
      shiny::updateSelectInput(session, inputId = "noise", choices = c("NONE", "NCIS"), selected = "NONE")
    }
  })
  
  
  shiny::observe({
    
    if(!is.null(metagene_object)) {
      
      m <- metagene_object$get_matrices()
      
      if(!is.null(m)){
        experiences_list <- c()
        for(region_name in input$plot_regions) {
          experiences_names <- names(m[[region_name]])
          experiences_list <- c(experiences_list, experiences_names)
        }
        
        experiences_list = unique(experiences_list)
        
        shiny::updateSelectInput(session, inputId = "plot_exps", choices = experiences_list, selected = experiences_list)
      }
    }
  })
  
  shiny::observeEvent(input$plot_regions, {
    
    m <- metagene_object$get_matrices()
    
    experiences_list <- c()
    for(region_name in input$plot_regions) {
      experiences_names <- names(m[[region_name]])
      experiences_list <- c(experiences_list, experiences_names)
    }
    
    experiences_list = unique(experiences_list)
    
    shiny::updateSelectInput(session, inputId = "plot_exps", choices = experiences_list, selected = experiences_list)
  })
  
  shiny::observeEvent(input$runMatrix, {
    if(!is.null(metagene_object)) {
      
      if(! input$use_design || nrow(metagene_object$get_design()) == 0) {
        used_design <- NA
      } else {
        used_design <- metagene_object$get_design()
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
        shinyBS::updateButton(session,inputId = "runMatrix", "Producing matrix...", disabled = TRUE)
        
        metagene_object$produce_matrices(design = used_design,
                                         bin_size = input$bin_size, 
                                         noise_removal = noise,
                                         normalization = norm, 
                                         flip_regions = input$flip)
        
        
        shinyBS::updateButton(session,inputId = "runMatrix", "Produce matrix", disabled = FALSE)
        shinyBS::updateButton(session,inputId = "go2plot", disabled = FALSE)
        shinyBS::updateButton(session,inputId = "updateMetagene", disabled = FALSE)
        
        shinyBS::createAlert(session, "run_matrix_alert", "run_matrix_success", title = "Produce matrix status",
                             content = "Matrix successfully produced", append = FALSE, dismiss = TRUE,
                             style="success")
        
        ### Construction visualisation de la matrice
        if(!is.null(metagene_object) && ! is.null(metagene_object$get_matrices())) {
          ## nom des regions
          m <- metagene_object$get_matrices()
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
          
          shinyBS::updateCollapse(session, "matrix", open = "Current matrix subset")
          
        }
        
      },error = function(e) {
        shinyBS::closeAlert(session,alertId = "load_d_fail_error")
        print(e)
        shinyBS::createAlert(session, "run_matrix_alert", "run_matrix_error", title = "Produce matrix status",
                             content = as.character(e), append = FALSE, dismiss = TRUE,
                             style="danger")
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
    if(!is.null(metagene_object) && !is.null(metagene_object$get_matrices()))  {
      
      
      
      if(! input$use_design){
        used_design <- NA
      } else {
        used_design <- metagene_object$get_design()
      }
      
      shinyBS::updateButton(session,inputId = "runPlot", "Plotting in progress...", disabled = TRUE)
      
      metagene_object$produce_data_frame(alpha = input$alpha, sample_count = input$sample_count)
      
      metagene_object$plot(region_names = input$plot_regions, exp_names = input$plot_exps, range = c(-1,1), title = input$plot_title)
      
      #       metagene_plot <- metagene_object$plot(design = used_design, regions_group = input$plot_regions,
      #                                      bin_size = input$bin_size, alpha = input$alpha, sample_count = input$sample_count,
      #                                      range = c(-1,1), title = input$plot_title, flip_regions = input$flip)
      
      pdf(tempfiles[1], onefile=T, paper="USr")
      print(metagene_object$get_plot())
      dev.off()
      
      png(tempfiles[2])
      print(metagene_object$get_plot())
      dev.off()
      
      shinyBS::updateButton(session,inputId = "runPlot", "Plotting", disabled = FALSE)
      
      output$mg_plot <- renderPlot({
        metagene_object$get_plot()
      })
    }
  })
  
  
  output$savePlotPDF <- shiny::downloadHandler(
    filename = function() {
      "metagene_plot.pdf"
    },
    content = function(file) {
      file.copy(tempfiles[1], file, overwrite = TRUE)
    }
  )
  
  output$savePlotPNG <- shiny::downloadHandler(
    filename = function() {
      "metagene_plot.png"
    },
    content = function(file) {
      file.copy(tempfiles[2], file, overwrite = TRUE)
    }
  )
  
  shinyBS::addPopover(session, id = "hloadMetagene", title = "", content =
                        "File must be a RData file (.RData, .Rda, .RDATA, .RDA) produced by Imetagene or metagene..")
  shinyBS::addPopover(session, id = "hbams", title = "", content =
                        "There is no hard limit in the number of BAM files that can be included in an analysis
                                             (but with too many BAM files, memory may become an issue). BAM files must be indexed
                                              For instance, if you use a file names file.bam, a file named file.bam.bai must be present in the same directory.")
  shinyBS::addPopover(session, id = "hbeds", title = "", content =
                        "To compare custom regions of interest, it is possible to use a list of one or more BED files.
                                                                       BED, narrowPeak and broadPeak format are supported.")
  shinyBS::addPopover(session, id = "hloadDesign", title = "", content =
                        "A design file is a tab-delimited file (.tsv) that describes one or more experiments. An experiment can contain one or more replicates and controls.
                                                                       The first column of a design file contains the list of bam names available in the current metagene object.
                                                                       The following columns correspond to the experiments. They must have a unique name and the possible values are 0 (ignore), 1 (chip) or 2 (control).")
  shinyBS::addPopover(session, id = "hbin_size", title = "", content =
                        "The size of each bin, in nucleotides")
  shinyBS::addPopover(session, id = "hnoise", title = "", content =
                        "Algorithm use to remove background. Requires a design file.")
  shinyBS::addPopover(session, id = "hnorm", title = "", content =
                        "Algorithm to use to normalize the samples.")
  shinyBS::addPopover(session, id = "hflip", title = "", content =
                        "Flip the regions that are on the negative strand.")
  shinyBS::addPopover(session, id = "huse_design", title = "", content =
                        "Use the design produced in the previous section.")
  shinyBS::addPopover(session, id = "halpha", title = "", content =
                        "The interval to display with a ribbon on the metagene plot. Corresponds to the confidence interval in the estimation of the mean using bootstrap.")
  shinyBS::addPopover(session, id = "hsample_count", title = "", content =
                        "The number of draw for the bootstrap analysis.")
})
