## app.R ##
library(shiny)
library(shinydashboard)

## working directory  ##
primersDesign_wd <- getwd() 

## Source for Primer Design  ##
source("generalDesign.R")
source("ReferenceGenome.R")

## libraries for primer design  ##
library(devtools)
library(rBLAST)
library(rtracklayer)
library(BSgenome)
library(Biostrings)

## libraries for ePCR  ##
library(tidyr)
library(dplyr)

## tooltips ##
library(shinyBS)
library(tippy)

## for ePCR ##
library(xml2)
library(httr)
library(httr)
library(stringr)
library(ggplot2)

dbHeader <- dashboardHeader(title = "EpiPrimer")

server <- function(input, output, session) {
  ### global variable for Primersettings
  def_settings <<- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, 5, 12, 0, 10, 0, 0, 0, 0, NA, NA,  NA, NA,    NA, 18, 25, 50, 60, 3,  200, 500, 0, 0, NA, NA,  30, "genomic")
  
  ### Data import:
  Dataset <- reactive({
    if (is.null(input$file)) {
      # User has not uploaded a file yet
      return(data.frame())
    }
    
    read_input_file <- read.delim(input$file$datapath)
    return(read_input_file)
  })
  
  ###### displaying input file #############
  
  output$table <- DT::renderDataTable({
    
    return(Dataset())
  },
  extensions = 'FixedHeader',
  options = list(fixedHeader = TRUE,
                 scrollY = TRUE),
  fillContainer = T,
  class = "display"
  )
  
  output$state <- eventReactive(input$action, {
    primerDesign_wd <- setwd(primersDesign_wd)
    
    if (is.null(input$file)){
      sprintf("Please upload an input file!")
    }
    
    else if(file.exists(paste(getwd(),input$name,sep="/"))){ #dir.exists(paste(getwd(),"Data",input$name,sep="/") 
      "Dataset already exists!"
      sprintf("Dataset %s already exists, please choose other ID!", input$name)
    } else {
      
      
      ww <-showModal(modalDialog(
        title = "We are currently computing your primers!",
        sprintf(paste0("You will be notified when the computation of your primers is finished. Please be patient, this might take a few minutes."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      
      # get current settings
      settings_for_pipeline <- def_settings
      
      #get input values for Forward and Reverse adapters
      if (!is.null(input$adapterForward) && !is.null(input$adapterReverse)){
        adaF <- input$adapterForward
        adaR <- input$adapterReverse
        # now set adapters, only if checkboxes were activated
        settings_for_pipeline[20] = adaF
        settings_for_pipeline[21] = adaR
      } else {
        settings_for_pipeline[20] = NA
        settings_for_pipeline[21] = NA
      }
      
      check_snps = FALSE
      if(!is.null(input$i_snps.amplicon) && !is.null(input$i_snps.primer1) && !is.null(input$i_snps.primer2)) {
        if (input$i_snps.amplicon != 0 || input$i_snps.primer1 != 0 || input$i_snps.primer2 != 0){
          check_snps = TRUE
        } else {
          check_snps = FALSE
        }
      }
      
      # now set checks for snps and repeats
      settings_for_pipeline[3] = check_snps
      settings_for_pipeline[4] = FALSE
      
      #call the primer design pipeline
      primer.design.pipeline.refactored(Dataset(),
                             path.out = paste(getwd(),input$name,sep="/"),
                             primer.type = input$i_primer_type,
                             low.complexity.primer.removal = as.logical(settings_for_pipeline[1]),
                             remove.primers.with.n = as.logical(settings_for_pipeline[2]),
                             check4snps = as.logical(settings_for_pipeline[3]),
                             check4repeats = as.logical(settings_for_pipeline[4]),
                             allow.repeats.in.primers = as.logical(settings_for_pipeline[5]),
                             allow.repeats.in.amplicon = as.logical(settings_for_pipeline[6]),
                             annotate.genes = as.logical(settings_for_pipeline[7]),
                             annotate.cpg.islands = as.logical(settings_for_pipeline[8]),
                             create.toplist = as.logical(settings_for_pipeline[9]),
                             max.bins.low.complexity = as.numeric(settings_for_pipeline[10]),
                             primer.align.binsize = as.numeric(settings_for_pipeline[11]),
                             min.snps.amplicon = as.numeric(settings_for_pipeline[12]),
                             max.snps.amplicon = as.numeric(settings_for_pipeline[13]),
                             min.snps.primer1 = as.numeric(settings_for_pipeline[14]),
                             max.snps.primer1 = as.numeric(settings_for_pipeline[15]),
                             min.snps.primer2 = as.numeric(settings_for_pipeline[16]),
                             max.snps.primer2 = as.numeric(settings_for_pipeline[17]),
                             hp.length.min = as.numeric(settings_for_pipeline[18]),
                             hp.length.max = as.numeric(settings_for_pipeline[19]),
                             add.ngs.adaptor.f = settings_for_pipeline[20],
                             add.ngs.adaptor.r = settings_for_pipeline[21],
                             strand = settings_for_pipeline[22],
                             min.length.primer = as.numeric(settings_for_pipeline[23]),
                             max.length.primer = as.numeric(settings_for_pipeline[24]),
                             min.Tm.primer = as.numeric(settings_for_pipeline[25]),
                             max.Tm.primer = as.numeric(settings_for_pipeline[26]),
                             max.Tm.difference.primer = as.numeric(settings_for_pipeline[27]),
                             min.length.amplicon = as.numeric(settings_for_pipeline[28]),
                             max.length.amplicon = as.numeric(settings_for_pipeline[29]),
                             min.number.gc.amplicon = as.numeric(settings_for_pipeline[30]),
                             min.number.cg.amplicon = as.numeric(settings_for_pipeline[31]),
                             min.C2T.primer1 = as.numeric(settings_for_pipeline[32]),
                             min.G2A.primer2 = as.numeric(settings_for_pipeline[33]),
                             chop.size = as.numeric(settings_for_pipeline[34])
      )
      
      sprintf("Finished Computation!!")
      ww <-showModal(modalDialog(
        title = "Primers are READY!",
        sprintf(paste0("You can find details concerning your run in the Menueitem 'Results of Primer Design' and further analyze them using the created files in the folder %s"),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      
      sprintf("Finished Computation!!")
    }})
  
  ################ download example files for sequence and region input ####################
  
  output$downloadSequenceFile <- downloadHandler(
    filename = function() {
      paste("ExampleSequenceFile", "txt", sep=".")
    },
    content = function(file){
      file_path <- file.path(primersDesign_wd, "input_sequences.txt", fsep=.Platform$file.sep)
      file.copy(file_path, file)
    },
    contentType = "txt"
  )
  
  output$downloadRegionsFile <- downloadHandler(
    filename = function() {
      paste("ExampleRegionsFile", "txt", sep=".")
    },
    content = function(file){
      file_path <- file.path(primersDesign_wd, "input.txt", fsep=.Platform$file.sep)
      file.copy(file_path, file)
    },
    contentType = "txt"
  )
  
  ################ download the results of Primers Design ####################
  
  output$downloadPrimer<- downloadHandler(
    filename = function() {
      paste(input$name,"zip",sep=".")
    },
    
    content <- function(con) {
      ss <- paste(getwd(),input$name,sep="/")
      tmpdir <- tempdir()
      setwd(tempdir())
      filesToSave <- c(ss) #List to hold paths to your files in shiny
      #Put all file paths inside filesToSave...
      zz <- zip(zipfile=con, files = filesToSave,flags = "-r9X", extras = "",zip ="zip")
      return(zz)
    },
    contentType = "results/zip"
  )
  
  ################ download the results of ePCR ####################
  
  output$downloadePCRresults<- downloadHandler(
    filename = function() {
      paste(input$blast_id, "zip", sep=".")
    },
    content <- function(con) {
      ss <- file.path(primersDesign_wd, "ePCR", input$blast_id, fsep=.Platform$file.sep)
      tmpdir <- tempdir()
      setwd(tempdir())
      filesToSave <- c(ss) #List to hold paths to your files in shiny
      #Put all file paths inside filesToSave...
      zz <- zip(zipfile=con, files = filesToSave, flags = "-r9X", extras = "", zip ="zip")
      return(zz)
    },
    contentType = "results/zip"
  )
  
  ############# display the top list of primers ########### 
  
  showToplist <- reactive({ if (!input$toplist) {return(NULL)}
    files <- data.frame(results=list.files(paste(primersDesign_wd,input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("toplist",files[["results"]])])
    if(length(file_path) == 0){
      ww <-showModal(modalDialog(
        title = "No Primers in the Top List!",
        sprintf(paste0("Unfortunateley, we could not find any primers in the Top List."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      stop("No primers contained in the Top List!")
    }
    
    print(file_path)
    read_file <- read.delim(file_path)
    
    print(read_file)
    return(read_file)
  })
  
  output$viewtoplist <- DT::renderDataTable({
    if (!input$toplist) {return(data.frame())}
    showToplist()
  },
  extensions = 'FixedHeader',
  options = list(fixedHeader = TRUE,
                 scrollY = TRUE),
  fillContainer = T,
  class = "display"
  )
  
  
  ############# display the whole list of primers ###########
  
  showWholelist <- reactive({ if (!input$wholelist) {return(NULL)}
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(results=list.files(paste(primersDesign_wd,input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("wholelist",files[["results"]])])
    if(length(file_path) == 0){
      ww <-showModal(modalDialog(
        title = "No Primers in the Whole List!",
        sprintf(paste0("Unfortunateley, we could not find any primers in the Whole List."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      stop("No primers contained in the Whole List!")
    }
    
    print(file_path)
    read_file <- read.delim(file_path)
    
    print(read_file)
    return(read_file)
  })
  
  output$viewwholelist <- DT::renderDataTable({
    if (!input$wholelist) {return(data.frame())}
    showWholelist()
  },
  extensions = 'FixedHeader',
  options = list(fixedHeader = TRUE,
                 scrollY = TRUE),
  fillContainer = T,
  class = "display"
  )
  
  ######### displaying selected primers table #######
  #also generate fasta files for forward and reverse primers 
  
  showSelectlist <- reactive({ if (!input$selectlist) {return(NULL)}
    s = input$viewwholelist_rows_selected
    whole_selected_list <- showWholelist()[s,]
    write.table(whole_selected_list,file=paste0(primersDesign_wd,"/",input$name,"/","SelectedPrimers", ".txt"),quote = FALSE,col.names = TRUE,row.names=FALSE, sep = "\t")
    
    short_selected_list <- showWholelist()[s,c("amplicon.id","primer1.sequence","primer2.sequence")]
    
    Fprimers_qc <- write.table(sapply(1:nrow(short_selected_list),function(x){
      paste0(paste0(">",short_selected_list[x,"amplicon.id"]),'\n',paste0(short_selected_list[x,"primer1.sequence"]),'\n')
    }),file=paste0(primersDesign_wd,"/",input$name,"/","Fprimers.fasta"),quote = FALSE,col.names = FALSE,row.names=FALSE, sep = "\t")
    
    Rprimers_qc <- write.table(sapply(1:nrow(short_selected_list),function(x){
      paste0(paste0(">",short_selected_list[x,"amplicon.id"]),'\n',paste0(short_selected_list[x,"primer2.sequence"]),'\n')
    }),file=paste0(primersDesign_wd,"/",input$name,"/","Rprimers.fasta"),quote = FALSE,col.names = FALSE,row.names=FALSE, sep = "\t")
    
    primers_list_ReadsExtraction <- short_selected_list
    colnames(primers_list_ReadsExtraction) <- c("names","Fprimer","Rprimer")
    PrimersTable_ReadsExtraction <- write.table(primers_list_ReadsExtraction,file=paste0(primersDesign_wd,"/",input$name,"/","PrList_ReadsExtraction.txt"),quote = FALSE,col.names = TRUE,row.names=FALSE, sep = "\t")
    
    return(whole_selected_list)
    
  })
  
  output$viewSelectlist <- DT::renderDataTable({
    return(showSelectlist())
  },
  extensions = 'FixedHeader',
  options = list(fixedHeader = TRUE,
                 scrollY = TRUE),
  fillContainer = T,
  class = "display"
  )
  
  ############# display the black list of primers ###########
  
  showBlacklist <- reactive({ if (!input$blacklist) {return(NULL)}
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(results=list.files(paste(primersDesign_wd,input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("blacklist",files[["results"]])])
    if(length(file_path) == 0){
      ww <-showModal(modalDialog(
        title = "No Primers in the Black List!",
        sprintf(paste0("Unfortunateley, we could not find any primers in the Black List."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      stop("No primers contained in the Black List!")
    }
    print(file_path)
    read_file <- read.delim(file_path)
    
    print(read_file)
    return(read_file)
  })
  
  output$viewblacklist <- DT::renderDataTable({
    if (!input$blacklist) {return(data.frame())}
    showBlacklist()
  },
  extensions = 'FixedHeader',
  options = list(fixedHeader = TRUE,
                 scrollY = TRUE),
  fillContainer = T,
  class = "display"
  )
  
  ############# display the white list of primers ###########
  showWhitelist <- reactive({ if (!input$whitelist) {return(NULL)}
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(results=list.files(paste(primersDesign_wd,input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("white",files[["results"]])])
    if(length(file_path) == 0){
      ww <-showModal(modalDialog(
        title = "No Primers Found on the White List!",
        sprintf(paste0("Unfortunateley, we could not find any primers in the White List."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      stop("No primers were found on the White List!")
    }
    print(file_path)
    read_file <- read.delim(file_path)
    
    print(read_file)
    return(read_file)
  })
  
  output$viewwhitelist <- DT::renderDataTable({
    if (!input$whitelist) {return(data.frame())}
    showWhitelist()
  },
  extensions = 'FixedHeader',
  options = list(fixedHeader = TRUE,
                 scrollY = TRUE),
  fillContainer = T,
  class = "display"
  )
  
  ############# display the log file of the primer design ###########
  
  showLogfile <- reactive({ if (!input$logfile) {return(NULL)}
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(results=list.files(paste(primersDesign_wd,input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("logfile",files[["results"]])])
    if(length(file_path) == 0){
      ww <-showModal(modalDialog(
        title = "No Logfile Found!",
        sprintf(paste0("Unfortunateley, we could not find any Logfile."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      stop("No Logfile found!")
    }
    print(file_path)
    read_file <- read.delim(file_path)
    
    print(read_file)
    return(read_file)
  })
  
  output$viewlogfile <- DT::renderDataTable({
    if (!input$logfile) {return(data.frame())}
    showLogfile()
  })
  
  ############# display the primer design by sequence  ###########
  
  showPrimerdesigns.by.sequence <- reactive({ if (!input$primerdesigns.by.sequence) {return(NULL)}
    files <- data.frame(results=list.files(paste(primersDesign_wd,input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("primerdesigns.by.sequence",files[["results"]])])
    if(length(file_path) == 0){
      ww <-showModal(modalDialog(
        title = "No Primers Found!",
        sprintf(paste0("Unfortunateley, we could not find any primers for your current job. Feel free to check your settings and try again."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      stop("No primers were found!")
    }
    print(file_path)
    read_file <- read.delim(file_path)
    
    print(read_file)
    return(read_file)
  })
  
  output$viewprimerdesigns <- DT::renderDataTable({
    if (!input$primerdesigns.by.sequence) {return(data.frame())}
    showPrimerdesigns.by.sequence()
  })
  
  ############# display the summary of the primer design ###########
  
  showSummary <- reactive({ if (!input$Summary) {return(NULL)}
    files <- data.frame(results=list.files(paste(primersDesign_wd,input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("summary",files[["results"]])])
    if(length(file_path) == 0){
      ww <-showModal(modalDialog(
        title = "No Summary Found!",
        sprintf(paste0("Unfortunateley, we could not find any Summary for your Primer Design Job."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      stop("No Summary found!")
    }
    print(file_path)
    read_file <- read.delim(file_path)
    
    print(read_file)
    return(read_file)
  })
  
  output$viewSummary <- DT::renderDataTable({
    if (!input$Summary) {return(data.frame())}
    showSummary()
  })
  
  ############# display the settings of the primer design ###########
  
  showSettings <- reactive({ if (!input$settings) {return(NULL)}
    files <- data.frame(results=list.files(paste(primersDesign_wd,input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("settings",files[["results"]])])
    if(length(file_path) == 0){
      ww <-showModal(modalDialog(
        title = "No Settings Found!",
        sprintf(paste0("Unfortunateley, we could not find any settings for your current job. "),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      stop("No Settings were found!")
    }
    print(file_path)
    read_file <- read.delim(file_path)
    
    print(read_file)
    return(read_file)
  })
  
  output$viewsettings <- DT::renderDataTable({
    if (!input$settings) {return(data.frame())}
    showSettings()
  })
  
  #################### display graphics of primer design#############
  
  showGraphics <-eventReactive(input$graphics, {
    files <- data.frame(graphs=list.files(paste(getwd(),input$name,"PrimerAnalysis","graphics",sep="/"),full.names=TRUE, pattern =".png"))
  })
  
  output$plot3 <-renderUI({
    #check, if there are graphics that can be displayed, if not inform the user 
    if(is.data.frame(showGraphics()) && nrow(showGraphics())== 0){
      ww <-showModal(modalDialog(
        title = "No Primers Found!",
        sprintf(paste0("Unfortunateley, we could not find any primers for your current job. Feel free to check your settings and try again."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
    } else {
      ss <- lapply(1:nrow(showGraphics()),function(x){
        out_ui <- paste0("image", x)
        imageOutput(out_ui, height = "800px")
      })
      do.call(tagList,ss)
    }
  })
  
  observe({
    for(x in 1:nrow(showGraphics()))
    {
      local({
        my_x <- x
        out_ui <- paste0("image",my_x)
        output[[out_ui]] <- renderImage({
          list(src = paste(showGraphics()$graphs[my_x]),
               alt = "Image failed to render", style="width: 800px; height: 700px")
        }, deleteFile = FALSE)
      })
    }
  })
  
  
  
  ################ download the selected primers list ####################
  
  output$downloadSelectedPrimers<- downloadHandler(
    filename = function() {
      paste0("SelectedPrimers", ".zip")
    },
    
    content <- function(con) {
      ss <-  paste0(primersDesign_wd,"/",input$name,"/","SelectedPrimers",".txt")
      tmpdir <- tempdir()
      setwd(tempdir())
      filesToSave <- c(ss) #List to hold paths to your files in shiny
      #Put all file paths inside filesToSave...
      zz <- zip(zipfile=con, files = filesToSave,flags = "-r9X", extras = "",zip ="zip")
      return(zz)
      print(getwd())
    }
  )
  
  ################ displaying the results of Primer QC #############
  
  FP <- reactive({
    if (is.null(input$Fprimers)){
      # User has not uploaded a file yet
      return(data.frame())
    }
    # use BioStrings package to display .fasta file properly
    read_Fprimers_file <- readDNAStringSet(as.character(input$Fprimers$datapath))
    Primer.Name = names(read_Fprimers_file)
    Sequence = paste(read_Fprimers_file)
    return(data.frame(Primer.Name, Sequence))
  })
  
  RP <- reactive({
    if (is.null(input$Rprimers)){
      # User has not uploaded a file yet
      return(data.frame())
    }
    #read_Rprimers_file <- read.table(input$Rprimers$datapath)
    # use BioStrings package to display .fasta file properly
    read_Rprimers_file <- readDNAStringSet(as.character(input$Rprimers$datapath))
    Primer.Name = names(read_Rprimers_file)
    Sequence = paste(read_Rprimers_file)
    return(data.frame(Primer.Name, Sequence))
  })
  
  import_Fprimers <- eventReactive(input$loadprimers, {
    # vF <- read.delim(paste0(primersDesign_wd,"/",input$name,"/","Fprimers.fasta"))
    # print(vF)
    # return(vF)
    
    # do this also using BioStrings, to display table nicely
    file_path_previous_run <- paste0(primersDesign_wd,"/",input$name,"/","Fprimers.fasta")
    read_primer_file <- readDNAStringSet(as.character(file_path_previous_run))
    FPrimer.Name = names(read_primer_file)
    Sequence = paste(read_primer_file)
    return(data.frame(FPrimer.Name, Sequence))
  })
  
  #user can import primers from previous step to display them
  import_Rprimers <- eventReactive(input$loadprimers, {
    # vR <- read.delim(paste0(primersDesign_wd,"/",input$name,"/","Rprimers.fasta"))
    # print(vR)
    # return(vR)
    # do this also using BioStrings, to display table nicely
    file_path_previous_run <- paste0(primersDesign_wd,"/",input$name,"/","Rprimers.fasta")
    read_primer_file <- readDNAStringSet(as.character(file_path_previous_run))
    RPrimer.Name = names(read_primer_file)
    Sequence = paste(read_primer_file)
    return(data.frame(RPrimer.Name, Sequence))
  })
  
  output$forward.primers <- DT::renderDataTable({
    if(is.null(input$Fprimers)){
      return(import_Fprimers())
    }
    else{
      return(FP())
    }
  },
  extensions = 'FixedHeader',
  options = list(fixedHeader = TRUE,
                 scrollY = TRUE),
  fillContainer = T,
  class = "display"
  )
  
  output$reverse.primers <- DT::renderDataTable({
    if(is.null(input$Rprimers)){
      return(import_Rprimers())
    }
    else{
      return(RP())
    }
  },
  extensions = 'FixedHeader',
  options = list(fixedHeader = TRUE,
                 scrollY = TRUE),
  fillContainer = T,
  class = "display"
  )
  
  output$primer_qc_start <- eventReactive(input$computePQC, {
    # inform the user that the virtual PCR has started
    showModal(modalDialog(
      title = "ePCR has started!",
      paste0("The ePCR for your primers is being computed. Your results will be available in a few moments. You can find them in the 'Results of ePCR' tab when they are ready."),
      easyClose = FALSE,
      footer = modalButton("Close")))
    
    # make result folder to hold results of this run
    result_folder <- paste(primersDesign_wd, "/ePCR/", input$blast_id, sep="")
    dir.create(result_folder)
    # make a logfile for this run
    logfile_path <- file.path(primersDesign_wd, "ePCR", input$blast_id, "logfile.txt", fsep=.Platform$file.sep)
    file.create(logfile_path)
    logfile <- file(logfile_path, open="wt")
    
    writeLines("ePCR has started!", logfile)
    writeLines(paste0("Resultfolder is: ", as.character(result_folder)), logfile)
    writeLines(paste0("Logfile is: ", logfile_path), logfile)
    writeLines(paste0("Analysis started on ", Sys.Date(), "at ", Sys.time()), logfile)
    
    # initialize graphs directory
    graphs_path <- file.path(primersDesign_wd, "ePCR", input$blast_id, "graphs")
    dir.create(graphs_path)
    writeLines(paste0("Graph directory has been made: ", graphs_path), logfile)
    
    #initialize summary file
    summary_file_path <- paste(getwd(), "/ePCR/", input$blast_id, "/Summary.txt", sep="")
    file.create(summary_file_path)
    
    summary_file <- file(summary_file_path, open="wt")
    writeLines(paste0("Summary file was created and opened for writing!"), logfile)
    
    # first check if the user has uploaded primers: check if the upload field is empty
    exists_upload = !is.null(input$Fprimers) && !is.null(input$Rprimers)
    # was there was a previous primer design job from which we can import the generated primers?
    exists_previous_job = file.exists(file.path(primersDesign_wd, input$name, "Fprimers.fasta", fsep=.Platform$file.sep)) && file.exists(file.path(primersDesign_wd, input$name, "Rprimers.fasta", fsep=.Platform$file.sep))
    # inform the user, that no primers are available
    if(!exists_upload && !exists_previous_job){
      # no file was uploaded, inform the user
      print("No file uploaded!")
      writeLines(paste0("The user did not upload an input file: "), logfile)
      writeLines(paste0("Previous job existence: ", exists_previous_job), logfile)
      writeLines(paste0("Upload existence: ", exists_upload), logfile)
      return (sprintf("No file was uploaded, please provide primers as input!"))
    }
    
    #get the sequences for the forward and reverse primers to be used == Fseq/Rseq  
    Fseq <- (if(is.null(input$Fprimers)) {
      vF <- readDNAStringSet(paste0(primersDesign_wd,"/",input$name,"/","Fprimers.fasta"))
    }
    else{
      #validate(need(readDNAStringSet(input$Fprimers$datapath)), "no upload")
      vF <- readDNAStringSet(input$Fprimers$datapath)
    } )
    Rseq <- (if(is.null(input$Rprimers)) {
      vR <- readDNAStringSet(paste0(primersDesign_wd,"/",input$name,"/","Rprimers.fasta"))
      
    }
    else{
      vR <- readDNAStringSet(input$Rprimers$datapath)
    } )
    
    # output of the used primers
    print(Fseq)
    writeLines(paste0("Used forward primers: ", Fseq), logfile)
    print(Rseq)
    writeLines(paste0("Used reverse primers: ", Rseq), logfile)
    
    #check if this is a bisulfite blast
    is_bis = input$is_bisulfite
    if (is_bis){
      print ("Starting bisulfite ePCR")
      writeLines(paste0("Bisulfite ePCR is conducted!"), logfile)
      
      # get a new reference genome to the organism in question using the ReferenceGenome class
      refgen <- new("ReferenceGenome", genome=getBSgenome(paste0("BSgenome.", gsub("\\.", "\\.UCSC\\.", input$genome))), name=paste0("BSgenome.", gsub("\\.", "\\.UCSC\\.", input$genome)), wd=file.path(primersDesign_wd, "database", paste0("BSgenome.", gsub("\\.", "\\.UCSC\\.", input$genome)), fsep=.Platform$file.sep))
      writeLines(paste0("Generated new reference genome: ", refgen@name), logfile)
      # build/get CT and GA databases to the organism in question --> this is handles by the ReferenceGenome.R class, that needs to be available
      dbList <- getBlastDB(refgen, input$is_bisulfite)
      writeLines(paste0("Database for new reference genome was fetched!"), logfile)
      
      # arguments for conducting the blast
      blast_args <- "-task blastn -evalue %s"
      # the e-value for calculation of bisulfite primers is set lower than the commonly used value of 10, to avoid too many (unnecessary) results
      costumized_BLAST_args <- sprintf(blast_args, 5)
      #print(costumized_BLAST_args)
      writeLines(paste0("performed blast using the arguments: ", costumized_blast_args), logfile)
      
      #blast forward and reverse primer against CT and GA converted genome
      F_CTblast <- predict(dbList$CTdb, Fseq, BLAST_args = costumized_BLAST_args)
      R_CTblast <- predict(dbList$CTdb, Rseq, BLAST_args = costumized_BLAST_args)
      F_GAblast <- predict(dbList$GAdb, Fseq, BLAST_args = costumized_BLAST_args)
      R_GAblast <- predict(dbList$GAdb, Rseq, BLAST_args = costumized_BLAST_args)
      writeLines(paste0("Blast has been performed and per primer results were written to results folder!"), logfile)
      
      # filter out hits with too many mismatches
      F_CTblast <- subset(F_CTblast, 
                              Mismatches <= input$primer_mismatches)
      R_CTblast <- subset(R_CTblast, 
                              Mismatches <= input$primer_mismatches)
      F_GAblast <- subset(F_GAblast, 
                         Mismatches <= input$primer_mismatches)
      R_GAblast <- subset(R_GAblast, 
                          Mismatches <= input$primer_mismatches)
      
      writeLines(paste0("Filtering blast results with too many mismatches done!"), logfile)
      
      #write blast results for both primers to seperate tables
      write.table(F_CTblast, paste0(result_folder, "\\Blast_Hits_Forward_Primers_C_to_T_converted_refgen"), col.names=T,row.names=F,sep="\t",dec=".",quote=F) 
      write.table(R_CTblast, paste(result_folder, "\\Blast_Hits_Reverse_Primers_C_to_T_converted_refgen"), col.names=T,row.names=F,sep="\t",dec=".",quote=F)
      write.table(F_GAblast, paste0(result_folder, "\\Blast_Hits_Forward_Primers_G_to_A_converted_refgen"), col.names=T,row.names=F,sep="\t",dec=".",quote=F) 
      write.table(R_GAblast, paste0(result_folder, "\\Blast_Hits_Reverse_Primers_G_to_A_converted_refgen"), col.names=T,row.names=F,sep="\t",dec=".",quote=F) 
      writeLines(paste0("Results of primer blasts have been written to the result tables!"), logfile)
      
      # write to summary file 
      writeLines(paste("Parameters \t Settings \n"), summary_file)
      writeLines(paste("Analysis start \t", Sys.time(), "\n"), summary_file)
      writeLines(paste("Result folder \t", result_folder, "\n"), summary_file)

      # total number of primers
      writeLines(paste("Total number of Primerpairs \t", length(Fseq), "\n"), summary_file)
      
      vec_matches_fprimers_CT <- vector()
      vec_matches_rprimers_CT <- vector()
      vec_matches_fprimers_GA <- vector()
      vec_matches_rprimers_GA <- vector()
      vec_seqnames <- vector()
      
      writeLines(paste0("Starting Statistic of Blast Matches per Primerpair!"), logfile)
      # loop over all primers and count the blast matches
      for (i in names(Fseq)){
        vec_seqnames <- c(vec_seqnames, i)
        primer_subsetCT = subset(F_CTblast, F_CTblast$QueryID == i)
        matchesCT = nrow(primer_subsetCT) 
        perfect_matchesCT = nrow(subset(primer_subsetCT, primer_subsetCT$Perc.Ident == 100))
        primer_subsetGA = subset(F_GAblast, F_GAblast$QueryID == i)
        matchesGA = nrow(primer_subsetGA) 
        perfect_matchesGA = nrow(subset(primer_subsetGA, primer_subsetGA$Perc.Ident == 100))
        vec_matches_fprimers_CT <- c(vec_matches_fprimers_CT, matchesCT)
        vec_matches_fprimers_GA <- c(vec_matches_fprimers_GA, matchesGA)
        writeLines(paste("Total number matches forward primer ", i, " to CT-reference genome \t", matchesCT,  "\n"), summary_file)
        writeLines(paste("Total number perfect matches forward primer", i, "  to CT-reference genome \t", perfect_matchesCT, "\n"), summary_file)
        writeLines(paste("Total number matches forward primer", i, " to GA-reference genome \t", matchesGA,  "\n"), summary_file)
        writeLines(paste("Total number perfect matches forward primer", i, " to GA-reference genome \t", perfect_matchesGA, "\n"), summary_file)
      }
      writeLines(paste0("Statistic of forward primers finished!"), logfile)
      
      for (i in names(Rseq)){
        primer_subsetCT = subset(R_CTblast, R_CTblast$QueryID == i)
        matchesCT = nrow(primer_subsetCT) 
        perfect_matchesCT = nrow(subset(primer_subsetCT, primer_subsetCT$Perc.Ident == 100))
        primer_subsetGA = subset(R_GAblast, R_GAblast$QueryID == i)
        matchesGA = nrow(primer_subsetGA) 
        perfect_matchesGA = nrow(subset(primer_subsetGA, primer_subsetGA$Perc.Ident == 100))
        vec_matches_rprimers_CT <- c(vec_matches_rprimers_CT, matchesCT)
        vec_matches_rprimers_GA <- c(vec_matches_rprimers_GA, matchesGA)
        writeLines(paste("Total number matches reverse primer ", i, " to CT-reference genome \t", matchesCT,  "\n"), summary_file)
        writeLines(paste("Total number perfect matches reverse primer ", i, " to CT-reference genome \t", perfect_matchesCT, "\n"), summary_file)
        writeLines(paste("Total number matches reverse primer ", i, " to GA-reference genome \t", matchesGA,  "\n"), summary_file)
        writeLines(paste("Total number perfect matches reverse primer ", i, " to GA-reference genome \t", perfect_matchesGA, "\n"), summary_file)
      }
      writeLines(paste0("Statistic of reverse primers finished!"), logfile)
      
      # create plot to visualize number of blast hits for bisulfite primers
      data_plot2 <- matrix(c(as.numeric(vec_matches_fprimers_CT), as.numeric(vec_matches_fprimers_GA), as.numeric(vec_matches_rprimers_CT), as.numeric(vec_matches_rprimers_GA)), nrow=length(Fseq), ncol=4)
      rownames(data_plot2) <- vec_seqnames
      png(file.path(graphs_path, "Blasthits_Bisulfite_Blast.png"), height=1000, width=1200, pointsize=24)
      barplot(t(data_plot2),
              main = "Blasthits per Primer to CT and GA converted genome",
              xlab = "Primer",
              col = c("black", "white", "#3c8dbc", "#f39c12")
      )
      legend("topleft",
             c("Blasthits forward primer C-to-T converted reference genome", "Blasthits forward primer G-to-A converted reference genome", "Blasthits reverse primer C-to-T converted reference genome", "Blasthits reverse primer G-to-A converted reference genome"),
             fill = c("black", "white", "#3c8dbc","#f39c12")
      )
      dev.off()
      writeLines(paste0("Created overview plot of blast hits per primerpair!"), logfile)
      
      # finding genomic ranges for all hits 
      hits<- c(
        GRanges(
          Source="F_CT", 
          AmpliconID = F_CTblast[["QueryID"]],
          seqnames = F_CTblast[["SubjectID"]],
          ranges = IRanges(
            start=pmin(F_CTblast[["S.start"]], F_CTblast[["S.end"]]),
            end = pmax(F_CTblast[["S.start"]], F_CTblast[["S.end"]])
          ), 
          strand=ifelse(F_CTblast[["S.start"]]>F_CTblast[["S.end"]],"-","+"), 
          length=F_CTblast[["Alignment.Length"]], 
          mismatches=F_CTblast[["Mismatches"]], 
          bit_score=F_CTblast[["Bits"]], 
          e_value=F_CTblast[["E"]]
        ),
        GRanges(
          Source="F_GA",
          AmpliconID = F_GAblast[["QueryID"]],
          seqnames = F_GAblast[["SubjectID"]],
          ranges = IRanges(
            start=pmin(F_GAblast[["S.start"]], F_GAblast[["S.end"]]),
            end=pmax(F_GAblast[["S.start"]], F_GAblast[["S.end"]])
          ),
          strand=ifelse(F_GAblast[["S.start"]]>F_GAblast[["S.end"]],"+","-"), 
          length=F_GAblast[["Alignment.Length"]], 
          mismatches=F_GAblast[["Mismatches"]],
          bit_score=F_GAblast[["Bits"]],
          e_value=F_GAblast[["E"]]
        ),
        GRanges(
          Source="R_CT",
          AmpliconID = R_CTblast[["QueryID"]],
          seqnames = R_CTblast[["SubjectID"]],
          ranges = IRanges(
            start=pmin(R_CTblast[["S.start"]], R_CTblast[["S.end"]]),
            end = pmax(R_CTblast[["S.start"]], R_CTblast[["S.end"]])
          ),
          strand=ifelse(R_CTblast[["S.start"]]>R_CTblast[["S.end"]],"-","+"),
          length=R_CTblast[["Alignment.Length"]],
          mismatches=R_CTblast[["Mismatches"]],
          bit_score=R_CTblast[["Bits"]],
          e_value=R_CTblast[["E"]]
        ),
        GRanges(
          Source="R_GA",
          AmpliconID = R_GAblast[["QueryID"]],
          seqnames = R_GAblast[["SubjectID"]],
          ranges = IRanges(
            start=pmin(R_GAblast[["S.start"]], R_GAblast[["S.end"]]),
            end = pmax(R_GAblast[["S.start"]], R_GAblast[["S.end"]])
          ),
          strand=ifelse(R_GAblast[["S.start"]]>R_GAblast[["S.end"]],"+","-"),
          length=R_GAblast[["Alignment.Length"]], 
          mismatches=R_GAblast[["Mismatches"]],
          bit_score=R_GAblast[["Bits"]], 
          e_value=R_GAblast[["E"]]
        )
      )
      writeLines(paste0("Generated GRanges objects for calculating overlaps of blast hits!"), logfile)
      
      #overlapping between genomic ranges 
      overlap_hits <- findOverlaps(hits, hits, maxgap=input$gap, ignore.strand=TRUE)
      writeLines(paste0("Calculated potential overlaps!"), logfile)
      
      df1 <- cbind(as.data.frame(hits[overlap_hits@from, ]),as.data.frame(hits[overlap_hits@to, ]))
      
      colnames(df1) <- paste(rep(c("F", "R"), each=11), colnames(df1), sep=".")
      writeLines(paste0("Wrote overlaps to dataframe!!"), logfile)
      
      df1 <-subset(df1,
                    F.strand == "+" &
                      R.strand == "-" &
                      as.character(F.AmpliconID) == as.character(R.AmpliconID) & 
                      as.character(F.seqnames) == as.character(R.seqnames) & 
                      abs(pmin(F.start, F.end) - pmax(R.start, R.end)) <= input$gap)
      writeLines(paste0("Filtered overlap results according to user input!"), logfile)
      writeLines(paste0("Finished ePCR part for bisulfite Primerpairs!"), logfile)
      
      print("Finished ePCR part for bisulfite Primerpairs")
      
    } else {
      print("Starting non-bisulfite ePCR")
      writeLines(paste0("Genomic ePCR is conducted!"), logfile)
      
      # get a new reference genome to the organism in question using the ReferenceGenome class
      refgen <- new("ReferenceGenome", genome=getBSgenome(paste0("BSgenome.", gsub("\\.", "\\.UCSC\\.", input$genome))), name=paste0("BSgenome.", gsub("\\.", "\\.UCSC\\.", input$genome)), wd=file.path(primersDesign_wd, "database", paste0("BSgenome.", gsub("\\.", "\\.UCSC\\.", input$genome)), fsep=.Platform$file.sep))
      writeLines(paste0("Generated new reference genome: ", refgen@name), logfile)
      
      #building databases
      dbList <- getBlastDB(refgen, input$is_bisulfite)
      writeLines(paste0("Database for new reference genome was fetched!"), logfile)
      
      #blasting 
      blast_args <- "-task blastn -evalue %s"
      costumized_BLAST_args <- sprintf(blast_args, 10)
      print(costumized_BLAST_args)
      writeLines(paste0("performed blast using the arguments: ", costumized_BLAST_args), logfile)
      
      primer1_blast <- predict(dbList$genomeDB, Fseq, BLAST_args = costumized_BLAST_args)
      primer2_blast <- predict(dbList$genomeDB, Rseq, BLAST_args = costumized_BLAST_args)
      writeLines(paste0("Blast has been performed and per primer results were written to results folder!"), logfile)
      
      # filter out hits with too many mismatches
      primer1_blast <- subset(primer1_blast, 
                              Mismatches <= input$primer_mismatches)
      primer2_blast <- subset(primer2_blast, 
                              Mismatches <= input$primer_mismatches)
      writeLines(paste0("Filtering blast results with too many mismatches done!"), logfile)
      
      #write blast results for both primers to table
      write.table(primer1_blast, paste0(result_folder, "\\Blast_Hits_Forward_Primers", sep=""), col.names=T,row.names=F,sep="\t",dec=".",quote=F) 
      write.table(primer2_blast, paste(result_folder, "\\Blast_Hits_Reverse_Primers", sep=""), col.names=T,row.names=F,sep="\t",dec=".",quote=F) 
      writeLines(paste0("Results of primer blasts have been written to the result tables!"), logfile)
      
      #write to summary file
      writeLines(paste("Parameters \t Settings \n"), summary_file)
      writeLines(paste("Analysis start \t", Sys.time(), "\n"), summary_file)
      writeLines(paste("Result folder \t", result_folder, "\n"), summary_file)
      
      
      #calculate number of perfect and imperfect matches
      perfect_matches_primer1 = primer1_blast[primer1_blast$Perc.Ident == 100, ]
      perfect_matches_primer2 = primer2_blast[primer2_blast$Perc.Ident == 100, ]
      num_imperfect_matches_primer1 = nrow(primer1_blast)
      num_imperfect_matches_primer2 = nrow(primer2_blast)
      num_perfect_matches_primer1 = nrow(perfect_matches_primer1)
      num_perfect_matches_primer2 = nrow(perfect_matches_primer2)
      
      # write to summary
      writeLines(paste("Total number of Primerpairs \t", length(Fseq), "\n"), summary_file)
      writeLines(paste("Total number perfect matches forward primers \t", num_perfect_matches_primer1, "\n"), summary_file)
      writeLines(paste("Total number perfect matches reverse primers \t", num_perfect_matches_primer2, "\n"), summary_file)
      writeLines(paste("Total number imperfect matches forward primers \t", num_imperfect_matches_primer1, "\n"), summary_file)
      writeLines(paste("Total number imperfect matches reverse primers \t", num_imperfect_matches_primer2, "\n"), summary_file)
      
      # calculate primer blast matches according to primer pair and make also a nice barplot from these numbers
      # plot 2: shows Blast Hits per Primerpair, divided up to Forward and Reverse Primerhits
      fw_primer_vec <- vector()
      rw_primer_vec <- vector()
      num_fw_blast_hits <- vector()
      num_rw_blast_hits <- vector()
      
      writeLines(paste0("Starting Statistic of Blast Matches per Primerpair!"), logfile)
      for (i in names(Fseq)){
        # append names to fw_primer_vec
        fw_primer_vec <- c(fw_primer_vec, i)
        sub_table <- subset(primer1_blast, primer1_blast$QueryID == i)
        num_fw_blast_hits <- c(num_fw_blast_hits, nrow(sub_table))
        writeLines(paste("Total blast hits forward primer ", i,  "\t", nrow(sub_table), "\n"), summary_file)
      }
      writeLines(paste0("Statistic of forward primers finished!"), logfile)
      
      for (i in names(Rseq)){
        #rw_primer_vec <- c(rw_primer_vec, i)
        sub_table <- subset(primer2_blast, primer2_blast$QueryID == i)
        num_rw_blast_hits <- c(num_rw_blast_hits, nrow(sub_table))
        writeLines(paste("Total blast hits reverse primer ", i,  "\t", nrow(sub_table), "\n"), summary_file)
      }
      writeLines(paste0("Statistic of reverse primers finished!"), logfile)
      
      # create plot to visualize number of blast hits
      data_plot2 <- matrix(c(as.numeric(num_fw_blast_hits), as.numeric(num_rw_blast_hits)), nrow=length(fw_primer_vec), ncol=2)
      rownames(data_plot2) <- fw_primer_vec
      png(file.path(graphs_path, "Blasthits_per_Primerpair.png"), height=1000, width=1200, pointsize=24)
      barplot(t(data_plot2),
              main = "Blasthits per Primerpair",
              xlab = "Primerpair",
              col = c("#3c8dbc","#f39c12")
      )
      legend("topleft",
             c("Forward Primer Blasthits","Reverse Primer Blasthits"),
             fill = c("#3c8dbc","#f39c12")
      )
      dev.off()
      writeLines(paste0("Created overview plot of blast hits per primerpair!"), logfile)
      
      #to check for close regions, always choose the same id and the same chromosome
      #to be able to do this, create GRanges Object from the perfect hits
      # finding genomic ranges for all hits 
      hits<- c(
        GRanges(
          Source="ForwardPrimerBlast", #name of used ReferenceGenome
          AmpliconID = primer1_blast[["QueryID"]],
          seqnames = primer1_blast[["SubjectID"]],
          ranges = IRanges(
            start=pmin(primer1_blast[["S.start"]], primer1_blast[["S.end"]]),
            end = pmax(primer1_blast[["S.start"]], primer1_blast[["S.end"]])
          ), 
          strand=ifelse(primer1_blast[["S.start"]]>primer1_blast[["S.end"]],"-","+"), 
          length=primer1_blast[["Alignment.Length"]], 
          mismatches=primer1_blast[["Mismatches"]], 
          bit_score=primer1_blast[["Bits"]], 
          e_value=primer1_blast[["E"]]
        ),
        GRanges(
          Source="ReversePrimerBlast",
          AmpliconID = primer2_blast[["QueryID"]],
          seqnames = primer2_blast[["SubjectID"]],
          ranges = IRanges(
            start=pmin(primer2_blast[["S.start"]], primer2_blast[["S.end"]]),
            end = pmax(primer2_blast[["S.start"]], primer2_blast[["S.end"]])
          ),
          strand=ifelse(primer2_blast[["S.start"]]>primer2_blast[["S.end"]],"-","+"),
          length=primer2_blast[["Alignment.Length"]],
          mismatches=primer2_blast[["Mismatches"]],
          bit_score=primer2_blast[["Bits"]],
          e_value=primer2_blast[["E"]]
        )
      )
      writeLines(paste0("Generated GRanges objects for calculating overlaps of blast hits!"), logfile)

      #overlapping between genomic ranges and processing results
      overlap_hits <- findOverlaps(hits, hits, maxgap=input$gap, ignore.strand=TRUE)
      writeLines(paste0("Calculated potential overlaps!"), logfile)
      df1 <-cbind(as.data.frame(hits[overlap_hits@from,]), as.data.frame(hits[overlap_hits@to,]))
      colnames(df1) <- paste(rep(c("F", "R"), each=11), colnames(df1), sep=".")
      writeLines(paste0("Wrote overlaps to dataframe!!"), logfile)
      
      df1 <-subset(df1,
                   F.strand == "+" &
                     R.strand == "-" &
                     as.character(F.AmpliconID) == as.character(R.AmpliconID) & 
                     as.character(F.seqnames) == as.character(R.seqnames) & 
                     abs(pmin(F.start, F.end) - pmax(R.start, R.end)) <= input$gap)
      
      writeLines(paste0("Filtered overlap results according to user input!"), logfile)
      writeLines(paste0("Finished ePCR part for bisulfite Primerpairs!"), logfile)
      
    } # end of non-bisulfite ePCR
    
    # Check, if there were overlaps - if not, inform the user
    if(nrow(df1) == 0){
      writeLines(paste("Analysis end \t", Sys.time(), "\n"), summary_file)
      close(summary_file)
      writeLines(paste0("No potential Amplicons were found for your ePCR! Analysis end \t", Sys.time(), "\n"), logfile)
      close(logfile)
      
      # make primer_qc_results_all file, but empty, to avoid errors for display.
      results_path_empty <- file.path(primersDesign_wd, "ePCR", input$blast_id, "primer_qc_results_all.txt", fsep=.Platform$file.sep)
      file.create(results_path_empty)
      
      # make a settings file of the job executed
      settings_file_path <- file.path(primersDesign_wd, "ePCR", input$blast_id, "Settings.txt", fsep=.Platform$file.sep)
      file.create(settings_file_path)
      
      settings_file <- file(settings_file_path, open="wt")
      writeLines(paste("Parameters \t", "Settings \n", sep= ""), settings_file)
      writeLines(paste("AnalysisID \t", input$blast_id, "\n", sep = ""), settings_file)
      writeLines(paste("Referencegenome \t", input$genome, "\n", sep = ""), settings_file)
      writeLines(paste("Bisulfiteanalysis \t", input$is_bisulfite, "\n", sep = ""), settings_file)
      writeLines(paste("Maximum size of reported products \t", input$gap, "\n", sep = ""), settings_file)
      writeLines(paste("Number of mismatches allowed in Primerblast \t", input$primer_mismatches, "\n", sep = ""), settings_file)
      close(settings_file)
      
      showModal(modalDialog(
        title = "No potential PCR products found!",
        paste0("The Quality Control for your Primers is finished. We could not find any potential PCR fragments for your Primerpairs. Details of your results are available in the Primer Design Quality Control tab."),
        easyClose = FALSE,
        footer = modalButton("Close")))
      return("No PCR fragments found.")
    }
    
    writeLines(paste("Start preparing resulttable for output"), logfile)
    # add used ReferenceGenome and PCR Productsize
    df1$assembly <- getAssemblyName(refgen)
    writeLines(paste("Wrote assembly to resulttable"), logfile)
    df1$Productsize <- ifelse(df1$F.start <= df1$R.end, df1$R.end-df1$F.start+1, df1$F.start-df1$R.end+1)
    writeLines(paste0("Determined productsize for resulttable"), logfile)
    
    # use sequeces only for fragments longer than 50 bp 
    df1 <- subset(df1,
                  Productsize >= 50)
    writeLines(paste("Discarded fragments, that are too short to be of importance in the upcoming analyis (Only fragments larger than 50 bp are analyzed)"), logfile)
    
    #url example for UCSC sequence retireval: http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr1:100000,200000
    # initialize result dataframe
    res_table <- data.frame(F.seqnames=character(), 
                            F.start=integer(), 
                            F.end=integer(), 
                            F.width=integer(), 
                            F.strand=character(),
                            F.Source=character(),
                            F.AmpliconID=character(),
                            F.Length=integer(),
                            F.mismatches=integer(),
                            F.bit_score=numeric(),
                            F.e_value=numeric(),
                            R.seqnames=character(), 
                            R.start=integer(), 
                            R.end=integer(), 
                            R.width=integer(), 
                            R.strand=character(),
                            R.Source=character(),
                            R.AmpliconID=character(),
                            R.Length=integer(),
                            R.mismatches=integer(),
                            R.bit_score=numeric(),
                            R.e_value=numeric(),
                            assembly=character(),
                            Productsize=integer(),
                            Productsequence=character(),
                            CpGs=integer())
    
    writeLines(paste0("Start fetching sequences for the first 100 results of eacht primerjob"), logfile)
    for (i in names(Fseq)){
      # get subset of df1
      primer_subset = subset(df1, df1$F.AmpliconID == i)
      # write Number of created potential amplicons for this primer to summary_file
      #writeLines(paste ("Total number of potential amplicons for primer ", i, "\t", nrow(primer_subset)), summary_file)
      
      # now get sequences for the first 100 entries
      # keep in mind: not every entry of Fseq must be contained in df1!
      if(nrow(primer_subset) != 0){
        chr <- primer_subset$F.seqnames
        start <- ifelse(primer_subset$F.start <= primer_subset$R.end, primer_subset$F.start, primer_subset$R.end)
        end <- ifelse(!primer_subset$F.start <= primer_subset$R.end, primer_subset$F.start, primer_subset$R.end)
        assembly <- getAssemblyName(refgen)
      
        # fetch sequences only for first 100 results for each primer pair
      
        url.full<-paste("http://genome.ucsc.edu/cgi-bin/das/",assembly,"/dna?segment=",chr,":",formatC(start,format="f",digits=0),",",formatC(end,format="f",digits=0),sep="")
      
        writeLines(paste0("Start fetching the first sequence for the results of primer ", i), logfile)
        if(length(url.full) != 0){
          r <- GET(url.full[1])
          is_err <- tryCatch(
            s <- content(r),
            error = function(e){
              sequences <- append (sequences, "NotFound")
              next()
              print("Error occured when fetching DNA sequence!")
            }
          )
          if(!inherits(is_err, "error")){
            seq <- xml_find_all(s, ".//DNA")
            split1 <- strsplit(as.character(seq), ">")[[1]][2]
            split2 <- strsplit(as.character(split1), "<")[[1]][1]
            sequences <- c(as.character(gsub("[\r\n]", "", split2)))
          } else {
            # ERROR HANDLING
            print("handeled error!")
            sequences <- append (sequences, "NotFound")
            next()
          }
        }
        writeLines(paste0("First sequence for the results fetched succesfully"), logfile)
      
        # This is reduced to fetching only the first 100 sequences to avoid too much runtime
        print(length(url.full))
        end_vec = length(url.full)
        if (end_vec > 100){
          end_vec = 100
          writeLines(paste0("Over 100 potential amplicons were present, set end for fetching sequences to: ", end_vec), logfile)
        }
      
        writeLines(paste0("Fetching the other sequences for created amplicons"), logfile)
        if(2 < end_vec){
            for (i in 2:end_vec){ 
              #r <- GET(url.full[i]) 
              # filter errors from content function
              is_err <- tryCatch(
                r <- GET(url.full[i]),
                s <- content(r),
                error = function(e){
                  print("Error occured when fetching DNA sequence!")
                }
              )
              # handle occurance of errors from content function
              if(!inherits(is_err, "error")){
                #REAL WORK
                seq <- xml_find_all(s, ".//DNA")
                split1 <- strsplit(as.character(seq), ">")[[1]][2]
                split2 <- strsplit(as.character(split1), "<")[[1]][1]
                sequences <- append(sequences, as.character(gsub("[\r\n]", "", split2)), i-1)
              } else {
                # ERROR HANDLING
                sequences <- append (sequences, "NotFound")
                next()
              }
            }
        }
        writeLines(paste0("Fetched sequences for amplicons succesfully"), logfile)
      
        # before adding to the dataframe, make sure, vector and dataframe are of the same length
        if (end_vec == 100){
          primer_subset$Productsequence <- sequences
          writeLines(paste0("Vector length of sequences is of size 100"), logfile)
        } else {
          # append enough NAs to be able to merge vector and dataframe
          #print(nrow(primer_subset))
          #print(100-as.numeric(nrow(primer_subset)))
          #seq_vec <- append(sequences, rep(NA, times = (100-as.numeric(nrow(primer_subset)))), after = 100)
          #print_seq_vec
          print(sequences)
          print(nrow(sequences))
          print(nrow(primer_subset))
          primer_subset$Productsequence <- sequences
          writeLines(paste0("Adjusted vector length for sequences to size 100"), logfile)
        }
        # now count CpGs in Sequences by couting them
        primer_subset$CpGs <- str_count(primer_subset$Productsequence, pattern = "cg")
        writeLines(paste0("Calculated CpG content of sequences"), logfile)
  
        # put these subsequence tables into a single resulttable.
        res_table = rbind(res_table, primer_subset)
        writeLines(paste0("Bound result for primertype ", i ," to total results"), logfile)
        }
      }

    # write Number of created amplicons to summary_file
    writeLines(paste ("Total number of potential amplicons \t", nrow(df1)), summary_file)
    
    writeLines(paste("Analysis end \t", Sys.time(), "\n"), summary_file)
    close(summary_file)
    
    # make an overview file of the job executed
    overview_file_path <- file.path(primersDesign_wd, "ePCR", input$blast_id, "Overview.txt", fsep=.Platform$file.sep)
    file.create(overview_file_path)
    overview_file <- file(overview_file_path, open="wt")
    F.PrimerID <- df1[, "F.AmpliconID"]
    table_df1_FprimerOccurances <- table(F.PrimerID)
    write.table(table_df1_FprimerOccurances, overview_file, row.names = FALSE)
    close(overview_file)
    writeLines(paste0("Overview file was created succesfully"), logfile)
    
    # make a settings file of the job executed
    settings_file_path <- file.path(primersDesign_wd, "ePCR", input$blast_id, "Settings.txt", fsep=.Platform$file.sep)
    file.create(settings_file_path)
    settings_file <- file(settings_file_path, open="wt")
    writeLines(paste("Parameters \t", "Settings \n", sep= ""), settings_file)
    writeLines(paste("AnalysisID \t", input$blast_id, "\n", sep = ""), settings_file)
    writeLines(paste("Referencegenome \t", input$genome, "\n", sep = ""), settings_file)
    writeLines(paste("Bisulfiteanalysis \t", input$is_bisulfite, "\n", sep = ""), settings_file)
    writeLines(paste("Maximum size of reported products \t", input$gap, "\n", sep = ""), settings_file)
    writeLines(paste("Number of mismatches allowed in Primerblast \t", input$primer_mismatches, "\n", sep = ""), settings_file)
    close(settings_file)
    writeLines(paste0("Settings file created succesfully"), logfile)
    
    write.table(res_table, file = file.path(primersDesign_wd, "ePCR", input$blast_id, "primer_qc_results_all.txt", fsep=.Platform$file.sep),
                col.names = TRUE, row.names=FALSE, sep="\t", dec=".") 
    writeLines(paste0("Wrote results to results table!"), logfile)
    
    # create barplot for results and put it in the results folder using ggplot2
    # visualize overview results: potential amplicons per primerpair
    require(ggplot2)
    # plot 1: shows amplicon frequency per primerpair
    df_plot1 <- as.data.frame(table_df1_FprimerOccurances)
    png(file.path(graphs_path, "PotentialAmpliconFrequencies.png"), height=1000, width=1200, pointsize=24)
    barplot(df_plot1$Freq, main="Potential Amplicons per Primerpair", names=df_plot1$F.PrimerID, col="#3c8dbc", xlab="Primerpair")
    dev.off()
    writeLines(paste0("Plot of created amplicons was created!"), logfile)
    
    showModal(modalDialog(
      title = "Computation of your virtual PCR has finished!",
      paste0("The Quality Control for your Primers is finished. Your results are available in the Primer Design Quality Control tab."),
      easyClose = FALSE,
      footer = modalButton("Close")))
  
    writeLines(paste0("ePCR function is now finished!"), logfile)
    close(logfile)
    print("Finished virtual PCR!")
    return ("Finished virtual PCR!")
  })
  
  preparePQC <- reactive({
    if (!input$refreshPQC) {return(data.frame())}
    wd <- primersDesign_wd
    
    ePCR_table <- read.delim(file.path(wd, "ePCR", input$blast_id, "primer_qc_results_all.txt", fsep=.Platform$file.sep ))
   
    # filter for certain columns of the result, we are not interested in displaying E-value and Bitscore
    ePCR_table_sub = subset(ePCR_table, select = -c(F.bit_score, R.bit_score, F.e_value, R.e_value, F.width, R.width))
    
    return(ePCR_table_sub)
  })
  
  output$pQC.results <- renderUI({
    if (!input$refreshPQC) {return(data.frame())}
    # select variables to display by selectInput
    
    selectedRange <- input$test_select
    if (length(selectedRange) == 0){
      return (data.frame())
    }
    
    file_name = file.path(primersDesign_wd, "ePCR", input$blast_id, "primer_qc_results_all.txt", fsep=.Platform$file.sep)
    # empty = (nrow(file_name) == 0)
    # 
    # if(length(empty) == 0){
    #   return(data.frame())
    # }
    
    table <- read.delim(file.path(primersDesign_wd, "ePCR", as.character(input$blast_id), "primer_qc_results_all.txt", fsep=.Platform$file.sep))
    selTable <- subset(table, F.AmpliconID == as.character(selectedRange))
    selTable <- subset(selTable, select = -c(F.bit_score, R.bit_score, F.e_value, R.e_value, F.width, R.width))
    
    output$out <- DT::renderDataTable({selTable}, extensions = 'FixedHeader',
                                                  options = list(fixedHeader = FALSE,
                                                                #scrollY = "200px",
                                                                scrollX = TRUE),
                                                  #fillContainer = TRUE,
                                                  class = "display")
    DT::dataTableOutput("out")
  })
  
  # display selectInput / Dropdownmenue to filter Results for one primerpair analyzed
  observeEvent(input$refreshPQC, {
    removeUI(
      selector= "div:has(>> #test_select)",
      
      immediate = TRUE
    )
    insertUI(
      selector= "#refreshPQC",
      where = "afterEnd",
      ui = selectInput(inputId = "test_select", label = "Filter results by primerpair: ", choices = preparePQC()[6], selected = preparePQC()[6], multiple = FALSE)
    )
  }#, once = TRUE
  )
  
  
  
  ############# display the overview of the ePCR ###########
  
  showOverviewePCR <- reactive({ if (!input$overviewePCR) {return(NULL)}
    files <- data.frame(results=list.files(file.path(getwd(), "ePCR", input$blast_id, fsep=.Platform$file.sep),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("Overview",files[["results"]])])
    print(file_path)
    if(length(file_path) == 0){
      ww <-showModal(modalDialog(
        title = "No Overview Found For ePCR!",
        sprintf(paste0("Unfortunateley, we could not find any Overview for your ePCR Job."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      stop("No Overview found!")
    }
    print(file_path)
    read_file <- read.delim(file_path, header=TRUE, sep=" ")
    
    print(read_file)
    return(read_file)
  })
  
  output$ePCR.overview <- DT::renderDataTable({
    if (!input$overviewePCR) {return(data.frame())}
    showOverviewePCR()
  },
    extensions = 'FixedHeader',
    options = list(fixedHeader = FALSE,
                   scrollY = TRUE),
    fillContainer = TRUE,
    class = "display"
  )
  
  ############# display the settings of the ePCR ###########
  
  showSettingsePCR <- reactive({ if (!input$settingsePCR) {return(NULL)}
    files <- data.frame(results=list.files(file.path(getwd(), "ePCR", input$blast_id, fsep=.Platform$file.sep), full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("Settings",files[["results"]])])
    print(file_path)
    if(length(file_path) == 0){
      ww <-showModal(modalDialog(
        title = "No Settings Found For ePCR!",
        sprintf(paste0("Unfortunateley, we could not find any Settings for your ePCR Job."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      stop("No Settings found!")
    }
    print(file_path)
    read_file <- read.delim(file_path, header=TRUE, sep="\t")
    
    print(read_file)
    return(read_file)
  })
  
  output$ePCR.settings <- DT::renderDataTable({
    if (!input$settingsePCR) {return(data.frame())}
    showSettingsePCR()
  },
  extensions = 'FixedHeader',
  options = list(fixedHeader = FALSE,
                 scrollY = TRUE),
  fillContainer = TRUE,
  class = "display"
  )
  
  ############# display the logfile of the ePCR ###########
  
  showLogfileePCR <- reactive({ if (!input$logfileePCR) {return(NULL)}
    files <- data.frame(results=list.files(file.path(primersDesign_wd, "ePCR", input$blast_id, fsep=.Platform$file.sep), full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("logfile",files[["results"]])])
    print(file_path)
    if(length(file_path) == 0){
      ww <-showModal(modalDialog(
        title = "No Logfile Found For ePCR!",
        sprintf(paste0("Unfortunateley, we could not find any logfile for your ePCR Job."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      stop("No Logfile found!")
    }
    print(file_path)
    read_file <- read.delim(file_path, header=TRUE, sep="\t")
    
    print(read_file)
    return(read_file)
  })
  
  output$ePCR.logfile <- DT::renderDataTable({
    if (!input$logfileePCR) {return(data.frame())}
    showLogfileePCR()
  },
  extensions = 'FixedHeader',
  options = list(fixedHeader = FALSE,
                 scrollY = TRUE),
  fillContainer = TRUE,
  class = "display"
  )
  
  ############# display the summary of the ePCR ###########
  
  showSummaryePCR <- reactive({ if (!input$summaryePCR) {return(NULL)}
    files <- data.frame(results=list.files(file.path(getwd(), "ePCR", input$blast_id, fsep=.Platform$file.sep), full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("Summary",files[["results"]])])
    print(file_path)
    if(length(file_path) == 0){
      ww <-showModal(modalDialog(
        title = "No Summary Found For ePCR!",
        sprintf(paste0("Unfortunateley, we could not find any Summary for your ePCR Job."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      stop("No Summary found!")
    }
    print(file_path)
    read_file <- read.delim(file_path)
    
    print(read_file)
    return(read_file)
  })
  
  output$pQC.summary <- DT::renderDataTable({
    if (!input$summaryePCR) {return(data.frame())}
    showSummaryePCR()
  })
  
  # default values for the different primer types
  defaultPrimerSettings <- function(primer_type) {
    # order of the given parameters: 
    # low.complexity.primer.removal=TRUE,
    # remove.primers.with.n=input$i_remove.primers.with.n,
    # check4snps=check_snps,
    # check4repeats=check_repeats,
    # allow.repeats.in.primers=input$i_allow.repeats.in.primers,
    # allow.repeats.in.amplicon=input$i_allow.repeats.in.amplicon,
    # annotate.genes=FALSE,
    # annotate.cpg.islands=FALSE,
    # create.toplist = TRUE,
    # max.bins.low.complexity=input$i_max.bins,
    # primer.align.binsize=input$i_primer.align.binsize,
    # min.snps.amplicon=0,
    # max.snps.amplicon=input$i_snps.amplicon,
    # min.snps.primer1=0,
    # max.snps.primer1=input$i_snps.primer1,
    # min.snps.primer2=0,
    # max.snps.primer2=input$i_snps.primer2,
    # hp.length.min=input$i_hp.length[1],
    # hp.length.max=input$i_hp.length[2],
    # add.ngs.adaptor.f=adaF,
    # add.ngs.adaptor.r=adaR,
    # strand=input$i_strand,
    # min.length.primer = input$i_primerlength[1],
    # max.length.primer = input$i_primerlength[2],
    # min.Tm.primer = input$i_primertemp[1],
    # max.Tm.primer = input$i_primertemp[2],
    # max.Tm.difference.primer = input$i_meltdiff,
    # min.length.amplicon = input$i_lengthAmp[1],
    # max.length.amplicon = input$i_lengthAmp[2],
    # min.number.gc.amplicon = input$i_minGC,
    # min.number.cg.amplicon = input$i_minCG,
    # min.C2T.primer1 = input$i_minC2T,
    # min.G2A.primer2 = input$i_minG2A,
    # chop.size= input$i_chop.size
    switch(primer_type,
           "genomic" = return(      c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, 5, 12, 0, 10, 0, 0, 0, 0, NA, NA,  NA, NA,    NA, 18, 25, 50, 60, 3,  200, 500, 0, 0, NA, NA,  30, "genomic")),
           "bisulfite" = return(    c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, 7, 12, 0, 10, 0, 0, 0, 0, NA, NA,  NA, NA, "top", 20, 32, 48, 60, 5,  200, 400, 0, 5,  3,  3,  NA, "bisulfite")),
           "NOME" = return(         c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, 7, 12, 0, 10, 0, 0, 0, 0, NA, NA,  NA, NA, "top", 20, 32, 48, 60, 6,  200, 400, 5, 5,  3,  3,  NA, "NOME")),
           "CLEVER" = return(       c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, 5, 12, 0, 10, 0, 0, 0, 0, NA, NA,  NA, NA, "top", 18, 25, 50, 60, 3,  200, 500, 0, 5,  3,  3,  NA, "CLEVER")),
           "hp_bisulfite" = return( c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, 7, 12, 0, 10, 0, 0, 0, 0, 50, 200, NA, NA, "top", 20, 32, 48, 60, 5,  200, 400, 0, 1,  3,  3,  NA, "hp_bisulfite")),
           "hp_NOME" = return(      c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, 7, 12, 0, 10, 0, 0, 0, 0, 50, 200, NA, NA, "top", 20, 32, 48, 60, 6,  200, 400, 1, 1,  3,  3,  NA, "hp_NOME")),
           "hp_CLEVER" = return(    c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, 5, 12, 0, 10, 0, 0, 0, 0, 50, 200, NA, NA, "top", 18, 25, 50, 60, 3,  200, 500, 0, 1,  3,  3,  NA, "hp_CLEVER")),
           "CrispRCas9PCR" = return(c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, 5, 12, 0, 10, 0, 0, 0, 0, NA, NA,  NA, NA, "top", NA, NA, NA, NA, NA, 200, 500, 0, 1,  NA, NA, NA, "CrispRCas9PCR")))
  }
  
  ################ observing changes in settings for calling the Primerdesign Pipeline ####################
  
  observe({input$i_primer_type
    def_settings <<- defaultPrimerSettings(input$i_primer_type)
  })
  
  observe({input$i_remove.primers.with.n
    def_settings[2] <<- as.logical(input$i_remove.primers.with.n)
  })
  
  observe({input$i_allow.repeats.in.primers
    def_settings[5] <<- as.logical(input$i_allow.repeats.in.primers)
  })
  
  observe({input$i_allow.repeats.in.amplicon
    def_settings[6] <<- as.logical(input$i_allow.repeats.in.amplicon)
  })
  
  observe({input$i_max.bins_genomic
    if (is.null(input$i_max.bins_genomic)) {return(NULL)}
    def_settings[10] <<- input$i_max.bins_genomic
  })
  
  observe({input$i_max.bins_bis
    if (is.null(input$i_max.bins_bis)) {return(NULL)}
    def_settings[10] <<- input$i_max.bins_bis
  })
  
  observe({input$i_max.bins_hp_bis
    if (is.null(input$i_max.bins_hp_bis)) {return(NULL)}
    def_settings[10] <<- input$i_max.bins_hp_bis
  })
  
  observe({input$i_max.bins_hp_clever
    if (is.null(input$i_max.bins_hp_clever)) {return(NULL)}
    def_settings[10] <<- input$i_max.bins_hp_clever
  })
  
  observe({input$i_max.bins_hp_NOME
    if (is.null(input$i_max.bins_hp_NOME)) {return(NULL)}
    def_settings[10] <<- input$i_max.bins_hp_NOME
  })
  
  observe({input$i_max.bins_crispr
    if (is.null(input$i_max.bins_crispr)) {return(NULL)}
    def_settings[10] <<- input$i_max.bins_crispr
  })
  
  observe({input$i_primer.align.binsize
    if (is.null(input$i_primer.align.binsize)) {return(NULL)}
    def_settings[11] <<- input$i_primer.align.binsize
  })
  
  observe({input$i_snps.amplicon
    if (is.null(input$i_snps.amplicon)) {return(NULL)}
    def_settings[13] <<- input$i_snps.amplicon
  })
  
  observe({input$i_snps.primer1
    if (is.null(input$i_snps.primer1)) {return(NULL)}
    def_settings[15] <<- input$i_snps.primer1
  })
  
  observe({input$i_snps.primer2
    if (is.null(input$i_snps.primer2)) {return(NULL)}
    def_settings[17] <<- input$i_snps.primer2
  })
  
  observe({input$i_hp.length
    if (is.null(input$i_hp.length)) {return(NULL)}
    def_settings[18] <<- input$i_hp.length[1]
    def_settings[19] <<- input$i_hp.length[2]
    
  })
  
  observe({input$i_strand
    def_settings[22] <<- input$i_strand
  })
  
  observe({input$i_primerlength_genomic
    if (is.null(input$i_primerlength_genomic)) {return(NULL)}
    def_settings[23] <<- input$i_primerlength_genomic[1]
    def_settings[24] <<- input$i_primerlength_genomic[2]
  })
  
  observe({input$i_primerlength_bis
    if (is.null(input$i_primerlength_bis)) {return(NULL)}
    def_settings[23] <<- input$i_primerlength_bis[1]
    def_settings[24] <<- input$i_primerlength_bis[2]
  })
  
  observe({input$i_primerlength_NOME
    if (is.null(input$i_primerlength_NOME)) {return(NULL)}
    def_settings[23] <<- input$i_primerlength_NOME[1]
    def_settings[24] <<- input$i_primerlength_NOME[2]
  })
  
  observe({input$i_primerlength_clever
    if (is.null(input$i_primerlength_clever)) {return(NULL)}
    def_settings[23] <<- input$i_primerlength_clever[1]
    def_settings[24] <<- input$i_primerlength_clever[2]
  })
  
  observe({input$i_primerlength_hp_bis
    if (is.null(input$i_primerlength_hp_bis)) {return(NULL)}
    def_settings[23] <<- input$i_primerlength_hp_bis[1]
    def_settings[24] <<- input$i_primerlength_hp_bis[2]
  })
  
  observe({input$i_primerlength_hp_NOME
    if (is.null(input$i_primerlength_hp_NOME)) {return(NULL)}
    def_settings[23] <<- input$i_primerlength_hp_NOME[1]
    def_settings[24] <<- input$i_primerlength_hp_NOME[2]
  })
  
  observe({input$i_primerlength_hp_clever
    if (is.null(input$i_primerlength_hp_clever)) {return(NULL)}
    def_settings[23] <<- input$i_primerlength_hp_clever[1]
    def_settings[24] <<- input$i_primerlength_hp_clever[2]
  })
  
  observe({input$i_primertemp_genomic
    if (is.null(input$i_primertemp_genomic)) {return(NULL)}
    def_settings[25] <<- input$i_primertemp_genomic[1]
    def_settings[26] <<- input$i_primertemp_genomic[2]
  })
  
  observe({input$i_primertemp_bis
    if (is.null(input$i_primertemp_bis)) {return(NULL)}
    def_settings[25] <<- input$i_primertemp_bis[1]
    def_settings[26] <<- input$i_primertemp_bis[2]
  })
  
  observe({input$i_primertemp_clever
    if (is.null(input$i_primertemp_clever)) {return(NULL)}
    def_settings[25] <<- input$i_primertemp_clever[1]
    def_settings[26] <<- input$i_primertemp_clever[2]
  })
  
  observe({input$i_primertemp_NOME
    if (is.null(input$i_primertemp_NOME)) {return(NULL)}
    def_settings[25] <<- input$i_primertemp_NOME[1]
    def_settings[26] <<- input$i_primertemp_NOME[2]
  })
  
  observe({input$i_primertemp_hp_bis
    if (is.null(input$i_primertemp_hp_bis)) {return(NULL)}
    def_settings[25] <<- input$i_primertemp_hp_bis[1]
    def_settings[26] <<- input$i_primertemp_hp_bis[2]
  })
  
  observe({input$i_primertemp__hp_NOME
    if (is.null(input$i_primertemp_hp_NOME)) {return(NULL)}
    def_settings[25] <<- input$i_primertemp_hp_NOME[1]
    def_settings[26] <<- input$i_primertemp_hp_NOME[2]
  })
  
  observe({input$i_primertemp_hp_clever
    if (is.null(input$i_primertemp_hp_clever)) {return(NULL)}
    def_settings[25] <<- input$i_primertemp_hp_clever[1]
    def_settings[26] <<- input$i_primertemp_hp_clever[2]
  })
  
  observe({input$i_meltdiff_genomic
    if (is.null(input$i_meltdiff_genomic)) {return(NULL)}
    def_settings[27] <<- input$i_meltdiff_genomic
  })
  
  observe({input$i_meltdiff_bis
    if (is.null(input$i_meltdiff_bis)) {return(NULL)}
    def_settings[27] <<- input$i_meltdiff_bis
  })
  
  observe({input$i_meltdiff_NOME
    if (is.null(input$i_meltdiff_NOME)) {return(NULL)}
    def_settings[27] <<- input$i_meltdiff_NOME
  })
  
  observe({input$i_meltdiff_clever
    if (is.null(input$i_meltdiff_clever)) {return(NULL)}
    def_settings[27] <<- input$i_meltdiff_clever
  })
  
  observe({input$i_meltdiff_hp_NOME
    if (is.null(input$i_meltdiff_hp_NOME)) {return(NULL)}
    def_settings[27] <<- input$i_meltdiff_hp_NOME
  })
  
  observe({input$i_meltdiff_hp_clever
    if (is.null(input$i_meltdiff_hp_clever)) {return(NULL)}
    def_settings[27] <<- input$i_meltdiff_hp_clever
  })
  
  observe({input$i_meltdiff_hp_bis
    if (is.null(input$i_meltdiff_hp_bis)) {return(NULL)}
    def_settings[27] <<- input$i_meltdiff_hp_bis
  })
  
  observe({input$i_lengthAmp_genomic
    if (is.null(input$i_lengthAmp_genomic)) {return(NULL)}
    def_settings[28] <<- input$i_lengthAmp_genomic[1]
    def_settings[29] <<- input$i_lengthAmp_genomic[2]
  })
  
  observe({input$i_lengthAmp_bis
    if (is.null(input$i_lengthAmp_bis)) {return(NULL)}
    def_settings[28] <<- input$i_lengthAmp_bis[1]
    def_settings[29] <<- input$i_lengthAmp_bis[2]
  })
  
  observe({input$i_lengthAmp_NOME
    if (is.null(input$i_lengthAmp_NOME)) {return(NULL)}
    def_settings[28] <<- input$i_lengthAmp_NOME[1]
    def_settings[29] <<- input$i_lengthAmp_NOME[2]
  })
  
  observe({input$i_lengthAmp_clever
    if (is.null(input$i_lengthAmp_clever)) {return(NULL)}
    def_settings[28] <<- input$i_lengthAmp_clever[1]
    def_settings[29] <<- input$i_lengthAmp_clever[2]
  })
  
  observe({input$i_lengthAmp_hp_bis
    if (is.null(input$i_lengthAmp_hp_bis)) {return(NULL)}
    def_settings[28] <<- input$i_lengthAmp_hp_bis[1]
    def_settings[29] <<- input$i_lengthAmp_hp_bis[2]
  })
  
  observe({input$i_lengthAmp_hp_NOME
    if (is.null(input$i_lengthAmp_hp_NOME)) {return(NULL)}
    def_settings[28] <<- input$i_lengthAmp_hp_NOME[1]
    def_settings[29] <<- input$i_lengthAmp_hp_NOME[2]
  })
  
  observe({input$i_lengthAmp_hp_clever
    if (is.null(input$i_lengthAmp_hp_clever)) {return(NULL)}
    def_settings[28] <<- input$i_lengthAmp_hp_clever[1]
    def_settings[29] <<- input$i_lengthAmp_hp_clever[2]
  })
  
  observe({input$i_lengthAmp_crispr
    if (is.null(input$i_lengthAmp_crispr)) {return(NULL)}
    def_settings[28] <<- input$i_lengthAmp_crispr[1]
    def_settings[29] <<- input$i_lengthAmp_crispr[2]
  })
  
  observe({input$i_minGC_genomic
    if (is.null(input$i_minGC_genomic)) {return(NULL)}
    def_settings[30] <<- input$i_minGC_genomic
  })
  
  observe({input$i_minGC_bis
    if (is.null(input$i_minGC_bis)) {return(NULL)}
    def_settings[30] <<- input$i_minGC_bis
  })
  
  observe({input$i_minGC_NOME
    if (is.null(input$i_minGC_NOME)) {return(NULL)}
    def_settings[30] <<- input$i_minGC_NOME
  })
  
  observe({input$i_minGC_clever
    if (is.null(input$i_minGC_clever)) {return(NULL)}
    def_settings[30] <<- input$i_minGC_clever
  })
  
  observe({input$i_minGC_hp_NOME
    if (is.null(input$i_minGC_hp_NOME)) {return(NULL)}
    def_settings[30] <<- input$i_minGC_hp_NOME
  })
  
  observe({input$i_minGC_hp_clever
    if (is.null(input$i_minGC_hp_clever)) {return(NULL)}
    def_settings[30] <<- input$i_minGC_hp_clever
  })
  
  observe({input$i_minGC_crispr
    if (is.null(input$i_minGC_crispr)) {return(NULL)}
    def_settings[30] <<- input$i_minGC_crispr
  })
  
  observe({input$i_minCG_genomic
    if (is.null(input$i_minCG_genomic)) {return(NULL)}
    def_settings[31] <<- input$i_minCG_genomic
  })
  
  observe({input$i_minCG_bis
    if (is.null(input$i_minCG_bis)) {return(NULL)}
    def_settings[31] <<- input$i_minCG_bis
  })
  
  observe({input$i_minCG_NOME
    if (is.null(input$i_minCG_NOME)) {return(NULL)}
    def_settings[31] <<- input$i_minCG_NOME
  })
  
  observe({input$i_minCG_clever
    if (is.null(input$i_minCG_clever)) {return(NULL)}
    def_settings[31] <<- input$i_minCG_clever
  })
  
  observe({input$i_minCG_hp_bis
    if (is.null(input$i_minCG_hp_bis)) {return(NULL)}
    def_settings[31] <<- input$i_minCG_hp_bis
  })
  
  observe({input$i_minCG_hp_NOME
    if (is.null(input$i_minCG_hp_NOME)) {return(NULL)}
    def_settings[31] <<- input$i_minCG_hp_NOME
  })
  
  observe({input$i_minCG_hp_clever
    if (is.null(input$i_minCG_hp_clever)) {return(NULL)}
    def_settings[31] <<- input$i_minCG_hp_clever
  })
  
  observe({input$i_minCG_crispr
    if (is.null(input$i_minCG_crispr)) {return(NULL)}
    def_settings[31] <<- input$i_minCG_crispr
  })
  
  observe({input$i_minC2T
    if (is.null(input$i_minC2T)) {return(NULL)}
    def_settings[32] <<- input$i_minC2T
  })
  
  observe({input$i_minG2A
    if (is.null(input$i_minG2A)) {return(NULL)}
    def_settings[33] <<- input$i_minG2A
  })
  
  observe({input$i_chop.size
    if (is.null(input$i_chop.size)) {return(NULL)}
    def_settings[34] <<- input$i_chop.size
  })
  
  # observeEvent(input$adapterF,{
  #   #print("observed FAda")
  #   insertUI(
  #     selector= "#adapterF",
  #     where = "afterEnd",
  #     ui = textInput("adapterForward", "Add Forward adapter: ", "CTTTCCCTACACGACGCTCTTCCGATCT")
  #   )
  # })
  # 
  # output$ReverseAdapter <- renderUI({
  #   if (!input$adapterR) {return(NULL)}
  #   textInput("adapterReverse", "Add Reverse adapter: ", "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
  # })
  
  ######### switch tab when button is pressed ############
  # observeEvent(input$switch_to_graphs_ePCR, {
  #   textOutput("This is to be implemented")
  # 
  # })
  
  
  
  ######### display graphics for ePCR ############
  showGraphicsePCR <-eventReactive(input$graphics_ePCR, {
    files <- data.frame(graphs=list.files(file.path(primersDesign_wd, "ePCR", input$blast_id, "graphs"),full.names=TRUE, pattern =".png"))
    print(files)
    })
  
  output$ePCR_Graphs <-renderUI({
    #check, if there are graphics that can be displayed, if not inform the user 
    if(is.data.frame(showGraphicsePCR()) && nrow(showGraphicsePCR())== 0){
      return(data.frame())
    } else {
      ss <- lapply(1:nrow(showGraphicsePCR()),function(x){
        out_ui <- paste0("image", x)
        imageOutput(out_ui, height = "800px")
      })
      do.call(tagList,ss)
    }
  })
  
  observe({
    for(x in 1:nrow(showGraphicsePCR()))
    {
      local({
        my_x <- x
        out_ui <- paste0("image",my_x)
        output[[out_ui]] <- renderImage({
          list(src = paste(showGraphicsePCR()$graphs[my_x]),
               alt = "Image failed to render", style="width: 800px; height: 700px")
        }, deleteFile = FALSE)
      })
    }
  })
  
}