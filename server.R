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
    }
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
    s = input$viewtoplist_rows_selected
    whole_selected_list <- showToplist()[s,]
    write.table(whole_selected_list,file=paste0(primersDesign_wd,"/",input$name,"/","SelectedPrimers", ".txt"),quote = FALSE,col.names = TRUE,row.names=FALSE, sep = "\t")
    
    short_selected_list <- showToplist()[s,c("amplicon.id","primer1.sequence","primer2.sequence")]
    
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
      title = "Computation of your virtual PCR has started!",
      paste0("The Quality Control for your Primers is being computed. Your results will be available in a few moments. You can find them in the 'Results of ePCR' tab when they are ready."),
      easyClose = FALSE,
      footer = modalButton("Close")))
    
    #check if this is a bisulfite blast
    is_bis = input$is_bisulfite
    if (is_bis){
      print ("Starting bisulfite ePCR")
      
      #first see if the user has uploaded primers
      # check if upload field is empty
      exists_upload = !is.null(input$Fprimers) && !is.null(input$Rprimers)
      # was there was a previous primer design job?
      exists_previous_job = file.exists(paste0(primersDesign_wd,"/",input$name,"/","Fprimers.fasta")) && file.exists(paste0(primersDesign_wd,"/",input$name,"/","Rprimers.fasta"))
      
      if(!xor(exists_upload, exists_previous_job)){
        # no file was uploaded, inform the user
        print("No file uploaded!")
        return (sprintf("No file was uploaded, please provide primers as input!"))
      }
      
      # get a new reference genome to the organism in question using the ReferenceGenome class
      refgen <- new("ReferenceGenome", genome=getBSgenome(paste0("BSgenome.", gsub("\\.", "\\.UCSC\\.", input$genome))), name=paste0("BSgenome.", gsub("\\.", "\\.UCSC\\.", input$genome)), wd=file.path(primersDesign_wd, "database", paste0("BSgenome.", gsub("\\.", "\\.UCSC\\.", input$genome)), fsep=.Platform$file.sep))
      
      #building databases
      dbList <- getBlastDB(refgen, input$is_bisulfite)
      
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
      
      #output of the used primers
      print(Fseq)
      print(Rseq)
      
      #blasting 
      blast_args <- "-task blastn -evalue %s"
      # set e-value for calculation a little bit lower, to avoid too many (unnecessary) results
      costumized_BLAST_args <- sprintf(blast_args, 5)
      print(costumized_BLAST_args)
      
      #blast forward and reverse primer against CT and GA converted genome
      F_CTblast <- predict(dbList$CTdb, Fseq, BLAST_args = costumized_BLAST_args)
      R_CTblast <- predict(dbList$CTdb, Rseq, BLAST_args = costumized_BLAST_args)
      F_GAblast <- predict(dbList$GAdb, Fseq, BLAST_args = costumized_BLAST_args)
      R_GAblast <- predict(dbList$GAdb, Rseq, BLAST_args = costumized_BLAST_args)
      
      # filter out hits with too many mismatches
      F_CTblast <- subset(F_CTblast, 
                              Mismatches <= input$primer_mismatches)
      R_CTblast <- subset(R_CTblast, 
                              Mismatches <= input$primer_mismatches)
      F_GAblast <- subset(F_GAblast, 
                         Mismatches <= input$primer_mismatches)
      R_GAblast <- subset(R_GAblast, 
                          Mismatches <= input$primer_mismatches)
      
      # write this to table
      result_folder <- paste(primersDesign_wd, "/PrimerQC/", input$blast_id, sep="")
      dir.create(result_folder)
      #write blast results for both primers to table
      write.table(F_CTblast, paste0(result_folder, "\\Blast_Hits_Forward_Primers_C_to_T_converted_refgen", sep=""), col.names=T,row.names=F,sep="\t",dec=".",quote=F) 
      write.table(R_CTblast, paste(result_folder, "\\Blast_Hits_Reverse_Primers_C_to_T_converted_refgen", sep=""), col.names=T,row.names=F,sep="\t",dec=".",quote=F)
      write.table(F_GAblast, paste0(result_folder, "\\Blast_Hits_Forward_Primers_G_to_A_converted_refgen", sep=""), col.names=T,row.names=F,sep="\t",dec=".",quote=F) 
      write.table(R_GAblast, paste(result_folder, "\\Blast_Hits_Reverse_Primers_G_to_A_converted_refgen", sep=""), col.names=T,row.names=F,sep="\t",dec=".",quote=F) 
      
      #initialize summary file
      summary_file_path <- paste(getwd(), "/PrimerQC/", input$blast_id, "/Summary.txt", sep="")
      file.create(summary_file_path)
      summary_file <- file(summary_file_path, open="wt")
      writeLines(paste("Parameters \t Settings \n"), summary_file)
      writeLines(paste("Analysis start \t", Sys.time(), "\n"), summary_file)
      writeLines(paste("Result folder \t", result_folder, "\n"), summary_file)
      
      #calculate number of perfect and imperfect matches
      perfect_matches_primer1_C2T = F_CTblast[F_CTblast$Perc.Ident == 100, ]
      perfect_matches_primer2_C2T = R_CTblast[R_CTblast$Perc.Ident == 100, ]
      imperfect_matches_primer1_G2A = F_GAblast[F_GAblast$Perc.Ident != 100, ]
      imperfect_matches_primer2_G2A = R_GAblast[R_GAblast$Perc.Ident != 100, ]
      
      num_perfect_matches_primer1_C2T = nrow(perfect_matches_primer1_C2T)
      num_perfect_matches_primer2_C2T = nrow(perfect_matches_primer2_C2T)
      num_imperfect_matches_primer1_G2A = nrow(imperfect_matches_primer1_G2A)
      num_imperfect_matches_primer2_G2A = nrow(imperfect_matches_primer2_G2A)
      
      # write to summary
      writeLines(paste("Total number of Primerpairs \t", length(Fseq), "\n"), summary_file)
      #writeLines(paste("Total number perfect matches forward primer \t", perfect_matches_primer1_C2T, "\n"), summary_file)
      #writeLines(paste("Total number perfect matches reverse primer \t", num_perfect_matches_primer2_C2T, "\n"), summary_file)
      #writeLines(paste("Total number imperfect matches forward primer \t", num_imperfect_matches_primer1_G2A, "\n"), summary_file)
      #writeLines(paste("Total number imperfect matches reverse primer \t", num_imperfect_matches_primer2_G2A, "\n"), summary_file)
      
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
      
      #overlapping between genomic ranges 
      overlap_hits <- findOverlaps(hits,hits,maxgap=input$gap,ignore.strand=TRUE)
      
      df1 <-cbind(as.data.frame(hits[overlap_hits@from,]),as.data.frame(hits[overlap_hits@to,]))
      
      colnames(df1)<-paste(rep(c("F","R"),each=11), colnames(df1), sep=".")
      
      df1 <-subset(df1,
                    F.strand == "+" &
                      R.strand == "-" &
                      as.character(F.AmpliconID) == as.character(R.AmpliconID) & 
                      as.character(F.seqnames) == as.character(R.seqnames) & 
                      abs(pmin(F.start,F.end)-pmax(R.start,R.end)) <= input$gap)
      
      
      print("Finished ePCR part for bisulfite Primerpairs")
      
    } else {
      print("Starting non-bisulfite ePCR")
      
      #first see if the user has uploaded primers
      # check if upload field is empty
      exists_upload = !is.null(input$Fprimers) && !is.null(input$Rprimers)
      # was there was a previous primer design job?
      exists_previous_job = file.exists(paste0(primersDesign_wd,"/",input$name,"/","Fprimers.fasta")) && file.exists(paste0(primersDesign_wd,"/",input$name,"/","Rprimers.fasta"))
      
      if(!xor(exists_upload, exists_previous_job)){
        # no file was uploaded, inform the user
        print("No file uploaded!")
        return (sprintf("No file was uploaded, please provide primers as input!"))
      }
      
      # get a new reference genome to the organism in question using the ReferenceGenome class
      refgen <- new("ReferenceGenome", genome=getBSgenome(paste0("BSgenome.", gsub("\\.", "\\.UCSC\\.", input$genome))), name=paste0("BSgenome.", gsub("\\.", "\\.UCSC\\.", input$genome)), wd=file.path(primersDesign_wd, "database", paste0("BSgenome.", gsub("\\.", "\\.UCSC\\.", input$genome)), fsep=.Platform$file.sep))
      
      #building databases
      dbList <- getBlastDB(refgen, input$is_bisulfite)
      
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
      
      #blasting 
      blast_args <- "-task blastn -evalue %s"
      costumized_BLAST_args <- sprintf(blast_args, 50)
      print(costumized_BLAST_args)
      
      primer1_blast <- predict(dbList$genomeDB, Fseq, BLAST_args = costumized_BLAST_args)
      primer2_blast <- predict(dbList$genomeDB, Rseq, BLAST_args = costumized_BLAST_args)
      
      # filter out hits with too many mismatches
      primer1_blast <- subset(primer1_blast, 
                              Mismatches <= input$primer_mismatches)
      primer2_blast <- subset(primer2_blast, 
                              Mismatches <= input$primer_mismatches)
      
      # write this to table
      result_folder <- paste(primersDesign_wd, "/PrimerQC/", input$blast_id, sep="")
      print(result_folder)
      dir.create(result_folder)
      #write blast results for both primers to table
      write.table(primer1_blast, paste0(result_folder, "\\Blast_Hits_Forward_Primers", sep=""), col.names=T,row.names=F,sep="\t",dec=".",quote=F) 
      write.table(primer2_blast, paste(result_folder, "\\Blast_Hits_Reverse_Primers", sep=""), col.names=T,row.names=F,sep="\t",dec=".",quote=F) 
      
      #initialize summary file
      summary_file_path <- paste(getwd(), "/PrimerQC/", input$blast_id, "/Summary.txt", sep="")
      file.create(summary_file_path)
      summary_file <- file(summary_file_path, open="wt")
      writeLines(paste("Parameters \t Settings \n"), summary_file)
      writeLines(paste("Analysis start \t", Sys.time(), "\n"), summary_file)
      writeLines(paste("Result folder \t", result_folder, "\n"), summary_file)
      
      #calculate number of perfect and imperfect matches
      perfect_matches_primer1 = primer1_blast[primer1_blast$Perc.Ident == 100, ]
      perfect_matches_primer2 = primer2_blast[primer2_blast$Perc.Ident == 100, ]
      imperfect_matches_primer1 = primer1_blast[primer1_blast$Perc.Ident != 100, ]
      imperfect_matches_primer2 = primer2_blast[primer2_blast$Perc.Ident != 100, ]
      
      num_perfect_matches_primer1 = nrow(perfect_matches_primer1)
      num_perfect_matches_primer2 = nrow(perfect_matches_primer2)
      num_imperfect_matches_primer1 = nrow(imperfect_matches_primer1)
      num_imperfect_matches_primer2 = nrow(imperfect_matches_primer2)
      
      # write to summary
      writeLines(paste("Total number of Primerpairs \t", length(Fseq), "\n"), summary_file)
      writeLines(paste("Total number perfect matches forward primer \t", num_perfect_matches_primer1, "\n"), summary_file)
      writeLines(paste("Total number perfect matches reverse primer \t", num_perfect_matches_primer2, "\n"), summary_file)
      writeLines(paste("Total number imperfect matches forward primer \t", num_imperfect_matches_primer1, "\n"), summary_file)
      writeLines(paste("Total number imperfect matches reverse primer \t", num_imperfect_matches_primer2, "\n"), summary_file)
      
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

      #overlapping between genomic ranges and processing results
      overlap_hits <- findOverlaps(hits, hits, maxgap=input$gap, ignore.strand=TRUE)
      df1 <-cbind(as.data.frame(hits[overlap_hits@from,]), as.data.frame(hits[overlap_hits@to,]))
      colnames(df1) <- paste(rep(c("F", "R"), each=11), colnames(df1), sep=".")
      
    } # end of non-bisulfite ePCR
    
    # add used ReferenceGenome and PCR Productsize
    df1$assembly <- getAssemblyName(refgen)
    df1$Productsize <- ifelse(df1$F.start <= df1$R.end, df1$R.end-df1$F.start+1, df1$F.start-df1$R.end+1)
    
    # calculate sequeces only for fragments longer than 50 bp 
    df1 <- subset(df1,
                  Productsize >= 50)
    
    #url example for UCSC sequence retireval: http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr1:100000,200000
    chr <- df1$F.seqnames
    start <- ifelse(df1$F.start <= df1$R.end, df1$F.start, df1$R.end)
    end <- ifelse(!df1$F.start <= df1$R.end, df1$F.start, df1$R.end)
    
    assembly <- getAssemblyName(refgen)
    url.full<-paste("http://genome.ucsc.edu/cgi-bin/das/",assembly,"/dna?segment=",chr,":",formatC(start,format="f",digits=0),",",formatC(end,format="f",digits=0),sep="")
    
    if(length(url.full) != 0){
      r <- GET(url.full[1])
      is_err <- tryCatch(
        s <- content(r),
        error = function(e){
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
        sequences <- append (sequences, "NotFound")
        next()
      }
    }
    
    # This is reduced to fetching only the first 100 sequences to avoid too much runtime
    end = length(url.full)
    if (end > 100){
      end = 100
    }
    
    for (i in 2:end){  
      r <- GET(url.full[i]) 
      # filter errors from content function
      is_err <- tryCatch(
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
    
    df1$Productsequence <- sequences
    
    # now count CpGs in Sequences by couting themn 
    df1$CpGs <- str_count(df1$Productsequence, pattern = "cg")
    
    # write Number of created amplicons to summary_file
    writeLines(paste ("Total number of potential amplicons \t", nrow(df1)), summary_file)
    
    writeLines(paste("Analysis end \t", Sys.time(), "\n"), summary_file)
    close(summary_file)
    
    # make an overview file of the job executed
    overview_file_path <- paste(getwd(), "/PrimerQC/", input$blast_id, "/Overview.txt", sep="")
    file.create(overview_file_path)
    overview_file <- file(overview_file_path, open="wt")
    F.PrimerID <- df1[, "F.AmpliconID"]
    table_df1_FprimerOccurances <- table(F.PrimerID)
    write.table(table_df1_FprimerOccurances, overview_file, row.names = FALSE)
    close(overview_file)
    
    # make a settings file of the job executed
    settings_file_path <- paste(getwd(), "/PrimerQC/", input$blast_id, "/Settings.txt", sep="")
    file.create(settings_file_path)
    settings_file <- file(settings_file_path, open="wt")
    writeLines(paste("Parameters \t", "Settings \n", sep= ""), settings_file)
    writeLines(paste("AnalysisID \t", input$blast_id, "\n", sep = ""), settings_file)
    writeLines(paste("Referencegenome \t", input$genome, "\n", sep = ""), settings_file)
    writeLines(paste("Bisulfiteanalysis \t", input$is_bisulfite, "\n", sep = ""), settings_file)
    writeLines(paste("Maximum size of reported products \t", input$gap, "\n", sep = ""), settings_file)
    writeLines(paste("Number of mismatches allowed in Primerblast \t", input$primer_mismatches, "\n", sep = ""), settings_file)
    close(settings_file)
    
    # write a table instead of returning the results
    if (!dir.exists(paste(primersDesign_wd, "/PrimerQC/", input$blast_id, sep=""))){
      dir.create(paste(primersDesign_wd, "/PrimerQC/", input$blast_id, sep=""))
    }
    write.table(df1, file = paste(primersDesign_wd, "/PrimerQC/", input$blast_id, "/", "primer_qc_results_all.txt", sep=""),
                col.names = TRUE, row.names=FALSE, sep="\t", dec=".") 
    
    # create barplot for results and put it in the results folder using ggplot2
    # TODO
    #plot <- ggplot()
    
    showModal(modalDialog(
      title = "Computation of your virtual PCR has finished!",
      paste0("The Quality Control for your Primers is being finished Your results are available in the Primer Design Quality Control tab."),
      easyClose = FALSE,
      footer = modalButton("Close")))
  
    return ("Finished virtual PCR!")
  })
  
  preparePQC <- reactive({
    if (!input$refreshPQC) {return(data.frame())}
    wd <- primersDesign_wd
    primerQC_table <- read.delim(paste(wd, "/PrimerQC/", input$blast_id, "/", "primer_qc_results_all.txt", sep=""))
    
    # filter for certain columns of the result, we are not interested in displaying E-value and Bitscore
    primerQC_table_sub = subset(primerQC_table, select = -c(F.bit_score, R.bit_score, F.e_value, R.e_value, F.width, R.width))
    
    if(length(primerQC_table_sub) == 0){
      ww <-showModal(modalDialog(
        title = "No ePCR results found!",
        sprintf(paste0("Unfortunateley, we were unable to perform a Primer Quality Control for your input. Please check your settings and try again.")),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      stop("No ePCR results found!")
    }
    return(primerQC_table_sub)
  })
  
  output$pQC.results <- renderUI({
    if (!input$refreshPQC) {return(data.frame())}
    #preparePQC()
    # select variables to display by selectInput
    selectedRange <- input$test_select
    if (length(selectedRange) == 0){
      return (data.frame())
    }
    table <- read.delim(paste0(primersDesign_wd, "/PrimerQC/", as.character(input$blast_id), "/", "primer_qc_results_all.txt"), sep="")
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
    files <- data.frame(results=list.files(paste(getwd(), "PrimerQC", input$blast_id, sep="/"),full.names=TRUE, pattern =".txt"))
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
    files <- data.frame(results=list.files(paste(getwd(), "PrimerQC", input$blast_id, sep="/"),full.names=TRUE, pattern =".txt"))
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
  
  ############# display the summary of the ePCR ###########
  
  showSummaryePCR <- reactive({ if (!input$summaryePCR) {return(NULL)}
    files <- data.frame(results=list.files(paste(getwd(), "PrimerQC", input$blast_id, sep="/"),full.names=TRUE, pattern =".txt"))
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
  
}