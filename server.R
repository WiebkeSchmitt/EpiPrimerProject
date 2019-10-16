# merged with primer qc step
#remove one panel and merge it with evalation

library(shiny)
library(DT)

#### for primer design ########
primersDesign_wd <- getwd() #this used to be hardcoded: "C:\\Users\\Wiebk\\Desktop\\epiprimer" 
#pipeline_jil <- file.path(primersDesign_wd, "primer.design.pipeline_jil_v10.3_standalone.r", fsep=.Platform$file.sep)
#source("C:\\Users\\Wiebk\\Desktop\\epiprimer\\primer.design.pipeline_jil_v10.3_standalone.R")
source("generalDesign.R")

#### for primer QC ########
library(devtools)
library(rBLAST)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome)
library(Biostrings) #for primer QC and Reads Extraction 
#source file for ReferenceGenome
source(file.path(primersDesign_wd, "ReferenceGenome.r", fsep=.Platform$file.sep))

############ for Flowcell QC ########
mypath <- ("./flowcell_package/newStruct_trimmed")
choices_info <- data.frame(path=list.files(mypath,full.names=TRUE, pattern =".csv"))
print(choices_info[["path"]])
choices_info[["nameS"]]<-sapply(strsplit(as.character(choices_info[["path"]]),"/"),function(x) x[length(x)])
print(choices_info[["nameS"]])

######### for Reads Extraction ###########
# the script needs Trim galore package 
# HTML needs Trim galore as well 
flowcell_folders <- data.frame(path=list.dirs(mypath,full.names=TRUE,recursive = FALSE))
packages_names <- sapply(flowcell_folders$path,function(x) basename(as.character(x)))

####### For Reads Alignment #########
wrapper_file <- file.path(primersDesign_wd, "wrapper.r", fsep=.Platform$file.sep)
source(wrapper_file)
library(seqinr)

library(shiny)
library(shinyBS)

####### For adding tooltips #######
#to use this, the shinyBS package needs to be installed: packages.install("shinyBS")
list("bsTooltip")
list("bsPopover")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
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
  })
  
  output$state <- eventReactive(input$action, {
    #primerDesign_wd <- setwd(primersDesign_wd)
    #print(primersDesign_wd)
    
    if (is.null(input$file)){
      sprintf("No Input file Uploaded!!")
    }
    
    else if(file.exists(paste(getwd(),input$name,sep="/"))){ #dir.exists(paste(getwd(),"Data",input$name,sep="/") 
      "Dataset already exists!"
      sprintf("Dataset %s already exists, please choose other ID!", input$name)
    }else{

      showNotification("Computation started!!",duration = 15,type="message")
      
      #get input values for Forward and Reverse adapters
      adaF <- input$adapterForward
      adaR <- input$adapterReverse
      
      #check, if adapters should be used, if not, remove default adapters
      if (!input$adapterF){
        adaF <- NULL
      }
      
      if (!input$adapterR){
        adaR <- NULL
      }
      
      #call the primer design pipeline
      primer.design.pipeline(Dataset(), 
                             #input.type = input$inputtype, 
                             path.out = paste(getwd(),input$name,sep="/"),
                             primer.type = input$i_primer_type,
                             low.complexity.primer.removal=TRUE,
                             remove.primers.with.n=input$i_remove.primers.with.n,
                             #use.target.columns=input$i_use.target.columns,
                             check4snps=input$i_check4snps,
                             check4repeats=input$i_check4repeats,
                             allow.repeats.in.primers=input$i_allow.repeats.in.primers,
                             allow.repeats.in.amplicon=input$i_allow.repeats.in.amplicon,
                             annotate.genes=input$i_annotate.genes,
                             annotate.cpg.islands=input$i_annotate.cpg.islands,
                             #create.toplist=input$i_create.toplist,
                             #create.graphics=input$i_create.graphics,
                             max.bins.low.complexity=input$i_max.bins.low.complexity,
                             primer.align.binsize=input$i_primer.align.binsize,
                             min.snps.amplicon=input$i_snps.amplicon[1],
                             max.snps.amplicon=input$i_snps.amplicon[2],
                             min.snps.primer1=input$i_snps.primer1[1],
                             max.snps.primer1=input$i_snps.primer1[2],
                             min.snps.primer2=input$i_snps.primer2[1],
                             max.snps.primer2=input$i_snps.primer2[2],
                             hp.length.min=input$i_hp.length[1],
                             hp.length.max=input$i_hp.length[2],
                             add.ngs.adaptor.f=adaF,
                             add.ngs.adaptor.r=adaR,
                             strand=input$i_strand,
                             min.length.primer = input$i_primerlength[1],
                             max.length.primer = input$i_primerlength[2],
                             min.Tm.primer = input$i_primertemp[1],
                             max.Tm.primer = input$i_primertemp[2],
                             max.Tm.difference.primer = input$i_meltdiff,
                             min.length.amplicon = input$i_lengthAmp[1],
                             max.length.amplicon = input$i_lengthAmp[2],
                             min.number.gc.amplicon = input$i_minGC,
                             min.number.cg.amplicon = input$i_minCG,
                             min.C2T.primer1 = input$i_minC2T,
                             min.G2A.primer2 = input$i_minG2A,
                             chop.size= input$i_chop.size
      )
      
      sprintf("Finished Computation!!")
      ww <-showModal(modalDialog(
        title = "Primers are READY!",
        sprintf(paste0("Check the results in the folder %s"),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      
      sprintf("Finished Computation!!")
    }})
  
  ################ download the results of Primers Design ####################
  
  output$downloadPrimer<- downloadHandler(
    filename = function() {
      paste(input$name,"zip",sep=".")
      #paste0(input$name, ".zip")
    },
    
    content <- function(con) {
      #ss <- getwd()
      ss <- paste(getwd(),input$name,sep="/")
      tmpdir <- tempdir()
      setwd(tempdir())
      filesToSave <- c(ss) #List to hold paths to your files in shiny
      #Put all file paths inside filesToSave...
      
      zz <- zip(zipfile=con, files = filesToSave,flags = "-r9X", extras = "",zip ="zip")
      #x <- zip(paste0(input$name), file.path(input$name), flags = "-r9X", extras = "",zip = Sys.getenv("R_ZIPCMD", "zip"))
      return(zz)
      
    }
  )
  
  ############# display the top list of primers ########### 
  
  
  showToplist <- reactive({ if (!input$toplist) {return(NULL)}
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("toplist",files[["results"]])])
    print(file_path)
    read_file <- read.delim(file_path)
    
    print(read_file)
    return(read_file)
  })
  
  output$viewtoplist <- DT::renderDataTable({
    if (!input$toplist) {return(data.frame())}
    showToplist()
  })
  
  
  ############# display the whole list of primers ###########
  
  showWholelist <- reactive({ if (!input$wholelist) {return(NULL)}
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("wholelist",files[["results"]])])
    print(file_path)
    read_file <- read.delim(file_path)
    
    print(read_file)
    return(read_file)
  })
  
  output$viewwholelist <- DT::renderDataTable({
    if (!input$wholelist) {return(data.frame())}
    showWholelist()
  })
  
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
    
    # write.table(showWholelist()[s,c("amplicon.id","primer1.sequence","primer2.sequence")],file=paste0(primersDesign_wd,"/",input$name,"/","SelectedPrimers", ".txt"),quote = FALSE,col.names = TRUE,row.names=TRUE, sep = "\t")
    
  })
  
  output$viewSelectlist <- DT::renderDataTable({
    
    return(showSelectlist())
  })
  
  
  
  ############# display the black list of primers ###########
  
  showBlacklist <- reactive({ if (!input$blacklist) {return(NULL)}
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("blacklist",files[["results"]])])
    print(file_path)
    read_file <- read.delim(file_path)
    
    print(read_file)
    return(read_file)
  })
  
  output$viewblacklist <- DT::renderDataTable({
    if (!input$blacklist) {return(data.frame())}
    showBlacklist()
  })
  
  ############# display the white list of primers ###########
  showWhitelist <- reactive({ if (!input$whitelist) {return(NULL)}
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("white",files[["results"]])])
    print(file_path)
    read_file <- read.delim(file_path)
    
    print(read_file)
    return(read_file)
  })
 
  output$viewwhitelist <- DT::renderDataTable({
    if (!input$whitelist) {return(data.frame())}
    showWhitelist()
  })
  
  ############# display the log file of the primer design ###########
  
  showLogfile <- reactive({ if (!input$logfile) {return(NULL)}
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("logfile",files[["results"]])])
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
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("primerdesigns.by.sequence",files[["results"]])])
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
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("summary",files[["results"]])])
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
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
    print(files)
    file_path <- as.character(files[["results"]][grep("settings",files[["results"]])])
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
  
  showGraphics <-eventReactive(input$graphics,{
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(graphs=list.files(paste(getwd(),input$name,"PrimerAnalysis","graphics",sep="/"),full.names=TRUE, pattern =".png"))
    
  })
  
  output$plot3 <-renderUI({
    ss <- lapply(1:nrow(showGraphics()),function(x){
      out_ui <- paste0("image",x)
      imageOutput(out_ui)
      print(imageOutput(out_ui))
    })
    
    do.call(tagList,ss)
    
  })
  
  
  observe({
    for(x in 1:nrow(showGraphics()))
    {
      local({
        my_x <- x
        out_ui <- paste0("image",my_x)
        print(my_x)
        print(out_ui)
        print(as.character(showGraphics()$graphs[my_x]))
        print(paste(showGraphics()$graphs[my_x]))
        output[[out_ui]] <- renderImage({
          list(src = paste(showGraphics()$graphs[my_x]),
               alt = "Image failed to render",style="width: 600px; height: 400px")
        }, deleteFile = FALSE)
        
      })
      
    }
  })
  
  
  
  ################ download the selected primers list ####################
  
  output$downloadSelectedPrimers<- downloadHandler(
    filename = function() {
      # paste("SelectedPrimers","txt",sep=".")
      paste0("SelectedPrimers", ".zip")
    },
    
    content <- function(con) {
      #ss <- getwd()
      ss <-  paste0(primersDesign_wd,"/",input$name,"/","SelectedPrimers",".txt")
      tmpdir <- tempdir()
      setwd(tempdir())
      filesToSave <- c(ss) #List to hold paths to your files in shiny
      #Put all file paths inside filesToSave...
      
      zz <- zip(zipfile=con, files = filesToSave,flags = "-r9X", extras = "",zip ="zip")
      #x <- zip(paste0(input$name), file.path(input$name), flags = "-r9X", extras = "",zip = Sys.getenv("R_ZIPCMD", "zip"))
      return(zz)
      print(getwd())
      
    }
  )
  
  ################ displaying the results of Primers QC #############
  
  FP <- reactive({
    if (is.null(input$Fprimers)){
      # User has not uploaded a file yet
      return(data.frame())
    }
    
    read_Fprimers_file <- read.table(input$Fprimers$datapath)
    return(read_Fprimers_file)
  })
  
  RP <- reactive({
    if (is.null(input$Rprimers)){
      # User has not uploaded a file yet
      return(data.frame())
    }
    
    read_Rprimers_file <- read.table(input$Rprimers$datapath)
    return(read_Rprimers_file)
  })
  
  import_Fprimers <- eventReactive(input$loadprimers, {
    #wd <- setwd(primersDesign_wd)
    vF <- read.delim(paste0(primersDesign_wd,"/",input$name,"/","Fprimers.fasta"))
    print(vF)
    return(vF)
    #vR <- readDNAStringSet(paste0(primersDesign_wd,"/",input$name,"/","Rprimers.fasta"))
  })
  
  #user can import primers from previous step to display them
  import_Rprimers <- eventReactive(input$loadprimers, {
    #wd <- setwd(primersDesign_wd)
    vR <- read.delim(paste0(primersDesign_wd,"/",input$name,"/","Rprimers.fasta"))
    print(vR)
    return(vR)
    #vR <- readDNAStringSet(paste0(primersDesign_wd,"/",input$name,"/","Rprimers.fasta"))
  })
  
  output$forward.primers <- DT::renderDataTable({
    if(is.null(input$Fprimers)){
      return(import_Fprimers())
    }
    else{
      
      return(FP())
    }
  })
  
  output$reverse.primers <- DT::renderDataTable({
    
    if(is.null(input$Rprimers)){
      return(import_Rprimers())
    }
    else{
      
      return(RP())
    }
  })
  
  primer_qc <-  eventReactive(input$computePQC, {
    
    showModal(modalDialog(
      title = "Computation has already started!",
      paste0("Primers QC being computed"),
      easyClose = FALSE,
      footer = modalButton("Close")))
      refgen <- new("ReferenceGenome", genome=getBSgenome(input$genome), name=(input$genome), wd=file.path(primersDesign_wd, "database", (input$genome), fsep=.Platform$file.sep))
      
      #building databases
      dbList <- getBlastDB(refgen)
      
      #get the sequences for the forward and reverse primers to be used == Fseq/Rseq  
      Fseq <- (if(is.null(input$Fprimers)) {
        #wd <- setwd(primersDesign_wd)
        vF <- readDNAStringSet(paste0(primersDesign_wd,"/",input$name,"/","Fprimers.fasta"))
        
      }
      else{
        vF <- readDNAStringSet(input$Fprimers$datapath)
      } )
      
      Rseq <- (if(is.null(input$Rprimers)) {
        #wd <- setwd(primersDesign_wd)
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
      costumized_BLAST_args <- sprintf(blast_args,input$Evalue)
      print(costumized_BLAST_args)
      F_CTblast <- predict(dbList$CTdb, Fseq,BLAST_args=costumized_BLAST_args)
      R_CTblast <- predict(dbList$CTdb, Rseq,BLAST_args=costumized_BLAST_args)
      F_GAblast <- predict(dbList$GAdb, Fseq,BLAST_args=costumized_BLAST_args)
      R_GAblast <- predict(dbList$GAdb, Rseq,BLAST_args=costumized_BLAST_args)
      
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
      
      sub1 <-subset(df1,
                    F.strand == "+" &
                      R.strand == "-" &
                      as.character(F.AmpliconID) == as.character(R.AmpliconID) & 
                      as.character(F.seqnames) == as.character(R.seqnames) & 
                      abs(pmin(F.start,F.end)-pmax(R.start,R.end))<input$gap &
                      as.character(vector) == refgen@name)
      
      print(sub1)
      return(sub1)
    
  }) 
  
  output$installed <- eventReactive(input$InstallGenome, {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(tobeinstalled, version = "3.8")
  })
  
  output$pQC.results <- DT::renderDataTable({
    
    primerQC_table <- subset(primer_qc(),
                             F.bit_score>=input$FbitScore &
                               R.bit_score>=input$RbitScore  &
                               F.e_value<=input$Evalue  & 
                               R.e_value<=input$Evalue  &
                               F.mismatches <= input$FMismatches &
                               R.mismatches <= input$RMismatches
                             
    )
    
  })

  # showSelectedpQC <- reactive({ if (!input$extractregions) {return(NULL)}
  # s = input$pQC.results_rows_selected
  # whole_selected_list <- primer_qc()[s,]
  # write.table(whole_selected_list,file=paste0(primersDesign_wd,"/","ssssqw","/","SelectedPrimersQC", ".txt"),quote = FALSE,col.names = TRUE,row.names=FALSE, sep = "\t")
  
  # short_selected_list <- primer_qc()[s,c("amplicon.id","primer1.sequence","primer2.sequence")]
  # Fprimers <- write.table(sapply(1:nrow(short_selected_list),function(x){
  # paste0(paste0(short_selected_list[x,"amplicon.id"]),paste0(short_selected_list[x,"primer1.sequence"]),'\n')
  # }),file=paste0(primersDesign_wd,"/","ssssqw","/","qcFprimers.fasta"),quote = FALSE,col.names = FALSE,row.names=FALSE, sep = "\t")
  
  # Rprimers <- write.table(sapply(1:nrow(short_selected_list),function(x){
  # paste0(paste0(short_selected_list[x,"amplicon.id"]),paste0(short_selected_list[x,"primer2.sequence"]),'\n')
  # }),file=paste0(primersDesign_wd,"/","ssssqw","/","qcRprimers.fasta"),quote = FALSE,col.names = FALSE,row.names=FALSE, sep = "\t")
  
  # return(whole_selected_list)
  
  # write.table(showWholelist()[s,c("amplicon.id","primer1.sequence","primer2.sequence")],file=paste0(primersDesign_wd,"/",input$name,"/","SelectedPrimers", ".txt"),quote = FALSE,col.names = TRUE,row.names=TRUE, sep = "\t")
  
  # })
  
  # output$viewSelectedpQC <- DT::renderDataTable({
  # if (!input$extractregions) {return(data.frame())}
  
  # return(showSelectedpQC())
  # })
  
  
  ################ displaying overview of flowcell ########
  
  PackageInput <- eventReactive(input$update, {
    
    read_table <- read.table(as.character(choices_info$path[input$flowcellPackage==choices_info$nameS]))
    indices<-colnames(read_table)[-1]
    keep<-rowSums(read_table[,indices] >=input$freq)>0
    
    #rest<- colSums(read_table[!keep,indices])
    filter_table <- read_table[keep,]
    #final <- rbind(filter_table,rest)
    return(filter_table)
    
    
  }, ignoreNULL = FALSE)
  
  output$view <- DT::renderDataTable({
    PackageInput()
  })
  
  # Downloadable csv of overview for selected package ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(PackageInput(), file, row.names = FALSE)
    }
  )
  
  ############## functions of Reads Extraction #############
  
  output$checkbox <- renderUI({
    indices_choices <- basename(list.dirs(paste0(mypath,"/",input$selectPackage),recursive = FALSE))
    checkboxGroupInput("checkboxgroup","Select indices", choices = indices_choices,selected=indices_choices[1])
    
    
  })
  
  # Select All/ Deselect All option
  
  observe({
    indices_choices <- basename(list.dirs(paste0(mypath,"/",input$selectPackage),recursive = FALSE))
    if(input$selectall == 0) return(NULL) 
    else if (input$selectall%%2 == 0)
    {
      updateCheckboxGroupInput(session,"checkboxgroup","Select indices",choices=indices_choices)
    }
    else
    {
      updateCheckboxGroupInput(session,"checkboxgroup","Select indices",choices=indices_choices,selected=indices_choices)
    }
  })
  #other way for select all worked properly for example but not in the script 
  # observe({
  # indices_choices <- basename(list.dirs(paste0(mypath,"/",input$select1),recursive = FALSE))
  # updateCheckboxGroupInput(session, 'checkboxgroup', choices = indices_choices,
  # selected = if (input$bar) indices_choices)
  # })
  
  ###################start operations####################
  
  #print(getwd())
  #set_dir <- setwd(primersDesign_wd)
  #print(getwd())
  
  output$state1 <- eventReactive(input$excute,{
    
    # if (is.null(input$primers.file)){
    # showModal(modalDialog(
    # title = "No Primers file Uploaded!",
    # paste0("Please upload your Primers file before running the tool "),
    # easyClose = FALSE,
    # footer = modalButton("Close")
    # ))
    # sprintf("No Primers file Uploaded!!")
    # }
    
    if(input$merged.files==FALSE && input$separated.files==FALSE){ #dir.exists(paste(getwd(),"Data",input$name,sep="/") 
      
      showModal(modalDialog(
        title = "No output Mode selected !",
        paste0("User should decide either merged or separated or both! "),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      #sprintf("User should decide either merged or separated or both!")
    }
    
    else if(file.exists(paste0(getwd(),"/",input$results.folderName,".zip"))){ #dir.exists(paste(getwd(),"Data",input$name,sep="/") 
      
      showModal(modalDialog(
        title = "Results Folder already exists!",
        paste0("You have a folder with the same chosen name, please choose proper name for your output folder "),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      #sprintf("Results  %s already exists, please change the output folder name!", input$results.folderName)
    }else{
      
      #give values for merge/separate options
      change_merged_output <- (if(input$merged.files==TRUE){
        vm <- "T"
      }
      else{
        vm <- "F"
      }
      )
      
      change_separated_output <- (if(input$separated.files==TRUE){
        vs <- "T"
      }
      else{
        vs <- "F"
      }
      )
      
      change_fasta_output <- (if(input$results.format=="FASTA"){
        vf <- "T"
      }
      else{
        vf <- "F"
      }
      )
      
      #showNotification("Reads Extraction started!!",duration = 15,type="message")
      
      showModal(modalDialog(
        title = "Reads Extraction started!",
        paste0("Results will be found in the folder ",input$results.folderName,'.'),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      
      #create directory to store results and required files for procedure
      #create_results_folder <- dir.create(paste0(primersDesign_wd,"/",input$results.folderName))
      print(getwd())
      
      
      #reading Linker file 
      read_linker <- (if (is.null(input$linker.file)){
        v <- "NA"
        
      }else{
        linker_read <- read.delim(input$linker.file$datapath)
        prepare_linkerFile <- write.table(linker_read,file=paste0(mypath,"/",input$results.folderName,"/", "linker.txt"),quote = FALSE,col.names = FALSE,row.names=FALSE, sep = "\t")
        #linker_file <- paste0("-d ",input$linker.file)
      })
      
      
      #reading primers file
      read_input <- (if(is.null(input$primers.file)){
        #wd <- setwd(primersDesign_wd)
        v <- read.delim(paste0(primersDesign_wd,"/",input$name,"/","PrList_ReadsExtraction.txt"))
      }
      else{
        v <- read.delim(input$primers.file$datapath)
      }
      )
      
      #adding kmers for each primer
      read_input[["miniFprimer"]] <- paste(substr(read_input[["Fprimer"]],1,5))
      read_input[["miniRprimer"]] <-  paste(substr(read_input[["Rprimer"]],1,5))
      read_input[["minicomb"]] <- paste(substr(read_input[["Fprimer"]],1,5),substr(read_input[["Rprimer"]],1,5),sep="_")
      
      print(input$primers.file)
      print(input$linker.file)
      print(input$results.format)
      print(input$results.folderName)
      print(input$min.seq)
      print(input$error.rate)
      print(input$merged.files)
      print(input$separated.files)
      print(input$selection.mode)
      print(input$selectPackage)
      
      #we need other version of primers file without column and row names in the last string
      #prepare_primerFile <- write.table(read_input[,1:3],file=paste0(mypath,"/",input$select1,"/", "primers.txt"),quote = FALSE,col.names = FALSE,row.names=FALSE, sep = "\t")
      
      # splitting the file into separated files
      split_primers <- sapply(read_input[,1],function(x){
        write.table(read_input[x,1:3],file=paste0(primersDesign_wd,"/",read_input[x,1], ".txt"),quote = FALSE,col.names = FALSE,row.names=FALSE, sep = "\t")
      })
      print(read_input[,1:3])
      
      # what to do in multiple indices????
      get_index_files <- sapply(input$checkboxgroup , function(x) data.frame(indexList=list.files(paste0(mypath,"/",input$selectPackage,"/",x),full.names=TRUE, pattern =".fastq.gz")))
      index_files <- data.frame(path=unlist(get_index_files,use.names = FALSE))
      index_files[["files"]]<-sapply(as.character(index_files[["path"]]),function(x) basename(x))
      index_files[["index"]]<-sapply(strsplit(as.character(index_files[["path"]]),"/"),function(x) x[grepl("ad1",x)==TRUE])
      index_files[["read"]]<-ifelse(grepl("_R1_", index_files[["files"]]),"R1","R2")
      index_files[["sample"]]<-sapply(strsplit(as.character(index_files[["path"]]),"/"),function(x) x[grepl("samples",x)==TRUE])
      print(head(index_files))
      
      
      
      #listing paths for all tables of indices
      tables_path <-sapply(input$checkboxgroup, function(x){
        paths <-data.frame(tables=list.files(paste0(mypath,"/",input$selectPackage,"/",x),full.names=TRUE, pattern ="counts.csv"))
      })
      print(tables_path)
      
      #building dataframe for information of all indices 
      
      df_indices <- data.frame(path=unlist(tables_path,use.names = FALSE))
      df_indices[["package"]]<-sapply(strsplit(as.character(df_indices[["path"]]),"/"),function(x) x[grepl("samples",x)==TRUE])
      df_indices[["table"]]<-sapply(strsplit(as.character(df_indices[["path"]]),"/"),function(x) x[grepl("L001_counts",x)==TRUE])#df_indices[["index_path"]]<-sapply(strsplit(as.character(df_indices[["path"]]),"/"),function(x) substr(df_indices[["path"]],1,length(x)-1))
      print(head(df_indices))
      
      #finding correspnding files for each primer by matching 5-mers of pimers and 5-mers of reads files (default: mm=1)
      reads_counts <- sapply(df_indices[["path"]],function(x){
        read_table <- as.data.frame(read.table(as.character(x)))
        list(read_table)})
      
      matching <- sapply(reads_counts,function(x){
        
        #choose selection mode
        if (input$selection.mode=="first"){
          sapply(read_input[["miniFprimer"]],function(y){
            mm <- isMatchingAt(y,as.character(x$R1),max.mismatch=input$mismatches)
            extract <- x$comb[mm==TRUE]
            
            
          }) 
          
        }
        
        else if(input$selection.mode=="second"){
          sapply(read_input[["miniRprimer"]],function(y){
            mm <- isMatchingAt(y,as.character(x$R2),max.mismatch=input$mismatches)
            extract <- x$comb[mm==TRUE]
            
          })
          
        }
        
        else if(input$selection.mode=="union"){
          
          m1 <- sapply(read_input[["miniFprimer"]],function(y){
            mm <- isMatchingAt(y,as.character(x$R1),max.mismatch=input$mismatches)
            extract1 <- x$comb[mm==TRUE]
            
          })
          
          m2 <- sapply(read_input[["miniRprimer"]],function(y){
            mm <- isMatchingAt(y,as.character(x$R2),max.mismatch=input$mismatches)
            extract1 <- x$comb[mm==TRUE]
            
          })
          
          files_union <- sapply(1:nrow(read_input),function(x){
            a <- union(m1[[x]],m2[[x]])
            list(a)
          })
          
        }
        
        else{
          m1 <- sapply(read_input[["miniFprimer"]],function(y){
            mm <- isMatchingAt(y,as.character(x$R1),max.mismatch=input$mismatches)
            extract1 <- x$comb[mm==TRUE]
            
          })
          
          m2 <- sapply(read_input[["miniRprimer"]],function(y){
            mm <- isMatchingAt(y,as.character(x$R2),max.mismatch=input$mismatches)
            extract1 <- x$comb[mm==TRUE]
            
          })
          
          files_intersection <- sapply(1:nrow(read_input),function(x){
            a <- intersect(m1[[x]],m2[[x]])
            list(a)
          })
        }
        
      })
      
      colnames(matching) <- c(input$checkboxgroup)
      print(ncol(matching))
      print(matching)
      
      allfiles_extracted <- 
        # (if(length(input$checkboxgroup) > 1){ #with more than one index we will get multiple lists
        #of list with certain number of rows and columns
        
        sapply(1:ncol(matching),function(x){
          sapply(1:nrow(matching),function(y){
            ss <- paste0("REST,",sapply(matching[,x][y],function(z) paste0(z ,collapse=",")),collapse=",")
          })
        })
      
      colnames(allfiles_extracted) <- c(input$checkboxgroup)
      
      print(allfiles_extracted)
      # else{ #with one index we will get only one list 
      
      # sapply(1:length(matching),function(x){
      # ss <- paste0("REST,",sapply(matching[x],function(z) paste0(z ,collapse=",")),collapse=",")
      
      # })
      # }
      # )
      
      #setwd(primersDesign_wd)
      
      #preparing Concatenation Strings
      
      concat_string_R1 <- sapply(1:ncol(allfiles_extracted),function(y) { #indices
        sapply(1:nrow(allfiles_extracted),function(x){ #primers
          read1 <- paste0("zcat ",mypath,"/",input$selectPackage,"/",colnames(allfiles_extracted)[y],"/","{",allfiles_extracted[x,y],"}","_R1_001.fastq.gz ","> ","primer",x,".files.", colnames(allfiles_extracted)[y],".R1_001.fastq")
          #list(read1)
        }) })
      
      concat_string_R2 <- sapply(1:ncol(allfiles_extracted),function(y) { #indices
        sapply(1:nrow(allfiles_extracted),function(x){ #primers
          read1 <- paste0("zcat ",mypath,"/",input$selectPackage,"/",colnames(allfiles_extracted)[y],"/","{",allfiles_extracted[x,y],"}","_R2_001.fastq.gz ","> ","primer",x,".files.", colnames(allfiles_extracted)[y],".R2_001.fastq")
          #list(read1)
        }) })
      
      colnames(concat_string_R1) <- c(input$checkboxgroup)
      colnames(concat_string_R2) <- c(input$checkboxgroup)
      
      print(concat_string_R1)
      print(concat_string_R2)
      
      run_concatenating_R1 <- sapply(concat_string_R1,function(x){
        system(x)})
      
      run_concatenating_R2 <- sapply(concat_string_R2,function(x){
        system(x)})
      
      #generating configure file 
      confiles <- list.files(primersDesign_wd,pattern=".fastq")
      print(confiles)
      confiquration_file <- data.frame(R1=confiles[grepl("R1_",confiles)])
      print(confiquration_file[["R1"]])
      confiquration_file[["R2"]] <- sub("R1","R2", confiquration_file[["R1"]])
      confiquration_file[["primers"]] <- sapply(strsplit(as.character(confiquration_file[["R1"]]),"[.]"),function(x) paste(x[grepl("primer",x)],"txt",sep="."))
      confiquration_file[["linker"]] <- rep(read_linker,nrow(confiquration_file))
      confiquration_file[["Output.fasta"]] <- rep(change_fasta_output,nrow(confiquration_file))
      confiquration_file[["split.files"]] <- rep(change_separated_output,nrow(confiquration_file))
      confiquration_file[["merge.files"]] <- rep(change_merged_output,nrow(confiquration_file))
      confiquration_file[["selection.mode"]] <- rep(input$selection.mode,nrow(confiquration_file))
      confiquration_file[["error.rate"]] <- rep(input$error.rate,nrow(confiquration_file))
      confiquration_file[["min.seq"]] <- rep(input$min.seq,nrow(confiquration_file))
      confiquration_file[["subfolder"]] <- sapply(strsplit(as.character(confiquration_file[["R1"]]),"[.]"),function(x) x[grepl("ad1",x)])
      
      print(confiquration_file)
      
      write_configure_file <- write.table(confiquration_file,file=paste0(primersDesign_wd,"/",input$results.folderName,".config"),quote = FALSE,col.names = FALSE,row.names=FALSE, sep = "\t")
      
      #arguments of the script
      
      configure_file <- paste0("-c ",input$results.folderName,".config")
      results_path <- paste0("-o ", primersDesign_wd)
      results_folder_name <- paste0("-n ",input$results.folderName)
      
      
      
      #Run the script to extract the reads 
      main_string <- "bash orgIllumina_config.sh %s %s %s"
      final_string <- sprintf(main_string,configure_file,results_path,results_folder_name)
      print(final_string)
      # run_concatenating_R1 <- system(concat_string_R1)
      # run_concatenating_R1 <- system(concat_string_R2)
      run_script <- system(final_string)
      
      showModal(modalDialog(
        title = "Computation Finished!",
        paste0("You can check the results in the folder ",input$results.folderName),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      
      print("Computation Finished!!")
      
    }})
  
  ################# Reads Alignemnt Functions ############
  
  ExReads <- reactive({
    if (is.null(input$reads.file)) {
      # User has not uploaded a file yet
      return(data.frame())
    }
    
    read_input_file <- readDNAStringSet(input$reads.file$datapath)
    print(read_input_file)
    return(read_input_file)
  })
  
  
  ExReference <- reactive({
    if (is.null(input$reference.file)) {
      # User has not uploaded a file yet
      return(data.frame())
    }
    
    read_input_file <- readDNAStringSet(input$reference.file$datapath)
    return(read_input_file)
  })
  
  
  
  ###### displaying Alignment Results  #############
  
  output$AlignmentResults <- DT::renderDataTable({
    
    return(gg())
  })
  
  
  gg <- eventReactive(input$excuteAlignment,{
    
    if(is.null(input$reads.file)|| is.null(input$reference.file)){ #dir.exists(paste(getwd(),"Data",input$name,sep="/") 
      
      showModal(modalDialog(
        title = "Not All Input Uploaded !",
        paste0("User should upload both Reads and reference files! "),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      #sprintf("User should decide either merged or separated or both!")
    }
    else {
      
      print(input$i_reference.sequence.id)
      
      print(input$reads.file)
      
      # primerDesign_wd <- setwd(primersDesign_wd)
      # print(primersDesign_wd)
      
      # if (is.null(input$file)){
      # sprintf("No Input file Uploaded!!")
      # }
      
      # else if(file.exists(paste(getwd(),input$name,sep="/"))){ #dir.exists(paste(getwd(),"Data",input$name,sep="/") 
      # "Dataset already exists!"
      # sprintf("Dataset %s already exists, please choose other ID!", input$name)
      # }else{
      
      showNotification("Comoutation started!!",duration = 15,type="message")
      
      process <- pairwiseSequenceAlignment.biqFormat(sequences.path=ExReads(),
                                                     reference.path=ExReference(),
                                                     mode=input$alignment.mode,
                                                     alignment.type=input$i_alignemnt.type,
                                                     # minimum.alignment.score=input$i_minimum.alignment.score,
                                                     # minimum.percentage.identity=input$i_minimum.percentage.identity,
                                                     reference.sequence.id=input$i_reference.sequence.id,
                                                     sample.id=input$i_sample.id
                                                     
      )
      
      # tw <- getwd()
      # print("got the directory")
      # setwd(tw)
      # print("switch to tw")
      
      # print("switch to tw")
      sprintf("Finished Computation!!")
      
      # ww <-showModal(modalDialog(
      # title = "Alignment is READY!",
      # sprintf(paste0("Check the results in the folder %s")),
      # easyClose = FALSE,
      # footer = modalButton("Close")
      # ))
      
      # sprintf("Finished Computation!!")
      return(process)
    }
  })
  
  ################ download the results of Reads Alignment ####################
  
  output$downloadAlignmentResults<-
    
    downloadHandler(
      filename = function() {
        paste("AlignmentResults",".txt",sep="")
        #paste0(input$name, ".zip")
      },
      
      content = function(file) {
        if(is.null(input$reads.file)|| is.null(input$reference.file)){ 
          
          showModal(modalDialog(
            title = "No results yet",
            paste0("User should upload both Reads and reference files! "),
            easyClose = FALSE,
            footer = modalButton("Close")
          ))
          
        }
        else {
          write.csv(gg(), file, row.names = FALSE)
        }
      }	
    )
  
})