## app.R ##
library(shiny)
library(shinydashboard)

## working directory  ##
primersDesign_wd <- getwd() 

## Source for Primer Design  ##
source("generalDesign.R")

## libraries for primer design  ##
library(devtools)
library(rBLAST)
library(rtracklayer)
library(BSgenome)
library(Biostrings)

## tooltips ##
library(shinyBS)
library(tippy)

## for Flowcell QC ##
mypath <- ("./flowcell_package/newStruct_trimmed")
choices_info <- data.frame(path=list.files(mypath,full.names=TRUE, pattern =".csv"))
print(choices_info[["path"]])
choices_info[["nameS"]]<-sapply(strsplit(as.character(choices_info[["path"]]),"/"),function(x) x[length(x)])
print(choices_info[["nameS"]])

## for Reads Extraction ##
# the script needs Trim galore package 
# HTML needs Trim galore as well 
flowcell_folders <- data.frame(path=list.dirs(mypath,full.names=TRUE,recursive = FALSE))
packages_names <- sapply(flowcell_folders$path,function(x) basename(as.character(x)))

## For Reads Alignment ##
wrapper_file <- file.path(primersDesign_wd, "wrapper.r", fsep=.Platform$file.sep)
source(wrapper_file)
library(seqinr)

dbHeader <- dashboardHeader(title = "EpiPrimer")

## UI using shiny dashboard ##
ui <- dashboardPage(skin = "yellow",
  dbHeader,
  dashboardSidebar(
    sidebarMenu(
      menuItem("Primer Design Start", tabName = "PrimerDesign", icon = icon("dna")),
      menuItem("Advanced Primer Settings", tabName = "AdvancedPrimerSettings", icon = icon("dashboard")),
      menuItem("Results of Primer Design", tabName = "PDresults", icon = icon("list-ol")),
      menuItem("Graphs of Primer Design", tabName = "PDgraphs", icon = icon("chart-bar")),
      menuItem("Primer Quality Control", tabName = "PrimerQC", icon = icon("check-circle")),
      menuItem("Results of Quality Control", tabName = "PrimerQCResults", icon = icon("list-ol")),
      menuItem("Imprint", tabName = "Imprint", icon = icon("paw"))
    )
  ),
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    tabItems(
      # Content of Primer Design Tab
      tabItem(tabName = "PrimerDesign",
              fluidRow(
                box(title = h2("Upload file for Primer Design"),
                    status = "primary", 
                    solidHeader = TRUE,
                    # Upload data:
                    helpText("Provided genomes are:"),
                    helpText("Human genome assemblies: hg18, hg19"),
                    helpText("Mouse genome assemblies: mm9, mm10"), 
                    helpText("If your primer design job was unsuccesful, you are recommended to check your advanced settings and to compute again."),
                    fileInput("file", "Upload regions or sequence file:", accept=c("txt", "text/plain")),
                    hr(),
                    DT::dataTableOutput("table"),
                    hr(),
                    textInput("name","Dataset name:", paste0("MyPrimerSet")),
                    bsTooltip("name", "Choose a name for your primer design folder", "top", "hover"),
                    hr(),
                    textOutput("state"),
                    tags$head(tags$style("#state{color: red;
                                         font-size: 30px;
                                         font-style: italic;
                                         }"
                                          )
                              ),
                    hr(),
                    helpText("You can view the results of your primer design job using the 'Results' and 'Graphs' tabs."),
                    helpText( h4("Download Example Files here: ")),
                    downloadButton("downloadSequenceFile", "Sequence File",
                                   style="margin-left:150px; margin-right:0px"),
                    downloadButton("downloadRegionsFile", "Regions File",
                                   style="margin-left:150px; margin-right:0px")
                    ),
                box(title = h2("Basic Primer Settings"), 
                    status = "primary",
                    solidHeader = TRUE,
                    radioButtons("i_primer_type", label = h3("Primer Type"),
                                 choices = list("Genomic"="genomic", "Bisulfite" = "bisulfite", "NOME" = "NOME", "CLEVER"="CLEVER",
                                                "Bisulfite (hairpin)" = "hp_bisulfite", "NOME (hairpin)" = "hp_NOME", "CLEVER (hairpin)"="hp_CLEVER", "CrispRCas9 Amplicon"="CrispRCas9PCR"),
                                 selected = "genomic"),
                    bsTooltip("i_primer_type", "What kind of primer do you want to create?", "left", "hover"),
                    hr(),
                    conditionalPanel(
                      condition = "input.i_primer_type != 'genomic' && input.i_primer_type != 'CrispRCas9PCR'",
                      radioButtons("i_strand", label = h3("Strand"),
                                   choices = list("Top" = "top", "Bottom" = "bottom", "Both" = "both"),
                                   selected = "top"),
                      bsTooltip("i_strand", "Choose the strand for which you want to create your primers!", "left", "hover")
                    ),
                    hr(),
                    helpText( h3("Choose some options for your primer design: ")),
                    hr(),
                    checkboxInput("i_remove.primers.with.n", label = h4("Remove primers that contain N bases"), TRUE),
                    checkboxInput("i_allow.repeats.in.primers", label = h4("Allow repeats in primers"), FALSE),
                    checkboxInput("i_allow.repeats.in.amplicon", label = h4("Allow repeats in Amplicons"), FALSE)
                  ),
                fluidRow(
                  actionButton("action", label="Compute Primers", icon("fas fa-calculator"), 
                               style="color: #fff; background-color: #3c8dbc; border-color: #337ab7; padding:25px; font-size:200%; width:1400px; margin-left:100px; margin-right:0px"),
                  bsTooltip("action", "The computation of your primers may take a few minutes, please wait until you receive a notifiication that your primers are finished.", "left")
                )
              )
                    
              ),
      tabItem(tabName = "AdvancedPrimerSettings",
              box(
                  sliderInput("i_snps.amplicon", label = h4("Number of SNPs allowed in the Amplicon"),
                              min = 0, max = 80, value = 10),
                  sliderInput("i_snps.primer1", label = h4("Number of SNPs allowed in the Forward Primer"),
                              min = 0, max = 80, value = 0),
                  sliderInput("i_snps.primer2", label = h4("Number of SNPs allowed in the Reverse Primer"),
                              min = 0, max = 80, value = 0),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'genomic'",
                    sliderInput("i_max.bins.low.complexity", label = h4("Maximum length of monomeric base stretches"),
                                min = 0, max = 10, value = 5)
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'bisulfite'",
                    sliderInput("i_max.bins.low.complexity", label = h4("Maximum length of monomeric base stretches"),
                                min = 0, max = 10, value = 7)
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'NOME'",
                    sliderInput("i_max.bins.low.complexity", label = h4("Maximum length of monomeric base stretches"),
                                min = 0, max = 10, value = 7)
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'CLEVER'",
                    sliderInput("i_max.bins.low.complexity", label = h4("Maximum length of monomeric base stretches"),
                                min = 0, max = 10, value = 5)
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'hp_bisulfite'",
                    sliderInput("i_max.bins.low.complexity", label = h4("Maximum length of monomeric base stretches"),
                                min = 0, max = 10, value = 7)
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'hp_NOME'",
                    sliderInput("i_max.bins.low.complexity", label = h4("Maximum length of monomeric base stretches"),
                                min = 0, max = 10, value = 7)
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'hp_CLEVER'",
                    sliderInput("i_max.bins.low.complexity", label = h4("Maximum length of monomeric base stretches"),
                                min = 0, max = 10, value = 5)
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'CrispRCas9PCR'",
                    sliderInput("i_max.bins.low.complexity", label = h4("Maximum length of monomeric base stretches"),
                                min = 0, max = 10, value = 5)
                  ),
                  sliderInput("i_primer.align.binsize", label = h4("Maximum length for Primer Self-Interaction"),
                              min = 0, max = 50, value = 12),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'genomic'",
                    sliderInput("i_primerlength", label = h4("Primer Length"),
                                min = 10, max = 80, value = c(18, 25))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'bisulfite'",
                    sliderInput("i_primerlength", label = h4("Primer Length"),
                                min = 10, max = 80, value = c(20, 32))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'NOME'",
                    sliderInput("i_primerlength", label = h4("Primer Length"),
                                min = 10, max = 80, value = c(20, 32))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'CLEVER'",
                    sliderInput("i_primerlength", label = h4("Primer Length"),
                                min = 10, max = 80, value = c(18, 25))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'hp_bisulfite'",
                    sliderInput("i_primerlength", label = h4("Primer Length"),
                                min = 10, max = 80, value = c(20, 32))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'hp_NOME'",
                    sliderInput("i_primerlength", label = h4("Primer Length"),
                                min = 10, max = 80, value = c(20, 32))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'hp_CLEVER'",
                    sliderInput("i_primerlength", label = h4("Primer Length"),
                                min = 10, max = 80, value = c(18, 25))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'genomic'",
                    sliderInput("i_primertemp", label = h4("Primer Melting Temperature"),
                                min = 40, max = 75, value = c(50, 60))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'bisulfite'",
                    sliderInput("i_primertemp", label = h4("Primer Melting Temperature"),
                                min = 40, max = 75, value = c(48, 60))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'NOME'",
                    sliderInput("i_primertemp", label = h4("Primer Melting Temperature"),
                                min = 40, max = 75, value = c(48, 60))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'CLEVER'",
                    sliderInput("i_primertemp", label = h4("Primer Melting Temperature"),
                                min = 40, max = 75, value = c(50, 60))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'hp_bisulfite'",
                    sliderInput("i_primertemp", label = h4("Primer Melting Temperature"),
                                min = 40, max = 75, value = c(48, 60))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'hp_NOME'",
                    sliderInput("i_primertemp", label = h4("Primer Melting Temperature"),
                                min = 40, max = 75, value = c(48, 60))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'hp_CLEVER'",
                    sliderInput("i_primertemp", label = h4("Primer Melting Temperature"),
                                min = 40, max = 75, value = c(50, 60))
                  )
                  ),
              box(
                conditionalPanel(
                  condition="input.i_primer_type == 'genomic'",
                  sliderInput("i_meltdiff", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                              min = 0, max = 10, value = 3)
                ),
                conditionalPanel(
                  condition="input.i_primer_type == 'bisulfite'",
                  sliderInput("i_meltdiff", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                              min = 0, max = 10, value = 5)
                ),
                conditionalPanel(
                  condition="input.i_primer_type == 'NOME'",
                  sliderInput("i_meltdiff", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                              min = 0, max = 10, value = 6)
                ),
                conditionalPanel(
                  condition="input.i_primer_type == 'CLEVER'",
                  sliderInput("i_meltdiff", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                              min = 0, max = 10, value = 3)
                ),
                conditionalPanel(
                  condition="input.i_primer_type == 'hp_bisulfite'",
                  sliderInput("i_meltdiff", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                              min = 0, max = 10, value = 5)
                ),
                conditionalPanel(
                  condition="input.i_primer_type == 'hp_NOME'",
                  sliderInput("i_meltdiff", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                              min = 0, max = 10, value = 6)
                ),
                conditionalPanel(
                  condition="input.i_primer_type == 'hp_CLEVER'",
                  sliderInput("i_meltdiff", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                              min = 0, max = 10, value = 3)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'genomic'",
                  sliderInput("i_lengthAmp", label = h4("Amplicon Length"),
                              min = 100, max = 800, value = c(200,500))
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'bisulfite'",
                  sliderInput("i_lengthAmp", label = h4("Amplicon Length"),
                              min = 100, max = 800, value = c(200,400))
                ),  
                conditionalPanel(
                  condition = "input.i_primer_type == 'NOME'",
                  sliderInput("i_lengthAmp", label = h4("Amplicon Length"),
                              min = 100, max = 800, value = c(200,400))
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'CLEVER'",
                  sliderInput("i_lengthAmp", label = h4("Amplicon Length"),
                              min = 100, max = 800, value = c(200,500))
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'hp_bisulfite'",
                  sliderInput("i_lengthAmp", label = h4("Amplicon Length"),
                              min = 100, max = 800, value = c(200,400))
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'hp_NOME'",
                  sliderInput("i_lengthAmp", label = h4("Amplicon Length"),
                              min = 100, max = 800, value = c(200,400))
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'hp_CLEVER'",
                  sliderInput("i_lengthAmp", label = h4("Amplicon Length"),
                              min = 100, max = 800, value = c(200,500))
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'CrispRCas9PCR'",
                  sliderInput("i_lengthAmp", label = h4("Amplicon Length"),
                              min = 100, max = 800, value = c(200,500))
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'genomic'",
                  sliderInput("i_minGC", h4("Minimum Number of GCs per amplicon"),
                              min = 0, max = 10, value = 0)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'bisulfite'",
                  sliderInput("i_minGC", h4("Minimum Number of GCs per amplicon"),
                              min = 0, max = 10, value = 0)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'NOME'",
                  sliderInput("i_minGC", h4("Minimum Number of GCs per amplicon"),
                              min = 0, max = 10, value = 5)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'CLEVER'",
                  sliderInput("i_minGC", h4("Minimum Number of GCs per amplicon"),
                              min = 0, max = 10, value = 0)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'hp_bisulfite'",
                  sliderInput("i_minGC", h4("Minimum Number of GCs per amplicon"),
                              min = 0, max = 10, value = 0)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'hp_NOME'",
                  sliderInput("i_minGC", h4("Minimum Number of GCs per amplicon"),
                              min = 0, max = 10, value = 1)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'hp_CLEVER'",
                  sliderInput("i_minGC", h4("Minimum Number of GCs per amplicon"),
                              min = 0, max = 10, value = 0)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'CrispRCas9PCR'",
                  sliderInput("i_minGC", h4("Minimum Number of GCs per amplicon"),
                              min = 0, max = 10, value = 0)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'genomic'",
                  sliderInput("i_minCG", h4("Minimum Number of CGs per amplicon"),
                              min = 0, max = 10, value = 0)
                ), 
                conditionalPanel(
                  condition = "input.i_primer_type == 'bisulfite'",
                  sliderInput("i_minCG", h4("Minimum Number of CGs per amplicon"),
                              min = 0, max = 10, value = 5)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'NOME'",
                  sliderInput("i_minCG", h4("Minimum Number of CGs per amplicon"),
                              min = 0, max = 10, value = 5)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'CLEVER'",
                  sliderInput("i_minCG", h4("Minimum Number of CGs per amplicon"),
                              min = 0, max = 10, value = 5)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'hp_bisulfite'",
                  sliderInput("i_minCG", h4("Minimum Number of CGs per amplicon"),
                              min = 0, max = 10, value = 1)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'hp_NOME'",
                  sliderInput("i_minCG", h4("Minimum Number of CGs per amplicon"),
                              min = 0, max = 10, value = 1)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'hp_CLEVER'",
                  sliderInput("i_minCG", h4("Minimum Number of CGs per amplicon"),
                              min = 0, max = 10, value = 1)
                ),
                conditionalPanel(
                  condition = "input.i_primer_type == 'CrispRCas9PCR'",
                  sliderInput("i_minCG", h4("Minimum Number of CGs per amplicon"),
                              min = 0, max = 10, value = 1)
                ),
                conditionalPanel(
                  condition="input.i_primer_type != 'genomic' && input.i_primer_type != 'CrispRCas9PCR'",
                  sliderInput("i_minC2T", h4("Minimum 'C' to 'T' conversions in forward primer"),
                              min = 0, max = 10, value = 3)
                ),
                  conditionalPanel(
                    condition="input.i_primer_type != 'genomic' && input.i_primer_type != 'CrispRCas9PCR'",
                    sliderInput("i_minG2A", h4("Minimum 'G' to 'A' conversions in reverse primer"),
                                min = 0, max = 10, value = 3)
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'hp_bisulfite' || input.i_primer_type == 'hp_NOME' || input.i_primer_type == 'hp_CLEVER'", 
                    sliderInput("i_hp.length", label = h4("length of one arm in the hairpin molecule"),
                                min = 0, max = 1000, value = c(50, 200))
                  ),
                  conditionalPanel(
                    condition = "input.i_primer_type == 'genomic'",
                    sliderInput("i_chop.size", label = h4("Input sequence slicing"),
                                min = 0, max = 50, value = 30)
                    
                  )
              ),
              box(title = h3("Add Adapters to my primers: "),
                  width = 12,
                checkboxInput("adapterF", label = h4("Add a specific sequence to the 5' end of the forward primer"), FALSE),
                conditionalPanel(
                  "input.adapterF == 1",
                  textInput("adapterForward", "Forward adapter: ", "TCTTTCCCTACACGACGCTCTTCCGATCT")
                ),
                checkboxInput("adapterR", label = h4("Add a specific sequence to the 5' end of the reverse primer"), FALSE),
                conditionalPanel(
                  "input.adapterR==1",
                  textInput("adapterReverse", "Reverse adapter: ", "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
                )
                )
              ),
      tabItem(tabName = "PDresults",
              fluidRow(
                  tabBox(title = "",
                        id = "PDresultsTabbox",
                        width = 12,
                         #height = "25x",
                         #width="1x",
                         tabPanel("Overview",
                                 actionButton("primerdesigns.by.sequence", label = "Overview", icon = icon("sync-alt")),
                                 hr(),
                                  DT::dataTableOutput("viewprimerdesigns")
                                 ),
                         tabPanel("Toplist",
                                  actionButton("toplist", label = "Toplist", icon = icon("sync-alt")),
                                  hr(),
                                 DT::dataTableOutput("viewtoplist")
                                 ),
                         tabPanel("Wholelist",
                                  actionButton("wholelist", label = "Wholelist", icon = icon("sync-alt")),
                                  hr(),
                                 DT::dataTableOutput("viewwholelist")
                                  ),
                        tabPanel("Whitelist",
                                 actionButton("whitelist", label = "Whitelist", icon = icon("sync-alt")),
                                  hr(),
                                  DT::dataTableOutput("viewwhitelist")
                                  ),
                        tabPanel("Blacklist",
                                  actionButton("blacklist", label = "Blacklist", icon = icon("sync-alt")),
                                 hr(),
                                  DT::dataTableOutput("viewblacklist")
                                  ),
                        tabPanel("Logfile",
                                 actionButton("logfile", label = "Logfile", icon = icon("sync-alt")),
                                 hr(),
                                 DT::dataTableOutput("viewlogfile")
                                 ),
                        tabPanel("Settings",
                                 actionButton("settings", label = "Settings", icon = icon("sync-alt")),
                                 hr(),
                                 DT::dataTableOutput("viewsettings")
                                  ),
                        tabPanel("Summary",
                                  actionButton("Summary", label = "Summary", icon = icon("sync-alt")),
                                  hr(),
                                  DT::dataTableOutput("viewSummary")
                                 )
                      )
              ),
              fluidRow(
                    box(title = "Selected List",
                        width = 12,
                        helpText("You can download primerpairs and analyze them further by adding them to your list of selected primers below."),
                        helpText("Primers added to the Selected List will automatically be downloaded to the Folder containing all details of your Primer Design run."),
                        helpText("To add primers to your Selected List, mark them on the Toplist before generating your Selected List."),
                        hr(),
                        actionButton("selectlist", label = "Generate Selected List"),
                        hr(),
                        DT::dataTableOutput("viewSelectlist"),
                        hr(),
                        downloadButton('downloadSelectedPrimers', 'Selected Primers')
                        )
              )
              ),
      tabItem(tabName = "PDgraphs",
              box(title = h2("Graphs"),
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  actionButton("graphics", label = "Graphs", icon = icon("sync-alt")),
                  uiOutput("plot3"))
              ),
      tabItem(tabName = "PrimerQC",
              box(title = h2("Upload files for Primer Quality Control"),
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 6,
                  helpText("Upload your own Primers for Quality Control or import them from the Selected List you created during Primer Design"),
                  actionButton("loadprimers", "Import Primers", icon = icon("file-import"),
                               style="margin-left:275px; margin-right:0px"),
                  bsTooltip("loadprimers", "Import primers you added to the Select List during Primer Design", "bottom", "hover"),
                  hr(),
                  textOutput("primer_qc"),
                  hr(),
                  helpText("Or Upload your own Primers (.fasta file format is needed)"),
                  fileInput("Fprimers", "Upload Forward Primers", multiple = TRUE ,accept = ".fasta"), 
                  DT::dataTableOutput("forward.primers"),
                  fileInput("Rprimers", "Upload Reverse Primers", multiple = TRUE, accept = ".fasta"), 
                  DT::dataTableOutput("reverse.primers")
                  ),
              box(title = h2("Settings for Primer Quality Control"),
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 6,
                  selectInput("genome", "Genome for Quality Control",choices=c(installed.genomes())),
                  bsTooltip("genome", "Select the genome against which you want to blast your primers!", "top", "hover"),
                  sliderInput("gap", "Maximum Fragment Size:", min = 0, max = 50000, value = 2000),
                  helpText("Please set an E-value for your quality control. This value can be used to filter your results by only returning results, that are equal or better than the E-Value."),
                  sliderInput("Evalue", "E Value:", min = 0, max = 100, value = 10),
                  h3("Parameters to filter Results"),
                  sliderInput("FMismatches", "Forward Primers Mismatches", min = 0, max = 5, value = 3),
                  sliderInput("RMismatches", "Reverse Primers Mismatches", min = 0, max = 5, value = 3),
                  sliderInput("FbitScore", "Forward Primers Bit Score", min = 0, max = 100, value = 25),
                  sliderInput("RbitScore", "Reverse Primers Bit Score", min = 0, max = 100, value = 25)
              ),
              actionButton("computePQC", "Primer Quality Control", icon("fas fa-flask"), 
                           style="color: #fff; background-color: #3c8dbc; border-color: #337ab7; padding:25px; font-size:200%; width:1400px; margin-left:75px; margin-right:0px")
              ),
      tabItem(tabName = "PrimerQCResults",
              box(title = h2("Results of Primer Quality Control"),
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  helpText("Find here the alignment of your primers to the selected Reference genome: "),
                  DT::dataTableOutput("pQC.results"),
                  helpText("tbd.")
              )
      ),
      tabItem(tabName = "Imprint",
              box(title = h2("Imprint"),
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 12,
                  helpText("This website was created by the Epigenetics department of Saarland University."),
                  tags$a(href="http://epigenetik.uni-saarland.de/en/home/", "Visit our Homepage here.")
              )
      )
              
    )
  )
)

server <- function(input, output) {
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
    #print(primersDesign_wd)
    
    if (is.null(input$file)){
      sprintf("Please upload an input file!")
    }
    
    else if(file.exists(paste(getwd(),input$name,sep="/"))){ #dir.exists(paste(getwd(),"Data",input$name,sep="/") 
      "Dataset already exists!"
      sprintf("Dataset %s already exists, please choose other ID!", input$name)
    }else{
      
      
      ww <-showModal(modalDialog(
        title = "We are currently computing your primers!",
        sprintf(paste0("You will be notified when the computation of your primers is finished. Please be patient, this might take a few minutes."),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      #showNotification("The computation of your primers has started!",duration = 25,type="message")
      
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
      
      if (input$i_allow.repeats.in.primers || input$i_allow.repeats.in.amplicon){
        check_repeats = TRUE
      } else {
        check_repeats = FALSE
      }
      
      if (input$i_snps.amplicon != 0 || input$i_snps.primer1 != 0 || input$i_snps.primer2 != 0){
        check_snps = TRUE
      } else {
        check_snps = FALSE
      }
      
      #call the primer design pipeline
      primer.design.pipeline(Dataset(), 
                             #input.type = input$inputtype, 
                             path.out = paste(getwd(),input$name,sep="/"),
                             primer.type = input$i_primer_type,
                             low.complexity.primer.removal=TRUE,
                             remove.primers.with.n=input$i_remove.primers.with.n,
                             #use.target.columns=input$i_use.target.columns,
                             check4snps=check_snps,
                             check4repeats=check_repeats,
                             allow.repeats.in.primers=input$i_allow.repeats.in.primers,
                             allow.repeats.in.amplicon=input$i_allow.repeats.in.amplicon,
                             annotate.genes=FALSE,
                             annotate.cpg.islands=FALSE,
                             #create.toplist=input$i_create.toplist,
                             #create.graphics=input$i_create.graphics,
                             max.bins.low.complexity=input$i_max.bins.low.complexity,
                             primer.align.binsize=input$i_primer.align.binsize,
                             min.snps.amplicon=0,
                             max.snps.amplicon=input$i_snps.amplicon,
                             min.snps.primer1=0,
                             max.snps.primer1=input$i_snps.primer1,
                             min.snps.primer2=0,
                             max.snps.primer2=input$i_snps.primer2,
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
        sprintf(paste0("You can find details concerning your run in the next tap and further analyze them using the created files in the folder %s"),input$name),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
      
      sprintf("Finished Computation!!")
    }})
  
  ################ download example files for sequence and region input ####################
  
  output$downloadSequenceFile <- downloadHandler(
    filename = function() {
      paste(input$name, "txt", sep=".")
    },
    content = function(file){
      file.copy("input_sequences.txt", file)
    }
  )
  
  output$downloadRegionsFile <- downloadHandler(
    filename = function() {
      paste("input", "txt", sep=".")
    },
    content = function(file){
      file.copy("input.txt", file)
    }
  )
  
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
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
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
    
    # write.table(showWholelist()[s,c("amplicon.id","primer1.sequence","primer2.sequence")],file=paste0(primersDesign_wd,"/",input$name,"/","SelectedPrimers", ".txt"),quote = FALSE,col.names = TRUE,row.names=TRUE, sep = "\t")
    
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
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
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
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
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
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
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
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
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
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
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
    #wd <- setwd(primersDesign_wd)
    #print(wd)
    files <- data.frame(results=list.files(paste(getwd(),input$name,"PrimerAnalysis",sep="/"),full.names=TRUE, pattern =".txt"))
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
  
  showGraphics <-eventReactive(input$graphics,{
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
        out_ui <- paste0("image",x)
        imageOutput(out_ui)
        print(imageOutput(out_ui))
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
    vF <- read.delim(paste0(primersDesign_wd,"/",input$name,"/","Fprimers.fasta"))
    print(vF)
    return(vF)
  })
  
  #user can import primers from previous step to display them
  import_Rprimers <- eventReactive(input$loadprimers, {
    vR <- read.delim(paste0(primersDesign_wd,"/",input$name,"/","Rprimers.fasta"))
    print(vR)
    return(vR)
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
  
  comp_started <- observeEvent(input$computePQC,{
    # inform the user that the virtual PCR has started
    showModal(modalDialog(
      title = "Computation of your virtual PCR has started!",
      paste0("The Quality Control for your Primers is being computed. Your results will be available in a few minutes. You van find them in the Primer Design Quality Control tab when they are ready."),
      easyClose = FALSE,
      footer = modalButton("Close")))
  })
  
  primer_qc <-  eventReactive(input$computePQC, {
    # get a new reference genome to the organism in question using the ReferenceGenome class
    refgen <- new("ReferenceGenome", genome=getBSgenome(input$genome), name=(input$genome), wd=file.path(primersDesign_wd, "database", (input$genome), fsep=.Platform$file.sep))
    
    #building databases
    dbList <- getBlastDB(refgen)
    
    #get the sequences for the forward and reverse primers to be used == Fseq/Rseq  
    Fseq <- (if(is.null(input$Fprimers)) {
      vF <- readDNAStringSet(paste0(primersDesign_wd,"/",input$name,"/","Fprimers.fasta"))
      
    }
    else{
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
                    abs(pmin(F.start,F.end)-pmax(R.start,R.end))<input$gap)
    
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
    
  },
  extensions = 'FixedHeader',
  options = list(fixedHeader = TRUE,
                 scrollY = TRUE),
  fillContainer = T,
  class = "display"
  )
  
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
  
  # observe({
  #   indices_choices <- basename(list.dirs(paste0(mypath,"/",input$selectPackage),recursive = FALSE))
  #   if(input$selectall == 0) return(NULL) 
  #   else if (input$selectall%%2 == 0)
  #   {
  #     updateCheckboxGroupInput(session,"checkboxgroup","Select indices",choices=indices_choices)
  #   }
  #   else
  #   {
  #     updateCheckboxGroupInput(session,"checkboxgroup","Select indices",choices=indices_choices,selected=indices_choices)
  #   }
  # })
  #other way for select all worked properly for example but not in the script 
  # observe({
  # indices_choices <- basename(list.dirs(paste0(mypath,"/",input$select1),recursive = FALSE))
  # updateCheckboxGroupInput(session, 'checkboxgroup', choices = indices_choices,
  # selected = if (input$bar) indices_choices)
  # })
  
  ###################start operations####################
  
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
  
  }

shinyApp(ui, server)