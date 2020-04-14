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
ui <- fluidPage(dashboardPage(skin = "yellow",
                              dashboardHeader(title = "EpiPrimer"),#dbHeader,
                              dashboardSidebar(
                                sidebarMenu(
                                  menuItem("Primer Design Start", tabName = "PrimerDesign", icon = icon("dna")),
                                  menuItem("Advanced Primer Settings", tabName = "AdvancedPrimerSettings", icon = icon("dashboard")),
                                  menuItem("Results of Primer Design", tabName = "PDresults", icon = icon("list-ol")),
                                  menuItem("Graphs of Primer Design", tabName = "PDgraphs", icon = icon("chart-bar")),
                                  menuItem("Primer Quality Control", tabName = "PrimerQC", icon = icon("check-circle")),
                                  menuItem("Advanced Primer QC Settings", tabName = "PrimerQCAdvanced", icon = icon("dashboard")),
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
                                          box(title = h3("Add Adapters: "),
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
                                              helpText("Or Upload your own Primers (.fasta file format is needed, we expect your primers to be in 5' to 3' orientation)"),
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
                                              checkboxInput("is_bisulfite", h5("These are bisulfite primers!"), FALSE),
                                              helpText("Primer Quality Control is currently only available for bisulfite primers of the human genome and for genomic primers of Mouse, C. elegans and Human!"),
                                              hr(),
                                              textOutput("primer_qc_start")
                                          ),
                                          actionButton("computePQC", label = "Primer Quality Control", icon("fas fa-flask"), 
                                                       style="color: #fff; background-color: #3c8dbc; border-color: #337ab7; padding:25px; font-size:200%; width:1400px; margin-left:75px; margin-right:0px")
                                  ),
                                  tabItem(tabName = "PrimerQCAdvanced",
                                          box(title = h2("Advanced Settings for Primer Quality Control"),
                                              status = "primary",
                                              solidHeader = TRUE,
                                              width = 12,
                                              sliderInput("gap", "Maximum Fragment Size:", min = 0, max = 50000, value = 2000),
                                              helpText("Please set an E-value for your quality control. This value can be used to filter your results by only returning results, that are equal or better than the E-Value."),
                                              sliderInput("Evalue", "E Value:", min = 0, max = 100, value = 10),
                                              h3("Parameters to filter Results"),
                                              sliderInput("FMismatches", "Forward Primers Mismatches", min = 0, max = 5, value = 3),
                                              sliderInput("RMismatches", "Reverse Primers Mismatches", min = 0, max = 5, value = 3),
                                              sliderInput("FbitScore", "Forward Primers Bit Score", min = 0, max = 100, value = 25),
                                              sliderInput("RbitScore", "Reverse Primers Bit Score", min = 0, max = 100, value = 25)
                                          )
                                  ),
                                  tabItem(tabName = "PrimerQCResults",
                                          box(title = h2("Results of Primer Quality Control"),
                                              status = "primary",
                                              solidHeader = TRUE,
                                              width = 12,
                                              helpText("Find here the alignment of your primers to the selected Reference genome: "),
                                              actionButton("refreshPQC", "Quality Control", icon = icon("sync-alt")),
                                              hr(),
                                              DT::dataTableOutput("pQC.results")
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
                )