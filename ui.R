library(shiny)
library(DT)

#### for primer design ########
primersDesign_wd <- getwd() 

#### for primer QC ########
library(devtools)
library(rBLAST)
library(rtracklayer)
library(BSgenome)
library(Biostrings) #for primer QC and Reads Extraction 

############ for Flowcell QC ########
choices_info <- data.frame(path=list.files(paste(primersDesign_wd, "\\flowcell_package\\newStruct_trimmed", sep=""),full.names=TRUE, pattern =".csv"))
print(choices_info[["path"]])
choices_info[["nameS"]]<-sapply(strsplit(as.character(choices_info[["path"]]),"/"),function(x) x[length(x)])
print(choices_info[["nameS"]])

######### for Reads Extraction ###########
# the script needs Trim galore package 
# HTML needs Trim galore as well 
flowcell_folders <- data.frame(path=list.dirs(primersDesign_wd,full.names=TRUE,recursive = FALSE))
packages_names <- sapply(flowcell_folders$path,function(x) basename(as.character(x)))

####### For Reads Alignment #########
wrapper_file <- file.path(primersDesign_wd, "wrapper.r", fsep=.Platform$file.sep)
source(wrapper_file)
library(seqinr)

library(shiny)

####### For Tooltips #########
library(shinyBS)
library(tippy)

####### For adding tooltips #######
#to use this, the shinyBS package needs to be installed: packages.install("shinyBS")
list("bsTooltip")
list("bsPopover")

# Define UI for application that draws a histogram
shinyUI(navbarPage(title=div(img(src="EpiPrimerLogo.png"), height="10", width="10", align="left"), #"EpiPrimer",
                   
                   ########## Overview of the Workflow ########
                   
                   ############  Procedure of Primers Design  #############
                   
                   tabPanel("Primers Design",
                            sidebarLayout(
                              sidebarPanel(
                                tippy("Upload the region for which you want to calculate potential primers. Your file should contain either a specification by region or sequence. Download example files here:","Your input must contain the fields 'chr', 'start', 'end', 'assembly' and 'sequenceID'"),
                                hr(),
                                downloadButton("downloadSequenceFile", "example sequence file"),
                                hr(),
                                downloadButton("downloadRegionsFile", "example regions file"),
                                hr(),
                                # Upload data:
                                fileInput("file", "Upload file:"),
                                hr(),
                                textInput("name","Dataset name:", paste0("PrimerSet",Sys.Date(),",",format(Sys.time(), "%X"))),
                                bsTooltip("name", "Choose a name for your primer design folder", "top", "hover"),
                                hr(),
                                #checkboxInput for primer options: 
                                #checkboxInput("i_low.complexity.primer.removal", label = h4("remove primers of low complexity"), TRUE),
                                checkboxInput("i_remove.primers.with.n", label = h4("Remove primers that contain N bases"), TRUE),
                                checkboxInput("i_check4snps", label = h4("Check for SNPs"), TRUE),
                                conditionalPanel(
                                  "input.i_check4snps ==1",
                                  sliderInput("i_snps.amplicon", label = h4("Number of SNPs allowed in the amplicon"),
                                              min = 0, max = 80, value = 0),
                                  sliderInput("i_snps.primer1", label = h4("Number of SNPs allowed in the forward primer"),
                                              min = 0, max = 80, value = 0),
                                  sliderInput("i_snps.primer2", label = h4("Number of SNPs allowed in the reverse primer"),
                                              min = 0, max = 80, value = 0)
                                ),
                                checkboxInput("i_check4repeats", label = h4("Check for repeats"), FALSE),
                                conditionalPanel(
                                  "input.i_check4repeats == 1",
                                  checkboxInput("i_allow.repeats.in.primers", label = h4("Allow repeats in primers"), FALSE),
                                  checkboxInput("i_allow.repeats.in.amplicon", label = h4("Allow repeats in Amplicons"), FALSE)
                                ),
                                checkboxInput("i_annotate.genes", label = h4("Annotate genes"), FALSE),
                                #checkboxInput("i_annotate.cpg.islands)", label = h4("Annotate CpG islands"), FALSE),
                                checkboxInput("adapterF", label = h4("Add a specific sequence to 5' end of forward primer"), FALSE),
                                conditionalPanel(
                                  "input.adapterF == 1",
                                  textInput("adapterForward", "Forward adapter: ", "TCTTTCCCTACACGACGCTCTTCCGATCT")
                                ),
                                checkboxInput("adapterR", label = h4("Add a specific sequence to 5' end of reverse primer"), FALSE),
                                conditionalPanel(
                                  "input.adapterR==1",
                                  textInput("adapterReverse", "Reverse adapter: ", "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
                                ),
                                
                                #compute primers
                                actionButton("action", label="Compute Primers", icon("fas fa-calculator"), 
                                             style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:25px; font-size:200%; width:400px"),
                                bsTooltip("action", "Computation takes approx. 4 minutes, please wait until you receive a notifiication"),
                                hr(),
                                #download primers
                                downloadButton('downloadPrimer', 'Download Primers',
                                              style="padding:25px; font-size:200%; width:400px"),
                                bsTooltip("downloadPrimer", "Download your primers to a local folder"),
                                br(),
                                textOutput("state"),
                                br()
                              ),
                              # Main:
                              mainPanel(
                                DT::dataTableOutput("table"),
                                h2(strong("Basic primer settings: ")),
                                fluidRow(
                                  column(6,
                                         radioButtons("i_primer_type", label = h3("Primer Type"),
                                                      choices = list("Genomic"="genomic", "Bisulfite" = "bisulfite", "NOME" = "NOME", "CLEVER"="CLEVER",
                                                                      "Bisulfite (hairpin)" = "hp_bisulfite", "NOME (hairpin)" = "hp_NOME", "CLEVER (hairpin)"="hp_CLEVER", "CrispRCas9 Amplicon"="CrispRCas9PCR"),
                                                      selected = "genomic"), 
                                         
                                         bsTooltip("i_primer_type", "What kind of primer do you want to create?", "bottom", "hover")
                                  ),
                                  
                                  column(6,
                                         conditionalPanel(
                                           condition = "input.i_primer_type != 'genomic' && input.i_primer_type != 'CrispRCas9PCR'",
                                           radioButtons("i_strand", label = h3("Strand"),
                                                        choices = list("Top" = "top", "Bottom" = "bottom", "Both" = "both"),
                                                        selected = "top"),
                                           bsTooltip("i_strand", "Choose the strand for which you want to create your primers!", "bottom", "hover")
                                         )
                                  )
              
                              ),
                              
                              h2(strong("Advanced primer settings: ")),
                                fluidRow(
                                    column(6,
                                            sliderInput("i_max.bins.low.complexity", label = h4("Maximum length of monomeric base stretches"),
                                                        min = 0, max = 10, value = 7)
                                          ),
                                    column(6,
                                           sliderInput("i_primer.align.binsize", label = h4("Maximum length for Primer Self-Interaction"),
                                                        min = 0, max = 50, value = 12)
                                          )
                                ),
                                fluidRow(
                                    column(6,
                                           sliderInput("i_primerlength", label = h4("Primer Length"),
                                                        min = 10, max = 80, value = c(23, 34))
                                           ),
                                    column(6,
                                           sliderInput("i_primertemp", label = h4("Primer Melting Temperature"),
                                                        min = 40, max = 75, value = c(48, 60))
                                           )
                                ),
                                fluidRow(
                                
                                column(6,
                                       sliderInput("i_meltdiff", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                                                    min = 0, max = 10, value = 2)
                                       ),
                                column(6,
                                       sliderInput("i_lengthAmp", label = h4("Amplicon Length"),
                                                    min = 100, max = 800, value = c(150,500))
                                       )
                                ),
                                fluidRow(
                                column(6,
                                       sliderInput("i_minGC", h4("Minimum Number of GCs per amplicon"),
                                                    min = 0, max = 10, value = 0)
                                       ),
                                column(6,
                                       sliderInput("i_minCG", h4("Minimum Number of CGs per amplicon"),
                                                    min = 0, max = 10, value = 5)
                                       )
                                ),
                                fluidRow(
                                  column(6,
                                       conditionalPanel(
                                         condition = "input.i_primer_type == 'hp_bisulfite' || input.i_primer_type == 'hp_NOME' || input.i_primer_type == 'hp_CLEVER'", 
                                         sliderInput("i_hp.length", label = h4("length of one arm in the hairpin molecule"),
                                                     min = 0, max = 1000, value = c(50, 300))
                                                        )
                                       ),
                                  column(6,
                                       conditionalPanel(
                                         condition = "input.i_primer_type == 'genomic'",
                                         sliderInput("i_chop.size", label = h4("input sequence slicing"),
                                                     min = 0, max = 50, value = 30)
              
                                       ))
                                ),
                                fluidRow(
                                column(6,
                                     conditionalPanel(
                                       condition="input.i_primer_type != 'genomic' && input.i_primer_type != 'CrispRCas9PCR'",
                                       sliderInput("i_minC2T", h4("Minimum 'C' to 'T' conversions in forward primer"),
                                                   min = 0, max = 10, value = 3)
                                     )
                                ),
                                column(6,
                                     conditionalPanel(
                                       condition="input.i_primer_type != 'genomic' && input.i_primer_type != 'CrispRCas9PCR'",
                                       sliderInput("i_minG2A", h4("Minimum 'G' to 'A' conversions in reverse primer"),
                                                   min = 0, max = 10, value = 3)
                                     )
                                )
                                )
                              )
                              )
                          ),
                   
                   ############  UI For Results of Primers Design  #############
                   
                   tabPanel("Results of Primers Design",
                            sidebarLayout(
                              sidebarPanel(
                                helpText(h4("To analyze the results of your primer design use the buttons below:")),
                                helpText("Find out for which sequences primers could be created by pushing this button"),
                                actionButton("primerdesigns.by.sequence", label = "results overview"),
                                #tippy_this("primersdesigns.by.sequence", "Tooltip", "right"),
                                #bsTooltip("primerdesigns.by.sequence", "which primers could be created?", "bottom", "hover"),
                                hr(),
                                helpText("A list of the best rated primers will be shown after clicking this button:"),
                                actionButton("toplist", label = "show toplist"),
                                hr(),
                                helpText("Generate a list of selected primers by marking primers from the toplist above and pushing this button. You can either download these primers using the button on the right or further analyze them in the Primer QC panel. (Use the 'Import Primers' button there)"),
                                actionButton("selectlist", label = "generate selected list"), downloadButton('downloadSelectedPrimers', 'Selected Primers'),
                                hr(),
                                helpText("Show the whole list of created primers with this button: "),
                                actionButton("wholelist", label = "show The whole list"),
                                hr(),
                                helpText("See the sequences, for which no primers  could be designed:"),
                                actionButton("blacklist", label = "show blacklist"),
                                hr(),
                                helpText("See the sequences for which primers could be designed:"),
                                actionButton("whitelist", label = "show whitelist"),
                                hr(),
                                helpText("View the computational details of this primer design"),
                                actionButton("logfile", label = "show logfile"),
                                hr(),
                                helpText("Review the settings of this primer design: "),
                                actionButton("settings", label = "show settings"),
                                hr(),
                                helpText("Show further details of the primer design, such as duration and number of created primers:"),
                                actionButton("Summary", label = "show Summary"),
                                hr()
                              ),
                              
                              # Main:
                              mainPanel(
                                h4("Primer Designs For Each Sequence"),DT::dataTableOutput("viewprimerdesigns"),br(),
                                h4("Top List"),DT::dataTableOutput("viewtoplist"),
                                br(),
                                h4("The List of selected primers"),
                                DT::dataTableOutput("viewSelectlist"),
                                br(),
                                h4("The Whole List"),DT::dataTableOutput("viewwholelist"),
                                br(),
                                h4("Black List"),DT::dataTableOutput("viewblacklist"),
                                br(),
                                h4("White List"),DT::dataTableOutput("viewwhitelist"),
                                br(),
                                h4("Log File"),
                                DT::dataTableOutput("viewlogfile"),br(),
                                h4("Settings"),DT::dataTableOutput("viewsettings"),br(),
                                h4("Summary of the Primer Design"),DT::dataTableOutput("viewSummary")
                              )
                            )
                   ),
                   
                   ########### UI For Graphs of primers design procedure #############
                   
                   tabPanel("Graphs of Primers Design",
                            sidebarLayout(
                              sidebarPanel(
                                helpText(h4("Graphs of the primer design will be shown after clicking this button:")),
                                actionButton("graphics", label = "show figures")
                              ),
                              # Main:
                              mainPanel(
                                h4("Graphs"),
                                uiOutput("plot3")
                              )
                            )
                   ),
                   
                   tabPanel("Primer QC",
                            sidebarLayout(
                              sidebarPanel(
                                # upload files
                                fileInput("Fprimers", "Upload Forward Primers", multiple = TRUE ,accept = ".fasta"), hr(),
                                fileInput("Rprimers", "Upload Reverse Primers", multiple = TRUE, accept = ".fasta"), hr(),
                                bsTooltip("Rprimers", "tooltip", "top", "hover"), 
                                helpText("Import primers created by Primer Design step:"),
                                actionButton("loadprimers", "Import Primers"),
                                bsTooltip("loadprimers", "Import primers you created from the Primer Design step", "top", "hover"), hr(),
                                #selectInput("genome", "Choose the genome",choices=c(data.frame(spec=unlist(lapply(strsplit(installed.genomes(), "[.]"), "[", 4))))),
                                selectInput("genome", "Choose the genome",choices=c(installed.genomes())),
                                bsTooltip("genome", "Select the genome against which you want to blast your primers!", "top", "hover"),
                                sliderInput("gap", "Maximum Fragment Size:", min = 0, max = 50000, value = 2000),
                                helpText("Please set an E-value for your quality control. This value can be used to filter your results by only returning results, that are equal or better than the E-Value."),
                                sliderInput("Evalue", "E Value:", min = 0, max = 100, value = 10),
                                actionButton("computePQC", "Start Primers QC"),
                                bsTooltip("computePQC", "A virtual PCR of your primers is being computed. The results show potential PCR products resulting form your choice of primers", "top", "hover"),
                                h3("Parameters to filter Results"),
                                sliderInput("FMismatches", "Forward Primers Mismatches", min = 0, max = 5, value = 3),
                                sliderInput("RMismatches", "Reverse Primers Mismatches", min = 0, max = 5, value = 3),
                                sliderInput("FbitScore", "Forward Primers Bit Score", min = 0, max = 100, value = 25),
                                sliderInput("RbitScore", "Reverse Primers Bit Score", min = 0, max = 100, value = 25) 
                              ),
                              
                              # Show the results
                              mainPanel(
                                h4("Forward Primers"),
                                DT::dataTableOutput("forward.primers"),
                                h4("Reverse Primers"),
                                DT::dataTableOutput("reverse.primers"),
                                h4("QC Overview"),
                                DT::dataTableOutput("pQC.results"),
                                textOutput("PrimerqcState"),
                                DT::dataTableOutput("viewSelectedpQC"),
                                actionButton("extractRegions", "Export Results to Reads Analysis"),
                                bsTooltip("extractRegions", "nothing is happening here yet.", "top", "hover")
                              )
                            )
                   ),
                   
                   ############## displaying the overview of the flowcell #########
                   
                   # App title ----
                   tabPanel("Flowcell QC",
                            
                            # Sidebar layout with input and output definitions ----
                            sidebarLayout(
                              
                              # Sidebar panel for inputs ----
                              sidebarPanel(
                                
                                # Input: Select a dataset ----
                                selectInput("flowcellPackage", "Choose a flow cell dataset:",
                                            choices = choices_info[["nameS"]]),
                                sliderInput("freq", "Occurance",
                                            min = 0, max = 150000, value = 6,step = 1),
                                
                                # Button
                                downloadButton("downloadData", "Download Overview"),
                                actionButton("update", "Update View")
                              ),
                              
                              
                              mainPanel(
                                
                                h4("Flow Cell Summary"),
                                DT::dataTableOutput("view")
                                
                              )
                              
                            )
                   ),
                   
                   ############# UI For Reads Extraction ###############
                   
                   tabPanel("Reads Extraction",
                            
                            # Sidebar layout with input and output definitions ----
                            sidebarLayout(
                              
                              # Sidebar panel for inputs ----
                              sidebarPanel(
                                
                                # Input: Select a dataset ----
                                # Upload primer list:
                                fileInput("primers.file", "Upload Primers File"),
                                # actionButton("import.primers", "Import Primers"),hr(),
                                fileInput("linker.file", "Upload Linker File"),
                                selectInput("results.format", "Format of Output Files", choices = c("FASTA","FASTQ"),selected="FASTQ"),
                                textInput("results.folderName","Name of the Output file:",paste0("Extracted.Reads")),
                                
                                sliderInput("min.seq", "Minimum sequences to report",
                                            min = 0, max = 2000, value = 20,step = 10),
                                sliderInput("mismatches", "Mismatch Value for Files Filtration",
                                            min = 0, max = 5, value = 1,step = 1),
                                sliderInput("error.rate", "Error Rate for Primer Matching",
                                            min = 0, max = 1, value = 0.1,step = 0.1),
                                #sliderInput("Quality.Cutoff", "Quality Cutoff Applied in The Trimming Step",
                                #min = 0, max = 40, value = 20,step = 2),

                                checkboxInput("merged.files", label = "Merge Output Files", TRUE),
                                checkboxInput("separated.files", label = "Separate Output Files", FALSE)
                              ),
                              
                              # Main panel for displaying outputs ----
                              mainPanel(
                                
                                fluidRow(
                                  column(4,selectInput("selectPackage", "choose flowcell package", choices = packages_names)),
                                  
                                  column(6,selectInput("selection.mode", "How to Select Reads", choices = list("Use Only First Primer" = "first", "Use Only Second Primer" = "second",
                                                                                                               "Fragments Matching Either Primer" = "union", "Fragments Matching Both Primers" = "intersection"),selected="first"))
                                ),
                                uiOutput("checkbox"),#checkboxInput('selectall', 'Select All/Deselect All'),
                                fluidRow(
                                  column(4,actionLink("selectall","Select All/Deselect All"))),
                                br(),
                                
                                textOutput("state1"),
                                actionButton("excute", "Run Reads Extraction"),
                                downloadButton("downloadData2", "Download Extracted Reads"),br(),br()
                                # helpText(h4(HTML("Tool Parameters: <br/> <br/> <br/>
                                # > Primer matching error rate: <br/>
                                # > Minimum sequences to report: Controls the minimum number of reads reported in the results <br/>
                                # > Mismatch Value: controls the number of reads files that will be used for reads extraction <br/>
                                # > Merge Output Files: Decides whether the reads in output files will be merged <br/>
                                # > Separate Output Files: Decides whether the reads in output files will be separated <br/>
                                # ")))
                              )
                            )
                   ),
                   
                   ###############Reads Alignemnt #############

                   # App title ----
                   tabPanel("Reads Alignment",
                            # Sidebar layout with input and output definitions ----
                            sidebarLayout(
                              # Sidebar panel for inputs ----
                              sidebarPanel(
                                # Input: Select a dataset ----
                                # Upload primer list:
                                fileInput("reads.file", "Upload Reads File(s)",multiple=TRUE),
                                fileInput("reference.file", "Upload Reference File(s)",multiple=TRUE),
                                selectInput("alignment.mode", "Alignment Mode", choices = c("bisulfite","Non-bisulfite"),selected="bisulfite"),
                                selectInput("i_alignment.type", "Alignment Type", choices = c("global-local","global","local","local-global"),selected="global-local"),
                                textInput("i_reference.sequence.id","Reference Sequence ID:",paste0("Extracted Ref")),
                                textInput("i_sample.id","Reads Sample ID:",paste0("Extracted Reads"))
                                # sliderInput("i_minimum.alignment.score", "Minimum Alignment Score",
                                # min = 0, max = 6000, value = 20,step = 10),
                                # sliderInput("i_minimum.percentage.identity", "Minimum Percentage Identity",
                                # min = 0, max = 100, value = 5,step = 1)
                              ),
                              # Main panel for displaying outputs ----
                              mainPanel(
                                DT::dataTableOutput("AlignmentResults"),
                                textOutput("Alignmentstate"),br(),
                                actionButton("excuteAlignment", "Excute Alignment"),
                                downloadButton("downloadAlignmentResults", "Download Extracted Reads"),br(),br()
                                # helpText(h4(HTML("Tool Parameters: <br/> <br/> <br/>
                                # > Primer matching error rate: <br/>
                                # > Minimum sequences to report: Controls the minimum number of reads reported in the results <br/>
                                # > Mismatch Value: controls the number of reads files that will be used for reads extraction <br/>
                                # > Merge Output Files: Decides whether the reads in output files will be merged <br/>
                                # > Separate Output Files: Decides whether the reads in output files will be separated <br/>
                                # ")))
                              )
                            )
                   ),
                   ################ Reads Analysis ######################
                   tabPanel("Reads Analysis",
                            
                            sidebarLayout(
                              sidebarPanel(
                                helpText(h4("Extracted Reads will be aligned against the extracted region of the reference Genome 
                                            that was exported from The Primers QC process. Results will be displayed with helpful plots and heats Maps."))
                                ),
                              # Main:
                              mainPanel(
                                helpText(h4("To be implemented Soon!"))
                              )
                              )
                   )
))
