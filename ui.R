## app.R ##
library(shiny)
library(shinydashboard)

## tooltips ##
library(shinyBS)
library(tippy)

dbHeader <- dashboardHeader(title = "EpiPrimer")

## to add CrispRCas9 as a primer type for primer design, simply extend the radio buttons i_primer_type with 'CrispRCas9PCR'! 
## Default settings for this primer type are already contained in the ui.R and server.R file!
## Computations of this type of primer are also possible (recomment the function contained in primer.design.r) BUT are NOT refactored.

## UI using shiny dashboard ##
ui <- fluidPage(dashboardPage(skin = "yellow",
                              dashboardHeader(title = "EpiPrimer"),#dbHeader,
                              dashboardSidebar(
                                sidebarMenu(id = "tabs",
                                  menuItem("Primer Design Start", tabName = "PrimerDesign", icon = icon("dna")),
                                  menuItem("Results of Primer Design", tabName = "PDresults", icon = icon("list-ol")),
                                  menuItem("ePCR Start", tabName = "PrimerQC", icon = icon("dna")), #check-circle
                                  menuItem("Results of ePCR", tabName = "PrimerQCResults", icon = icon("list-ol")),
                                  menuItem("Imprint", tabName = "Imprint", icon = icon("university"))
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
                                                helpText("Human genome assemblies: hg18, hg19, hg38"),
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
                                          helpText("You can view the results of your primer design job using the 'Results' tab."),
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
                                                                          "Bisulfite (hairpin)" = "hp_bisulfite", "NOME (hairpin)" = "hp_NOME", "CLEVER (hairpin)"="hp_CLEVER"), #"CrispRCas9 Amplicon"="CrispRCas9PCR"),
                                                           selected = "genomic"),
                                              bsTooltip("i_primer_type", "What kind of primer do you want to create?", "left", "hover"),
                                              hr(),
                                              conditionalPanel(
                                                condition = "input.i_primer_type != 'genomic' && input.i_primer_type != 'CrispRCas9PCR'",
                                                radioButtons("i_strand", label = h3("Strand"),
                                                             choices = list("Top" = "top", "Bottom" = "bottom", "Both" = "both"),
                                                             selected = "top"),
                                                bsTooltip("i_strand", "Choose the strand for which you want to create your primers!", "left", "hover")
                                              )
                                          ),
                                          hr(),
                                          box(
                                            fluidRow(column(width=6,
                                                            sliderInput("i_snps.amplicon", label = h4("Number of SNPs allowed in the Amplicon"),
                                                                        min = 0, max = 80, value = 10),
                                                            sliderInput("i_snps.primer1", label = h4("Number of SNPs allowed in the Forward Primer"),
                                                                        min = 0, max = 80, value = 0),
                                                            sliderInput("i_snps.primer2", label = h4("Number of SNPs allowed in the Reverse Primer"),
                                                                        min = 0, max = 80, value = 0),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'genomic'",
                                                              sliderInput("i_max.bins_genomic", label = h4("Maximum length of monomeric base stretches"),
                                                                          min = 0, max = 10, value = 5)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'bisulfite'",
                                                              sliderInput("i_max.bins_bis", label = h4("Maximum length of monomeric base stretches"),
                                                                          min = 0, max = 10, value = 7)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'NOME'",
                                                              sliderInput("i_max.bins_NOME", label = h4("Maximum length of monomeric base stretches"),
                                                                          min = 0, max = 10, value = 7)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'CLEVER'",
                                                              sliderInput("i_max.bins_CLEVER", label = h4("Maximum length of monomeric base stretches"),
                                                                          min = 0, max = 10, value = 5)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_bisulfite'",
                                                              sliderInput("i_max.bins_hp_bis", label = h4("Maximum length of monomeric base stretches"),
                                                                          min = 0, max = 10, value = 7)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_NOME'",
                                                              sliderInput("i_max.bins_hp_NOME", label = h4("Maximum length of monomeric base stretches"),
                                                                          min = 0, max = 10, value = 7)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_CLEVER'",
                                                              sliderInput("i_max.bins_hp_clever", label = h4("Maximum length of monomeric base stretches"),
                                                                          min = 0, max = 10, value = 5)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'CrispRCas9PCR'",
                                                              sliderInput("i_max.bins_crispr", label = h4("Maximum length of monomeric base stretches"),
                                                                          min = 0, max = 10, value = 5)
                                                            ),
                                                            sliderInput("i_primer.align.binsize", label = h4("Maximum length for Primer Self-Interaction"),
                                                                        min = 0, max = 50, value = 12),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'genomic'",
                                                              sliderInput("i_primerlength_genomic", label = h4("Primer Length"),
                                                                          min = 10, max = 80, value = c(18, 25))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'bisulfite'",
                                                              sliderInput("i_primerlength_bis", label = h4("Primer Length"),
                                                                          min = 10, max = 80, value = c(20, 32))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'NOME'",
                                                              sliderInput("i_primerlength_NOME", label = h4("Primer Length"),
                                                                          min = 10, max = 80, value = c(20, 32))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'CLEVER'",
                                                              sliderInput("i_primerlength_clever", label = h4("Primer Length"),
                                                                          min = 10, max = 80, value = c(18, 25))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_bisulfite'",
                                                              sliderInput("i_primerlength_hp_bis", label = h4("Primer Length"),
                                                                          min = 10, max = 80, value = c(20, 32))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_NOME'",
                                                              sliderInput("i_primerlength_hp_NOME", label = h4("Primer Length"),
                                                                          min = 10, max = 80, value = c(20, 32))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_CLEVER'",
                                                              sliderInput("i_primerlength_hp_clever", label = h4("Primer Length"),
                                                                          min = 10, max = 80, value = c(18, 25))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'genomic'",
                                                              sliderInput("i_primertemp_genomic", label = h4("Primer Melting Temperature"),
                                                                          min = 40, max = 75, value = c(50, 60))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'bisulfite'",
                                                              sliderInput("i_primertemp_bis", label = h4("Primer Melting Temperature"),
                                                                          min = 40, max = 75, value = c(48, 60))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'NOME'",
                                                              sliderInput("i_primertemp_NOME", label = h4("Primer Melting Temperature"),
                                                                          min = 40, max = 75, value = c(48, 60))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'CLEVER'",
                                                              sliderInput("i_primertemp_clever", label = h4("Primer Melting Temperature"),
                                                                          min = 40, max = 75, value = c(50, 60))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_bisulfite'",
                                                              sliderInput("i_primertemp_hp_bis", label = h4("Primer Melting Temperature"),
                                                                          min = 40, max = 75, value = c(48, 60))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_NOME'",
                                                              sliderInput("i_primertemp_hp_NOME", label = h4("Primer Melting Temperature"),
                                                                          min = 40, max = 75, value = c(48, 60))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_CLEVER'",
                                                              sliderInput("i_primertemp_hp_clever", label = h4("Primer Melting Temperature"),
                                                                          min = 40, max = 75, value = c(50, 60))
                                                            )
                                            ),
                                            column(width = 6,
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'genomic'",
                                                     sliderInput("i_meltdiff_genomic", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                                                                 min = 0, max = 10, value = 3)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'bisulfite'",
                                                     sliderInput("i_meltdiff_bis", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                                                                 min = 0, max = 10, value = 5)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'NOME'",
                                                     sliderInput("i_meltdiff_NOME", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                                                                 min = 0, max = 10, value = 6)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'CLEVER'",
                                                     sliderInput("i_meltdiff_clever", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                                                                 min = 0, max = 10, value = 3)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'hp_bisulfite'",
                                                     sliderInput("i_meltdiff_hp_bis", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                                                                 min = 0, max = 10, value = 5)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'hp_NOME'",
                                                     sliderInput("i_meltdiff_hp_NOME", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                                                                 min = 0, max = 10, value = 6)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'hp_CLEVER'",
                                                     sliderInput("i_meltdiff_hp_clever", label = h4("Maximum Difference in Primer Melting Temperature \n (degrees Celsius)"),
                                                                 min = 0, max = 10, value = 3)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'genomic'",
                                                     sliderInput("i_lengthAmp_genomic", label = h4("Amplicon Length"),
                                                                 min = 100, max = 800, value = c(200,500))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'bisulfite'",
                                                     sliderInput("i_lengthAmp_bis", label = h4("Amplicon Length"),
                                                                 min = 100, max = 800, value = c(200,400))
                                                   ),  
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'NOME'",
                                                     sliderInput("i_lengthAmp_NOME", label = h4("Amplicon Length"),
                                                                 min = 100, max = 800, value = c(200,400))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'CLEVER'",
                                                     sliderInput("i_lengthAmp_clever", label = h4("Amplicon Length"),
                                                                 min = 100, max = 800, value = c(200,500))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_bisulfite'",
                                                     sliderInput("i_lengthAmp_hp_bis", label = h4("Amplicon Length"),
                                                                 min = 100, max = 800, value = c(200,400))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_NOME'",
                                                     sliderInput("i_lengthAmp_hp_NOME", label = h4("Amplicon Length"),
                                                                 min = 100, max = 800, value = c(200,400))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_CLEVER'",
                                                     sliderInput("i_lengthAmp_hp_clever", label = h4("Amplicon Length"),
                                                                 min = 100, max = 800, value = c(200,500))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'CrispRCas9PCR'",
                                                     sliderInput("i_lengthAmp_crispr", label = h4("Amplicon Length"),
                                                                 min = 100, max = 800, value = c(200,500))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'genomic'",
                                                     sliderInput("i_minGC_genomic", h4("Minimum Number of GCs per amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'bisulfite'",
                                                     sliderInput("i_minGC_bis", h4("Minimum Number of GCs per amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'NOME'",
                                                     sliderInput("i_minGC_NOME", h4("Minimum Number of GCs per amplicon"),
                                                                 min = 0, max = 10, value = 5)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'CLEVER'",
                                                     sliderInput("i_minGC_clever", h4("Minimum Number of GCs per amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_bisulfite'",
                                                     sliderInput("i_minGC_hp_bis", h4("Minimum Number of GCs per amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_NOME'",
                                                     sliderInput("i_minGC_hp_NOME", h4("Minimum Number of GCs per amplicon"),
                                                                 min = 0, max = 10, value = 1)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_CLEVER'",
                                                     sliderInput("i_minGC_hp_clever", h4("Minimum Number of GCs per amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'CrispRCas9PCR'",
                                                     sliderInput("i_minGC_crispr", h4("Minimum Number of GCs per amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'genomic'",
                                                     sliderInput("i_minCG_genomic", h4("Minimum Number of CGs per amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ), 
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'bisulfite'",
                                                     sliderInput("i_minCG_bis", h4("Minimum Number of CGs per amplicon"),
                                                                 min = 0, max = 10, value = 5)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'NOME'",
                                                     sliderInput("i_minCG_NOME", h4("Minimum Number of CGs per amplicon"),
                                                                 min = 0, max = 10, value = 5)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'CLEVER'",
                                                     sliderInput("i_minCG_clever", h4("Minimum Number of CGs per amplicon"),
                                                                 min = 0, max = 10, value = 5)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_bisulfite'",
                                                     sliderInput("i_minCG_hp_bis", h4("Minimum Number of CGs per amplicon"),
                                                                 min = 0, max = 10, value = 1)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_NOME'",
                                                     sliderInput("i_minCG_hp_NOME", h4("Minimum Number of CGs per amplicon"),
                                                                 min = 0, max = 10, value = 1)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_CLEVER'",
                                                     sliderInput("i_minCG_hp_clever", h4("Minimum Number of CGs per amplicon"),
                                                                 min = 0, max = 10, value = 1)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'CrispRCas9PCR'",
                                                     sliderInput("i_minCG_crispr", h4("Minimum Number of CGs per amplicon"),
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
                                            )
                                            ),
                                            hr(),
                                            h2("Add Adapters to my Primers: "),
                                            # adapter box
                                            fluidRow(column(12,
                                                            textInput("adapterForward", h4("Add a specific sequence to the 5' end of the forward primer: "), ""), #CTTTCCCTACACGACGCTCTTCCGATCT
                                                            textInput("adapterReverse", h4("Add a specific sequence to the 5' end of the reverse primer: "), "") #GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
                                            )),
                                            hr(),
                                            h2("Further Settings: "),
                                            fluidRow(column(12,
                                                            checkboxInput("i_remove.primers.with.n", label = h4("Remove primers that contain N bases"), TRUE),
                                                            checkboxInput("i_allow.repeats.in.primers", label = h4("Allow repeats in primers"), FALSE),
                                                            checkboxInput("i_allow.repeats.in.amplicon", label = h4("Allow repeats in Amplicons"), FALSE)
                                                            )
                                                     )
                                            
                                            
                                            ,
                                            title = h2("Advanced Primer Settings"),
                                            width = 12,
                                            status = "primary",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            collapsed = TRUE
                                          ),
                                          fluidRow(
                                            actionButton("action", label="Compute Primers", icon("fas fa-calculator"), align="center",
                                                         style="color: #fff; background-color: #3c8dbc; border-color: #337ab7; padding:25px; font-size:200%; width:1000px; margin-left:200px; margin-right:0px"),
                                            bsTooltip("action", "The computation of your primers may take a few minutes, please wait until you receive a notifiication that your primers are finished.", "left")
                                          )
                                          
                                                )
                                          
                                          ),
                                  
                                  tabItem(tabName = "PDresults",
                                          fluidRow(
                                            tabBox(title = "",
                                                   id = "PDresultsTabbox",
                                                   width = 12,
                                                   tabPanel("Overview",
                                                            #actionButton("primerdesigns.by.sequence", label = "Overview", icon = icon("sync-alt")),
                                                            #hr(),
                                                            DT::dataTableOutput("viewprimerdesigns")
                                                   ),
                                                   tabPanel("Toplist",
                                                            #actionButton("toplist", label = "Toplist", icon = icon("sync-alt")),
                                                            #hr(),
                                                            DT::dataTableOutput("viewtoplist")
                                                   ),
                                                   tabPanel("Wholelist",
                                                            #actionButton("wholelist", label = "Wholelist", icon = icon("sync-alt")),
                                                            #hr(),
                                                            DT::dataTableOutput("viewwholelist")
                                                   ),
                                                   tabPanel("Whitelist",
                                                            #actionButton("whitelist", label = "Whitelist", icon = icon("sync-alt")),
                                                            #hr(),
                                                            DT::dataTableOutput("viewwhitelist")
                                                   ),
                                                   tabPanel("Blacklist",
                                                            #actionButton("blacklist", label = "Blacklist", icon = icon("sync-alt")),
                                                            #hr(),
                                                            DT::dataTableOutput("viewblacklist")
                                                   ),
                                                   tabPanel("Logfile",
                                                            #actionButton("logfile", label = "Logfile", icon = icon("sync-alt")),
                                                            #hr(),
                                                            DT::dataTableOutput("viewlogfile")
                                                   ),
                                                   tabPanel("Settings",
                                                            #actionButton("settings", label = "Settings", icon = icon("sync-alt")),
                                                            #hr(),
                                                            DT::dataTableOutput("viewsettings")
                                                   ),
                                                   tabPanel("Summary",
                                                            #actionButton("Summary", label = "Summary", icon = icon("sync-alt")),
                                                            #hr(),
                                                            DT::dataTableOutput("viewSummary")
                                                   )
                                            )
                                          ),
                                          fluidRow(
                                            box(title = "Selected List",
                                                width = 12,
                                                helpText("You can download primer pairs or hand them over for ePCR by adding them to your list of selected primers below. To add primers to your Selected List, mark them on the Wholelist!"),
                                                hr(),
                                                #actionButton("selectlist", label = "Generate Selected List"),
                                                #hr(),
                                                DT::dataTableOutput("viewSelectlist"),
                                                hr(),
                                                downloadButton('downloadSelectedPrimers', 'Selected Primers')
                                            )
                                          ),
                                          box(#actionButton("graphics", label = "Graphs", icon = icon("sync-alt")),
                                              uiOutput("plot3"),
                                              title = h2("Graphs"),
                                              status = "primary",
                                              solidHeader = TRUE,
                                              collapsible = TRUE,
                                              width = 12
                                          )
                                          
                                  ),
                                  tabItem(tabName = "PrimerQC",
                                          box(title = h2("Upload files for ePCR"),
                                              status = "primary", 
                                              solidHeader = TRUE,
                                              width = 6,
                                              helpText("To check the quality of your primer pairs, we perform BLAST search against the Reference Genome and provide potential fragments during PCR.", align = "left"),
                                              hr(),
                                              h4("Import primer pairs from Primer Design that were added to your Selected List:", align = "center"),
                                              actionButton("loadprimers", "Import Primers", icon = icon("file-import"), align="center",
                                                           style="margin-left:275px; margin-right:0px"),
                                              bsTooltip("loadprimers", "Import primers you added to the Selected List during Primer Design", "bottom", "hover"),
                                              hr(),
                                              h4("Upload your own Primers for ePCR below:", align = "left"),
                                              helpText("(.fasta file format is required, we expect your primers to be in 5'-3' orientation)", align = "left"),
                                              fileInput("Fprimers", "Upload Forward Primers", multiple = TRUE, accept = ".fasta"), 
                                              DT::dataTableOutput("forward.primers"),
                                              fileInput("Rprimers", "Upload Reverse Primers", multiple = TRUE, accept = ".fasta"), 
                                              DT::dataTableOutput("reverse.primers")
                                          ),
                                          box(title = h2("Settings for ePCR"),
                                              status = "primary", 
                                              solidHeader = TRUE,
                                              width = 6,
                                              textInput("blast_id","Name your Primerblast:", paste0("MyBlast")),
                                              selectInput("genome", "Genome for Quality Control",choices = gsub("UCSC.", "", gsub("BSgenome.", "", c("BSgenome.Hsapiens.UCSC.hg18", "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm9", "BSgenome.Mmusculus.UCSC.mm10")))), #gsub("UCSC.", "", gsub("BSgenome.", "", installed.genomes()))), #(names(a) = strsplit(installed.genomes(), ".", fixed=TRUE)[1][4])),
                                              bsTooltip("genome", "Select the genome against which you want to blast your primers!", "top", "hover"),
                                              checkboxInput("is_bisulfite", h5("These are bisulfite primers!"), FALSE),
                                              sliderInput("gap", "Maximum size of reported product", min = 500, max = 10000, value = 2000, step = 500),
                                              sliderInput("primer_mismatches", "Number of Mismatches allowed in Primer Blast", min=0, max=25, value=7),
                                              hr(),
                                              textOutput("primer_qc_start"),
                                              tags$head(tags$style("#primer_qc_start{color: red;
                                                                     font-size: 30px;
                                                                   font-style: italic;
                                                                   }"
                                                        )
                                              ),
                                              hr()
                                          ),
                                          actionButton("compute_ePCR", label = "Start ePCR", icon("fas fa-flask"), align="center",
                                                       style="color: #fff; background-color: #3c8dbc; border-color: #337ab7; padding:25px; font-size:200%; width:1000px; margin-left:150px; margin-right:0px")
                                  
                                  ),
                                  tabItem(tabName = "PrimerQCResults",
                                          fluidRow(
                                            tabBox(title = "",
                                                   id = "PQCTabbox",
                                                   width = 12,
                                                   tabPanel("Results",
                                                            actionButton("refreshPQC", label = "Results", icon = icon("sync-alt")),
                                                            selectInput("test_select", label = "Filter results by primerpair: ", choices ="", multiple = FALSE),
                                                            #hr(),
                                                            tags$div(id="txt_for_selector", helpText("")),
                                                            uiOutput("pQC.results")
                                                   ),
                                                   tabPanel("Overview",
                                                            #actionButton("overviewePCR", label = "Overview", icon = icon("sync-alt")),
                                                            #hr(),
                                                            DT::dataTableOutput("ePCR.overview")
                                                            ),
                                                   tabPanel("Settings",
                                                            #actionButton("settingsePCR", label = "Settings", icon = icon("sync-alt")),
                                                            #hr(),
                                                            DT::dataTableOutput("ePCR.settings")
                                                            ),
                                                   tabPanel("Logfile",
                                                            #actionButton("logfileePCR", label = "Logfile", icon = icon("sync-alt")),
                                                            #hr(),
                                                            DT::dataTableOutput("ePCR.logfile")),
                                                   tabPanel("Summary",
                                                            #actionButton("summaryePCR", label = "Summary", icon = icon("sync-alt")),
                                                            #hr(),
                                                            DT::dataTableOutput("pQC.summary")
                                                   )#,
                                                   #hr()
                                                   #downloadButton('downloadePCRresults', 'Results')
                                            )
                                            ),
                                            
                                            fluidRow(box(#actionButton("graphics_ePCR", label = "Graphs", icon = icon("sync-alt")),
                                              uiOutput("ePCR_Graphs"),
                                              title = h2("Graphs for ePCR"),
                                              status = "primary",
                                              solidHeader = TRUE,
                                              collapsible = TRUE,
                                              collapsed = FALSE,
                                              width = 12
                                            )
                                          )
                                          
                                  ),
                                  tabItem(tabName = "Imprint",
                                          box(title = h2("Imprint"),
                                              status = "primary", 
                                              solidHeader = TRUE,
                                              width = 12,
                                              helpText("This website was created by the Epigenetics department of Saarland University."),
                                              tags$a(href="http://epigenetik.uni-saarland.de/en/home/", "Visit our Homepage here."),
                                              hr(),
                                              helpText("In case you encounter any problems during use of this website, please contact us!")
                                          )
                                  )
                                  
                                )
                          )
                      )
                )