## app ##
library(shiny)
library(shinydashboard)
library(shinycssloaders)

## tooltips ##
library(shinyBS)
library(tippy)

## to add CrispRCas9 as a primer type for primer design, simply extend the radio buttons i_primer_type with 'CrispRCas9PCR'! 
## Default settings for this primer type are already contained in the ui.R and server.R file!
## Computations of this type of primer are also possible (recomment the function contained in primer.design.r) BUT are NOT refactored.

## UI using shiny dashboard ##
ui <- fluidPage(dashboardPage(skin = "yellow",
                              dashboardHeader(title = "EpiPrimer"),
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
                                            box(title = h2("Upload File for Primer Design"),
                                                status = "primary", 
                                                solidHeader = TRUE,
                                                # Upload data:
                                                strong(helpText("Provided Genomes:")),
                                                helpText("Human Genome Assemblies: hg18, hg19, hg38"),
                                                helpText("Mouse Genome Assemblies: mm9, mm10"),
                                                hr(),
                                                fileInput("file", "Upload Regions or Sequence File:", accept=c("txt", "text/plain")),
                                                hr(),
                                                DT::dataTableOutput("table"),
                                                hr(),
                                                textInput("name","Dataset Name:", paste0("MyPrimerSet")),
                                                bsTooltip("name", "Choose a name for your primer design job", "top", "hover"),
                                                hr(),
                                                textOutput("state"),
                                                tags$head(tags$style("#state{color: red;
                                                                     font-size: 30px;
                                                                     font-style: italic;
                                                                     }"
                                          )
                                                ),
                                          hr(),
                                          helpText( h4("Download Example Files Here: ")),
                                          downloadButton("downloadSequenceFile", "Sequence File",
                                                         style="margin-left:100px; margin-right:0px"),
                                          downloadButton("downloadRegionsFile", "Regions File",
                                                         style="margin-left:100px; margin-right:0px")
                                                ),
                                          box(title = h2("Basic Primer Settings"), 
                                              status = "primary",
                                              solidHeader = TRUE,
                                              radioButtons("i_primer_type", label = h3("Primer Type"),
                                                           choices = list("Genomic"="genomic", "Bisulfite" = "bisulfite", "NOME" = "NOME", "CLEVER"="CLEVER",
                                                                          "Bisulfite (Hairpin)" = "hp_bisulfite", "NOME (Hairpin)" = "hp_NOME", "CLEVER (Hairpin)"="hp_CLEVER"), #"CrispRCas9 Amplicon"="CrispRCas9PCR"),
                                                           selected = "genomic"),
                                              bsTooltip("i_primer_type", "What Kind of Primer Do You Want to Create?", "left", "hover"),
                                              hr(),
                                              radioButtons("i_strand", label = h3("DNA Strand"),
                                                           choices = list("Top" = "top", "Bottom" = "bottom", "Both" = "both"),
                                                           selected = "top"),
                                              bsTooltip("i_strand", "Choose the Strand for Which You Want to Create Your Primers!", "left", "hover"),
                                              hr(),
                                              helpText("In Case Your Primer Design Job Was Unsatisfactory, You Are Recommended to Check the Advanced Primer Settings Below and to Compute Again.")
                                          ),
                                          hr(),
                                          box(
                                            fluidRow(column(width=6,
                                                            sliderInput("i_snps.amplicon", label = h4("Number of SNPs Allowed in the Amplicon"),
                                                                        min = 0, max = 80, value = 10),
                                                            sliderInput("i_snps.primer1", label = h4("Number of SNPs Allowed in the Forward Primer"),
                                                                        min = 0, max = 80, value = 0),
                                                            sliderInput("i_snps.primer2", label = h4("Number of SNPs Allowed in the Reverse Primer"),
                                                                        min = 0, max = 80, value = 0),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'genomic'",
                                                              sliderInput("i_max.bins_genomic", label = h4("Maximum Length of Monomeric Base Stretches"),
                                                                          min = 0, max = 10, value = 5)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'bisulfite'",
                                                              sliderInput("i_max.bins_bis", label = h4("Maximum Length of Monomeric Base Stretches"),
                                                                          min = 0, max = 10, value = 7)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'NOME'",
                                                              sliderInput("i_max.bins_NOME", label = h4("Maximum Length of Monomeric Base Stretches"),
                                                                          min = 0, max = 10, value = 7)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'CLEVER'",
                                                              sliderInput("i_max.bins_CLEVER", label = h4("Maximum Length of Monomeric Base Stretches"),
                                                                          min = 0, max = 10, value = 5)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_bisulfite'",
                                                              sliderInput("i_max.bins_hp_bis", label = h4("Maximum Length of Monomeric Base Stretches"),
                                                                          min = 0, max = 10, value = 7)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_NOME'",
                                                              sliderInput("i_max.bins_hp_NOME", label = h4("Maximum Length of Monomeric Base Stretches"),
                                                                          min = 0, max = 10, value = 7)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_CLEVER'",
                                                              sliderInput("i_max.bins_hp_clever", label = h4("Maximum Length of Monomeric Base Stretches"),
                                                                          min = 0, max = 10, value = 5)
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'CrispRCas9PCR'",
                                                              sliderInput("i_max.bins_crispr", label = h4("Maximum Length of Monomeric Base Stretches"),
                                                                          min = 0, max = 10, value = 5)
                                                            ),
                                                            sliderInput("i_primer.align.binsize", label = h4("Maximum Length for Primer Self-Interaction"),
                                                                        min = 0, max = 50, value = 12),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'genomic'",
                                                              sliderInput("i_primerlength_genomic", label = h4("Primer Length (bp)"),
                                                                          min = 10, max = 80, value = c(18, 25))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'bisulfite'",
                                                              sliderInput("i_primerlength_bis", label = h4("Primer Length (bp)"),
                                                                          min = 10, max = 80, value = c(20, 32))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'NOME'",
                                                              sliderInput("i_primerlength_NOME", label = h4("Primer Length (bp)"),
                                                                          min = 10, max = 80, value = c(20, 32))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'CLEVER'",
                                                              sliderInput("i_primerlength_clever", label = h4("Primer Length (bp)"),
                                                                          min = 10, max = 80, value = c(18, 25))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_bisulfite'",
                                                              sliderInput("i_primerlength_hp_bis", label = h4("Primer Length (bp)"),
                                                                          min = 10, max = 80, value = c(20, 32))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_NOME'",
                                                              sliderInput("i_primerlength_hp_NOME", label = h4("Primer Length (bp)"),
                                                                          min = 10, max = 80, value = c(20, 32))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_CLEVER'",
                                                              sliderInput("i_primerlength_hp_clever", label = h4("Primer Length (bp)"),
                                                                          min = 10, max = 80, value = c(18, 25))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'genomic'",
                                                              sliderInput("i_primertemp_genomic", label = h4("Primer Melting Temperature (\u00B0C)"),
                                                                          min = 40, max = 75, value = c(50, 60))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'bisulfite'",
                                                              sliderInput("i_primertemp_bis", label = h4("Primer Melting Temperature (\u00B0C)"),
                                                                          min = 40, max = 75, value = c(48, 60))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'NOME'",
                                                              sliderInput("i_primertemp_NOME", label = h4("Primer Melting Temperature (\u00B0C)"),
                                                                          min = 40, max = 75, value = c(48, 60))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'CLEVER'",
                                                              sliderInput("i_primertemp_clever", label = h4("Primer Melting Temperature (\u00B0C)"),
                                                                          min = 40, max = 75, value = c(50, 60))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_bisulfite'",
                                                              sliderInput("i_primertemp_hp_bis", label = h4("Primer Melting Temperature (\u00B0C)"),
                                                                          min = 40, max = 75, value = c(48, 60))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_NOME'",
                                                              sliderInput("i_primertemp_hp_NOME", label = h4("Primer Melting Temperature (\u00B0C)"),
                                                                          min = 40, max = 75, value = c(48, 60))
                                                            ),
                                                            conditionalPanel(
                                                              condition = "input.i_primer_type == 'hp_CLEVER'",
                                                              sliderInput("i_primertemp_hp_clever", label = h4("Primer Melting Temperature (\u00B0C)"),
                                                                          min = 40, max = 75, value = c(50, 60))
                                                            )
                                            ),
                                            column(width = 6,
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'genomic'",
                                                     sliderInput("i_meltdiff_genomic", label = h4("Maximum Difference in Primer Melting Temperature (\u00B0C)"),
                                                                 min = 0, max = 10, value = 3)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'bisulfite'",
                                                     sliderInput("i_meltdiff_bis", label = h4("Maximum Difference in Primer Melting Temperature (\u00B0C)"),
                                                                 min = 0, max = 10, value = 5)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'NOME'",
                                                     sliderInput("i_meltdiff_NOME", label = h4("Maximum Difference in Primer Melting Temperature (\u00B0C)"),
                                                                 min = 0, max = 10, value = 6)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'CLEVER'",
                                                     sliderInput("i_meltdiff_clever", label = h4("Maximum Difference in Primer Melting Temperature (\u00B0C)"),
                                                                 min = 0, max = 10, value = 3)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'hp_bisulfite'",
                                                     sliderInput("i_meltdiff_hp_bis", label = h4("Maximum Difference in Primer Melting Temperature (\u00B0C)"),
                                                                 min = 0, max = 10, value = 5)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'hp_NOME'",
                                                     sliderInput("i_meltdiff_hp_NOME", label = h4("Maximum Difference in Primer Melting Temperature (\u00B0C)"),
                                                                 min = 0, max = 10, value = 6)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type == 'hp_CLEVER'",
                                                     sliderInput("i_meltdiff_hp_clever", label = h4("Maximum Difference in Primer Melting Temperature (\u00B0C)"),
                                                                 min = 0, max = 10, value = 3)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'genomic'",
                                                     sliderInput("i_lengthAmp_genomic", label = h4("Amplicon Length (bp)"),
                                                                 min = 100, max = 800, value = c(200,500))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'bisulfite'",
                                                     sliderInput("i_lengthAmp_bis", label = h4("Amplicon Length (bp)"),
                                                                 min = 100, max = 800, value = c(200,400))
                                                   ),  
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'NOME'",
                                                     sliderInput("i_lengthAmp_NOME", label = h4("Amplicon Length (bp)"),
                                                                 min = 100, max = 800, value = c(200,400))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'CLEVER'",
                                                     sliderInput("i_lengthAmp_clever", label = h4("Amplicon Length (bp)"),
                                                                 min = 100, max = 800, value = c(200,500))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_bisulfite'",
                                                     sliderInput("i_lengthAmp_hp_bis", label = h4("Amplicon Length (bp)"),
                                                                 min = 100, max = 800, value = c(200,400))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_NOME'",
                                                     sliderInput("i_lengthAmp_hp_NOME", label = h4("Amplicon Length (bp)"),
                                                                 min = 100, max = 800, value = c(200,400))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_CLEVER'",
                                                     sliderInput("i_lengthAmp_hp_clever", label = h4("Amplicon Length (bp)"),
                                                                 min = 100, max = 800, value = c(200,500))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'CrispRCas9PCR'",
                                                     sliderInput("i_lengthAmp_crispr", label = h4("Amplicon Length (bp)"),
                                                                 min = 100, max = 800, value = c(200,500))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'genomic'",
                                                     sliderInput("i_minGC_genomic", h4("Minimum Number of GCs per Amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'bisulfite'",
                                                     sliderInput("i_minGC_bis", h4("Minimum Number of GCs per Amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'NOME'",
                                                     sliderInput("i_minGC_NOME", h4("Minimum Number of GCs per Amplicon"),
                                                                 min = 0, max = 10, value = 5)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'CLEVER'",
                                                     sliderInput("i_minGC_clever", h4("Minimum Number of GCs per Amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_bisulfite'",
                                                     sliderInput("i_minGC_hp_bis", h4("Minimum Number of GCs per Amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_NOME'",
                                                     sliderInput("i_minGC_hp_NOME", h4("Minimum Number of GCs per Amplicon"),
                                                                 min = 0, max = 10, value = 1)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_CLEVER'",
                                                     sliderInput("i_minGC_hp_clever", h4("Minimum Number of GCs per Amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'CrispRCas9PCR'",
                                                     sliderInput("i_minGC_crispr", h4("Minimum Number of GCs per Amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'genomic'",
                                                     sliderInput("i_minCG_genomic", h4("Minimum Number of CGs per Amplicon"),
                                                                 min = 0, max = 10, value = 0)
                                                   ), 
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'bisulfite'",
                                                     sliderInput("i_minCG_bis", h4("Minimum Number of CGs per Amplicon"),
                                                                 min = 0, max = 10, value = 5)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'NOME'",
                                                     sliderInput("i_minCG_NOME", h4("Minimum Number of CGs per Amplicon"),
                                                                 min = 0, max = 10, value = 5)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'CLEVER'",
                                                     sliderInput("i_minCG_clever", h4("Minimum Number of CGs per Amplicon"),
                                                                 min = 0, max = 10, value = 5)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_bisulfite'",
                                                     sliderInput("i_minCG_hp_bis", h4("Minimum Number of CGs per Amplicon"),
                                                                 min = 0, max = 10, value = 1)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_NOME'",
                                                     sliderInput("i_minCG_hp_NOME", h4("Minimum Number of CGs per Amplicon"),
                                                                 min = 0, max = 10, value = 1)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_CLEVER'",
                                                     sliderInput("i_minCG_hp_clever", h4("Minimum Number of CGs per Amplicon"),
                                                                 min = 0, max = 10, value = 1)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'CrispRCas9PCR'",
                                                     sliderInput("i_minCG_crispr", h4("Minimum Number of CGs per Amplicon"),
                                                                 min = 0, max = 10, value = 1)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type != 'genomic' && input.i_primer_type != 'CrispRCas9PCR'",
                                                     sliderInput("i_minC2T", h4("Minimum 'C' to 'T' Conversions in Forward Primer"),
                                                                 min = 0, max = 10, value = 3)
                                                   ),
                                                   conditionalPanel(
                                                     condition="input.i_primer_type != 'genomic' && input.i_primer_type != 'CrispRCas9PCR'",
                                                     sliderInput("i_minG2A", h4("Minimum 'G' to 'A' Conversions in Reverse Primer"),
                                                                 min = 0, max = 10, value = 3)
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'hp_bisulfite' || input.i_primer_type == 'hp_NOME' || input.i_primer_type == 'hp_CLEVER'", 
                                                     sliderInput("i_hp.length", label = h4("Length of One Arm in the Hairpin Molecule"),
                                                                 min = 0, max = 1000, value = c(50, 200))
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.i_primer_type == 'genomic'",
                                                     sliderInput("i_chop.size", label = h4("Input Sequence Slicing"),
                                                                 min = 0, max = 50, value = 30)
                                                     
                                                   )
                                            )
                                            ),
                                            hr(),
                                            h2("Add Adapters to my Primers: "),
                                            # adapter box
                                            fluidRow(column(12,
                                                            textInput("adapterForward", h4("Add a Specific Sequence to the 5' End of the Forward Primer: "), ""), #CTTTCCCTACACGACGCTCTTCCGATCT
                                                            textInput("adapterReverse", h4("Add a Specific Sequence to the 5' End of the Reverse Primer: "), "") #GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
                                            )),
                                            hr(),
                                            h2("Further Settings: "),
                                            fluidRow(column(12,
                                                            checkboxInput("i_remove.primers.with.n", label = h4("Remove Primers That Contain N Bases"), TRUE),
                                                            checkboxInput("i_allow.repeats.in.primers", label = h4("Allow Repeats in Primers"), FALSE),
                                                            checkboxInput("i_allow.repeats.in.amplicon", label = h4("Allow Repeats in Amplicons"), FALSE)
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
                                                         style="color: #fff; background-color: #3c8dbc; border-color: #337ab7; padding:25px; font-size:200%; width:1000px; margin-left:100px; margin-right:0px"),
                                            bsTooltip("action", "The Computation of Your Primers May Take a Few Minutes, Please Wait Until You Receive a Notification That Your Primers Are Finished.", "left")
                                          )
                                          
                                                )
                                          
                                          ),
                                  
                                  tabItem(tabName = "PDresults",
                                          fluidRow(
                                            tabBox(title = "",
                                                   id = "PDresultsTabbox",
                                                   width = 12,
                                                   tabPanel("Overview",
                                                            withSpinner(DT::dataTableOutput("viewprimerdesigns"), type = 1, color = "#3c8dbc", size = 2)
                                                   ),
                                                   tabPanel("Toplist",
                                                            withSpinner(DT::dataTableOutput("viewtoplist"), type = 1, color = "#3c8dbc", size = 2)
                                                   ),
                                                   tabPanel("Complete Primer List",
                                                            withSpinner(DT::dataTableOutput("viewwholelist"), type = 1, color = "#3c8dbc", size = 2)
                                                   ),
                                                   tabPanel("Whitelist",
                                                            withSpinner(DT::dataTableOutput("viewwhitelist"), type = 1, color = "#3c8dbc", size = 2)
                                                   ),
                                                   tabPanel("Blacklist",
                                                            withSpinner(DT::dataTableOutput("viewblacklist"), type = 1, color = "#3c8dbc", size = 2)
                                                   ),
                                                   tabPanel("Logfile",
                                                            withSpinner(DT::dataTableOutput("viewlogfile"), type = 1, color = "#3c8dbc", size = 2)
                                                   ),
                                                   tabPanel("Settings",
                                                            withSpinner(DT::dataTableOutput("viewsettings"), type = 1, color = "#3c8dbc", size = 2)
                                                   ),
                                                   tabPanel("Summary",
                                                            withSpinner(DT::dataTableOutput("viewSummary"), type = 1, color = "#3c8dbc", size = 2)
                                                   )
                                            )
                                          ),
                                          fluidRow(
                                            box(title = "Selected List",
                                                width = 12,
                                                helpText("You Can Download Your Primer Pairs or Hand Them Over for ePCR by Adding Them to Your List of Selected Primers Below. To Add Primers to Your Selected List, Mark Them on the Complete Primer List Above!"),
                                                hr(),
                                                DT::dataTableOutput("viewSelectlist"),
                                                hr(),
                                                downloadButton('downloadSelectedPrimers', 'Selected Primers')
                                            )
                                          ),
                                          box(withSpinner(uiOutput("plot3"), type = 1, color = "#3c8dbc", size = 2),
                                              title = h2("Graphs"),
                                              status = "primary",
                                              solidHeader = TRUE,
                                              collapsible = TRUE,
                                              width = 12
                                          )
                                          
                                  ),
                                  tabItem(tabName = "PrimerQC",
                                          box(title = h2("Upload Files for ePCR"),
                                              status = "primary", 
                                              solidHeader = TRUE,
                                              width = 6,
                                              helpText("To Check the Quality of Your Primer Pairs, We Perform BLAST Search Against the Reference Genome and Provide Potential Fragments During PCR.", align = "left"),
                                              hr(),
                                              h5(strong("Upload Your Own Primers for ePCR:"), align = "left"),
                                              helpText("(.fasta File Format is Required, We Expect Your Primers to Be in 5'-3' Orientation)", align = "left"),
                                              fileInput("Fprimers", "Upload Forward Primers", multiple = TRUE, accept = ".fasta"), 
                                              DT::dataTableOutput("forward.primers"),
                                              fileInput("Rprimers", "Upload Reverse Primers", multiple = TRUE, accept = ".fasta"), 
                                              DT::dataTableOutput("reverse.primers"),
                                              hr(),
                                              h3(strong("OR"), align = "center"),
                                              hr(),
                                              h5("Import Primer Pairs From Primer Design That Were Added to Your Selected List:", align = "center"),
                                              actionButton("loadprimers", "Import Primers", icon = icon("file-import"), align="center",
                                                           style="margin-left:225px; margin-right:0px"),
                                              bsTooltip("loadprimers", "Import Primers You Added to the Selected List", "bottom", "hover")
                                          ),
                                          box(title = h2("Settings for ePCR"),
                                              status = "primary", 
                                              solidHeader = TRUE,
                                              width = 6,
                                              textInput("blast_id","Name Your Primerblast:", paste0("MyBlast")),
                                              selectInput("genome", "Reference Genome for ePCR: ",choices = gsub("UCSC.", "", gsub("BSgenome.", "", c("BSgenome.Hsapiens.UCSC.hg18", "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm9", "BSgenome.Mmusculus.UCSC.mm10")))), #gsub("UCSC.", "", gsub("BSgenome.", "", installed.genomes()))), #(names(a) = strsplit(installed.genomes(), ".", fixed=TRUE)[1][4])),
                                              bsTooltip("genome", "Select the Genome Against Which You Want to BLAST Your Primer Pairs!", "top", "hover"),
                                              checkboxInput("is_bisulfite", h5("These Are Bisulfite Primers!"), FALSE),
                                              bsTooltip("is_bisulfite", "Check This Box if Your Primer Pairs Were Designed for a Bisulfite Converted Reference Genome", "top", "hover"),
                                              sliderInput("gap", "Maximum Size of Reported Product", min = 500, max = 10000, value = 2000, step = 500),
                                              sliderInput("primer_mismatches", "Number of Mismatches Allowed in Primer BLAST", min=0, max=25, value=7),
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
                                                       style="color: #fff; background-color: #3c8dbc; border-color: #337ab7; padding:25px; font-size:200%; width:1000px; margin-left:100px; margin-right:0px"),
                                          bsTooltip("compute_ePCR", "The Computation of Your ePCR May Take a Few Minutes, Please Wait Until You Receive a Notification That Computation Has Finished.", "left", "hover")
                                  
                                  ),
                                  tabItem(tabName = "PrimerQCResults",
                                          fluidRow(
                                            tabBox(title = "",
                                                   id = "PQCTabbox",
                                                   width = 12,
                                                   tabPanel("Results",
                                                            #actionButton("refreshPQC", label = "Results", icon = icon("sync-alt")),
                                                            #bsTooltip("refreshPQC", "Please Press This Button to View the Results of Your ePCR!", "right", "hover"),
                                                            selectInput("test_select", label = "Filter Results by Primerpair: ", choices ="", multiple = FALSE),
                                                            tags$div(id="txt_for_selector", helpText("")),
                                                            withSpinner(uiOutput("pQC.results"), type = 1, color = "#3c8dbc", size = 2)
                                                   ),
                                                   tabPanel("Overview",
                                                            withSpinner(DT::dataTableOutput("ePCR.overview"), type = 1, color = "#3c8dbc", size = 2)
                                                            ),
                                                   tabPanel("Settings",
                                                            withSpinner(DT::dataTableOutput("ePCR.settings"), type = 1, color = "#3c8dbc", size = 2)
                                                            ),
                                                   tabPanel("Logfile",
                                                            withSpinner(DT::dataTableOutput("ePCR.logfile"), type = 1, color = "#3c8dbc", size = 2)
                                                            ),
                                                   tabPanel("Summary",
                                                            withSpinner(DT::dataTableOutput("pQC.summary"), type = 1, color = "#3c8dbc", size = 2)
                                                   )#,
                                                   #hr()
                                                   #downloadButton('downloadePCRresults', 'Results')
                                            )
                                            ),
                                            
                                            fluidRow(box(#actionButton("graphics_ePCR", label = "Graphs", icon = icon("sync-alt")),
                                              withSpinner(uiOutput("ePCR_Graphs"), type = 1, color = "#3c8dbc", size = 2),
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
                                              helpText("This Website Was Created by the Epigenetics Department of Saarland University."),
                                              tags$a(href="http://epigenetik.uni-saarland.de/en/home/", "Visit Our Homepage Here."),
                                              hr(),
                                              helpText("In Case You Encounter Any Problems During Use of This Website, Please Contact Us!")
                                          )
                                  )
                                  
                                )
                          )
                      )
                )