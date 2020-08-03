# EpiPrimerProject

Description of the project: 

Table of Contents:

# 1. 	Using the Webapplication

A web application is currently available on the T7600 precision workstation and can be reached via: http://134.96.17.175:8888/
The usage of the web application requires access to the T7600 server. After a session on the server is started, the web application can be reached via the above URL.

## 1.1 	Installations required for local development

Requirements for local development of the project: 
- RStudio (or similar development environment for R) and R version 3.6.0 installed
- Git installed

## 1.2	Local folder structure required for development

The code of the epiprimer project may be clonded from its GitHub repository: 

https://github.com/WiebkeSchmitt/EpiPrimerProject.git

Additionally, two folders are required: 
1. A folder called "database" containing the reference genomes needed for ePCR. 
2. A folder called "ePCR" to hold results from the ePCR tool

Overview of the required folder structure: 

.
+-- _server.R
+-- _ui.R
+-- _generalDesign.R
+-- _HelperFunctions.R 
+-- _primer.design.R 
+-- _Primerpair.R 
+-- _ReferenceGenome.R 
+-- _ePCR
+-- _database
|   +-- Bis_Hsapiens.hg18
	|	+--	CTgenome.fa
	|	+--	CTgenome.fa.nhr
	|	+--	CTgenome.fa.nin
	|	+--	CTgenome.fa.nsq
	|	+--	GAgenome.fa	
	|	+--	GAgenome.fa.nhr
	|	+--	GAgenome.fa.nin
	|	+--	GAgenome.fa.nsq
|   +-- Bis_Hsapiens.hg19
	|	+--	CTgenome.fa
	|	+--	CTgenome.fa.nhr
	|	+--	CTgenome.fa.nin
	|	+--	CTgenome.fa.nsq
	|	+--	GAgenome.fa	
	|	+--	GAgenome.fa.nhr
	|	+--	GAgenome.fa.nin
	|	+--	GAgenome.fa.nsq
|	+-- Bis_Hsapiens.hg38
	|	+--	CTgenome.fa
	|	+--	CTgenome.fa.nhr
	|	+--	CTgenome.fa.nin
	|	+--	CTgenome.fa.nsq
	|	+--	GAgenome.fa	
	|	+--	GAgenome.fa.nhr
	|	+--	GAgenome.fa.nin
	|	+--	GAgenome.fa.nsq
|	+-- Bis_Mmusculus.mm9
	|	+--	CTgenome.fa
	|	+--	CTgenome.fa.nhr
	|	+--	CTgenome.fa.nin
	|	+--	CTgenome.fa.nsq
	|	+--	GAgenome.fa	
	|	+--	GAgenome.fa.nhr
	|	+--	GAgenome.fa.nin
	|	+--	GAgenome.fa.nsq
|	+-- Bis_Mmusculus.mm10
	|	+--	CTgenome.fa
	|	+--	CTgenome.fa.nhr
	|	+--	CTgenome.fa.nin
	|	+--	CTgenome.fa.nsq
	|	+--	GAgenome.fa	
	|	+--	GAgenome.fa.nhr
	|	+--	GAgenome.fa.nin
	|	+--	GAgenome.fa.nsq
|	+-- BSgenome.Hsapiens.UCSC.hg18
	|	+--	BSgenome.Hsapiens.UCSC.hg18.fasta
	|	+--	BSgenome.Hsapiens.UCSC.hg18.fasta.nhr
	|	+--	BSgenome.Hsapiens.UCSC.hg18.fasta.nin
	|	+--	BSgenome.Hsapiens.UCSC.hg18.fasta.nsq	
|	+-- BSgenome.Hsapiens.UCSC.hg19
	|	+--	BSgenome.Hsapiens.UCSC.hg19.fasta
	|	+--	BSgenome.Hsapiens.UCSC.hg19.fasta.nhr
	|	+--	BSgenome.Hsapiens.UCSC.hg19.fasta.nin
	|	+--	BSgenome.Hsapiens.UCSC.hg19.fasta.nsq	
|	+-- BSgenome.Hsapiens.UCSC.hg38
	|	+--	BSgenome.Hsapiens.UCSC.hg38.fasta
	|	+--	BSgenome.Hsapiens.UCSC.hg38.fasta.nhr
	|	+--	BSgenome.Hsapiens.UCSC.hg38.fasta.nin
	|	+--	BSgenome.Hsapiens.UCSC.hg38.fasta.nsq	
|	+-- BSgenome.Mmusculus.UCSC.mm9
	|	+--	BSgenome.Hsapiens.UCSC.mm9.fasta
	|	+--	BSgenome.Hsapiens.UCSC.mm9.fasta.nhr
	|	+--	BSgenome.Hsapiens.UCSC.mm9.fasta.nin
	|	+--	BSgenome.Hsapiens.UCSC.mm9.fasta.nsq	
|	+-- BSgenome.Mmusculus.UCSC.mm10
	|	+--	BSgenome.Hsapiens.UCSC.mm10.fasta
	|	+--	BSgenome.Hsapiens.UCSC.mm10.fasta.nhr
	|	+--	BSgenome.Hsapiens.UCSC.mm10.fasta.nin
	|	+--	BSgenome.Hsapiens.UCSC.mm10.fasta.nsq	
+-- _www
|   +-- custom.css
+-- _test

This list does not contain all files contained in the git repository, but all files minimally required to run EpiPrimer.
The files from the BSGenome package required for ePCR of genomic primer pairs (.nhr, .nin, .nsq) will be created automatically, when the project runs an ePCR for the genome for the first time. 
The BSGenome for the required organism and assembly has to first be installed manually via the RStudio command line in order to be used: 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")

Organism and assembly need to be adjusted to the organism and assembly in question.
The folders representing the bisulfite converted reference genomes may be adapted from EpiPrimer's folder on the t-7600 server. 

R packages, that need to be installed to run EpiPrimer locally: (This may not be a complete list)
- shiny
- shinydashboard
- devtools
- rBLAST
- BSGenome
- Biostrings
- tidyr
- dplyr
- shinybs
- tippy
- xml2
- httr
- stringr
- ggplot2

You can now open RStudio and herein open the project folder. Pushing the button "Run app" when opening the file ui.R or server.R will run EpiPrimer locally.

## 1.3. 	Live version folder structure on the t-7600 server

For working with the live version of EpiPrimer, access to the folder on the t-7600 server is needed. 

The folder of the live version of EpiPrimer is located here: 
/projects/epiprimer
The code for the live version is here: 
/projects/epiprimer/epiprimer_git/EpiPrimerProject
Logfiles of the live version required for debugging may be found here: 
/projects/epiprimer/epiprimer_git/log

To publish a new version of EpiPrimer on the t-7600 server contact Michael Scherer.

#	2. 	Usage

##	2.1	Primer Design Pipeline

The original pipeline for primerdesign is contained in the file primer.design.pipeline_jil_v3.0_standalone.
This pipeline can be run seperately via the commandline if desired, e. g.: 
primer.design.pipeline("./input_single.txt", "regions", "./folder_name", "job_name", "bisulfite", "fast", "top", 23, 34, 48, 60, 2, TRUE, 7, TRUE, 150, 500, 0, 5, FALSE, 12, 3, 3, TRUE, 0.01, 0, 20, 0, 0, 0, 0, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, NA, 50, 300, 30, "", "", TRUE)
Primers can be designed for mouse genome assemblies 8 and 9 as well as for human genome assemblies 18, 19 and 38. 

This primer design pipeline was refactored and is now contained in several files: 
_generalDesign.R
_HelperFunctions.R 
_primer.design.R 
_Primerpair.R 

During the primer design process, the generalDesign.R script designs primer pairs
according to the user input passed on by the user interface of EpiPrimer. The script
calls the required function from the script primer.design.R, which then first performs
common primer computations using its function "generalDesign()" and then further
primer type specific operations. The resulting primer pairs are stored as Primerpair
objects and returned as data structure dataframe to the generalDesign.R script,
which then prepares the primer pairs as output for the user interface.

##	2.2	ePCR (formerly: Primer Quality Control, Primer Blast)

The ePCR module allows the blast of user designed primers and primers created during the primer design against a reference genome and calculates potential PCR fragments resulting from these potential bindings sites.
EpiPrimer expect primers to be given in 5' to 3' direction. 

We support ePCR for: hg18, hg19, hg38, mm9, mm10

###	2.2.1	Creation of Bisulfite Converted Reference Genomes

In case EpiPrimer should be extended by further bisulfite converted reference genomes, the following commands may be used to create them: 

Assume, the genome to be converted is saved as a .fa file in location:
"./Genome Name folder/Genome.fa"

To create the C-to-T converted reference genome, all Cs present in the genome will
have to be converted to Ts. To do so, the following command may be used from the
folder "Genome Name":

cat Genome.fa j tr 'cC' 'tT' > Genome CTgenome.fa

Since this command will also convert the Cs contained in the keyword "chr" in the
.fa file to Ts, this step has to be reversed afterwards in order to maintain a valid .fa
file. To do so, the following command can be used:

sed -i 's/>t/>c/g' Genome CTgenome.fa

To create the G-to-A converted reference genome, all Gs present in the genome will
have to be converted to As. Similarly to the creation of the C-to-T converted reference
genome, the following command may be used from the folder "Genome Name":

cat Genome.fa j tr 'gG' 'aA' > Genome GAgenome.fa

No reversal steps are necessary when creating G-to-A converted reference genomes.
File and folder names need to be adjusted according to the individual reference
genome used.

# 3. 	Credits

This app was created at the Epigenetics department of Saarland University:
[Epigenetics Saarland University](http://epigenetik.uni-saarland.de/en/home/ "Epigenetics Homepage")

Questions concerning the project may be adressed to Dr. Gilles Gasparoni. 
