# EpiPrimerProject


Description of the project: 

Table of Contents:

# 1. 	Installation

Requirements for local development of the project: 
- RStudio and R version 3.6.0 installed
- Git installed

## 1.1 	Using the Webapplication

A web application will be available shortly.

## 1.2	Running the project locally for development

Clone the code for the epiprimer project from its GitHub repository: 

https://github.com/WiebkeSchmitt/EpiPrimerProject.git

Additionally, you need two folders: 
1. A folder called "database" containing the reference genomes needed for ePCR
2. A folder called "PrimerQC" to hold results from the ePCR tool

R packages, that need to be installed: 

//TODO

You can now open RStudio and herein open the project folder. Pushing the button "Run app" when opening the file ui_shinydahsboard.R will run the project.

#	2. 	Usage

##	2.1	Primer Design Pipeline

The pipeline for primerdesign is contained in the file primer.design.pipeline_jil_v3.0_standalone. 
This pipeline can be run seperately via the commandline if desired.

Primers can be designed for mouse genome assemblies 8 and 9 as well as for human genome assemblies 18 and 19. 

##	2.2	ePCR (formerly: Primer Quality Control, Primer Blast)

The ePCR module allows the blast of user designed primers and primers created during the primer design against a reference genome. 

We expect primers to be given in 5' to 3' direction. 

We support this for: hg18, hg19, hg38, mm9, mm10

//TODO

# 3. 	Credits

This app was created at the Epigenetics department of Saarland University:
[Epigenetics Saarland University](http://epigenetik.uni-saarland.de/en/home/ "Epigenetics Homepage")

Questions concerning the project may be adressed to Dr. Gilles Gasparoni. 
