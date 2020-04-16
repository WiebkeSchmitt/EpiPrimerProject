# EpiPrimerProject


Description of the project: 

Table of Contents:

# 1. 	Installation

Requirements for local development of the project: 
- RStudio and R version 3.6.0 installed
- Git installed

## 1.1 	Using the Webapplication

A web application will be available shortly.

## 1.2	Running the project locally

Clone the code for the epiprimer project from its GitHub repository: 

https://github.com/WiebkeSchmitt/EpiPrimerProject.git

Additionally, you need to create two new folders: 
1. A folder called "database"
2. A folder called "flowcell_package" (this folder should be given to you by an administrator of the project)
3. A folder called "PrimerQC" to hold results from Primer Quality Control

R packages, that need to be installed: 

//TODO

You can now open RStudio and herein open the project folder. Pushing the button "Run app" when opening the file ui_shinydahsboard.R will run the project.

#	2. 	Usage

##	2.1	Primer Design Pipeline

The pipeline for primerdesign is contained in the file primer.design.pipeline_jil_v3.0_standalone. 
This pipeline can be run seperately via the commandline if desired.

Primers can be designed for mouse genome assemblies 8 and 9 as well as for human genome assemblies 18 and 19. 

##	2.2	Primer Blast (formerly: Primer Quality Control)

The Primer Blast allows (as already contained in the name) the blast of user designed primers and primer created during the primer design against a reference genome. 

We expect primers to be given in 5' to 3' direction. 

We support blasting for primers against a bisulfite treated human genome.  

//TODO

# 3. 	Credits

This app was created at the Epigenetics department of Saarland University:
[Epigenetics Saarland University](http://epigenetik.uni-saarland.de/en/home/ "Epigenetics Homepage")

Questions concerning the project may be adressed to Dr. Gilles Gasparoni. 
