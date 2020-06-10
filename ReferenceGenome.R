library(BSgenome)
library(methods)
library(devtools)
library(rBLAST)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)


setClass("ReferenceGenome",
         #contains="BSgenome", #--> use this, if inheritance of BSgenome is indeed needed
         slots = list(genome="BSgenome",
                      name="character",
                      wd="character"))

#setting generic is not needed here, we are overwriting an already existing function
setMethod("show", signature(object="ReferenceGenome"), function(object){
  return (object@name)
})

setGeneric(
  "setReferenceGenome",
  function(x){
    standardGeneric("setReferenceGenome")
  }
)

setMethod("setReferenceGenome", "ReferenceGenome", function(x){
  cat("set the referencegenome here. This is to be implemented. ")
})

setGeneric(
  "getReferenceGenome",
  function(y){
    standardGeneric("getReferenceGenome")
  }
)

setMethod("getReferenceGenome", "ReferenceGenome", function(y){
  return (y@genome)
})

setGeneric(
  "getBlastDB",
  function(z, is_bisulfite){
    standardGeneric("getBlastDB")
  }
)

setMethod("getBlastDB", "ReferenceGenome", function(z, is_bisulfite){

  # get database for bisulfite primers 
  if (is_bisulfite == TRUE){
    print("bisulfite primer quality control")

    #now check if makeblastdb is needed for this genome
    # TODO: Adjust this for multiple reference genomes!
    CTdb <- (
      if(!file.exists(paste("C:/Users/Wiebk/Desktop/epiprimer/database", "/", "Bis_hg19", "/", "CTgenome", ".fa.nhr", sep=""))
         || !file.exists(paste("C:/Users/Wiebk/Desktop/epiprimer/database", "/", "Bis_hg19", "/", "CTgenome", ".fa.nin", sep=""))
         || !file.exists(paste("C:/Users/Wiebk/Desktop/epiprimer/database", "/", "Bis_hg19", "/", "CTgenome", ".fa.nsq", sep=""))
      ){
        # build the database and index for ref genome
        makeblastdb(file=(file.path(paste("C:/Users/Wiebk/Desktop/epiprimer/database", "/", "Bis_hg19", "/", "CTgenome", ".fa", sep=""))), dbtype = "nucl")
        blast(file.path(paste("C:/Users/Wiebk/Desktop/epiprimer/database", "/", "Bis_hg19", "/", "CTgenome", ".fa", sep="")))
      }
      else{
        #files already exist so no need to build the database again
        blast(file.path(paste("C:/Users/Wiebk/Desktop/epiprimer/database", "/", "Bis_hg19", "/", "CTgenome", ".fa", sep="")))
      }
    )
    
    GAdb <- (
      if(!file.exists(paste("C:/Users/Wiebk/Desktop/epiprimer/database", "/", "Bis_hg19", "/", "GAgenome", ".fa.nhr", sep=""))
         || !file.exists(paste("C:/Users/Wiebk/Desktop/epiprimer/database", "/", "Bis_hg19", "/", "GAgenome", ".fa.nin", sep=""))
         || !file.exists(paste("C:/Users/Wiebk/Desktop/epiprimer/database", "/", "Bis_hg19", "/", "GAgenome", ".fa.nsq", sep=""))
      ){
        # build the database and index for ref genome
        makeblastdb(file=(file.path(paste("C:/Users/Wiebk/Desktop/epiprimer/database", "/", "Bis_hg19", "/", "GAgenome", ".fa", sep=""))), dbtype = "nucl")
        blast(file.path(paste("C:/Users/Wiebk/Desktop/epiprimer/database", "/", "Bis_hg19", "/", "GAgenome", ".fa", sep = "")))
      }
      else{
        #files already exist so no need to build the database again
        blast(file.path(paste("C:/Users/Wiebk/Desktop/epiprimer/database", "/", "Bis_hg19", "/", "GAgenome", ".fa", sep = "")))
      }
    )
    
    #return CTdb and GAdb for further calculations in server.R as a list of CTdb and GAdb
    dbList <- list("CTdb"= CTdb, "GAdb"= GAdb)
    return (dbList)
    
  } else {
    
    #normal primer blast 
    
    print("non-bisulfite primer reference genome is created")
    
    #create a folder, if there is not yet one for this genome - this line throws a warning, if the dir exists already
    if(!dir.exists(paste(z@wd, fsep=.Platform$file.sep))){
      dir.create(file.path(z@wd, "database", fsep=.Platform$file.sep), recursive = TRUE)
      #get the fasta file for the corresponding genome
      fastafile <- export(getReferenceGenome(z), (file.path(z@wd, paste(z@name, ".fasta", sep=""), fsep=.Platform$file.sep)), compress="no", compression_level=NA, verbose=TRUE)
    }
    
    #now check if makeblastdb is needed for this genome
    #print(file.exists(paste(z@wd, "/", z@name, ".fasta", sep="")))
    genomeDB <- (
      if(!file.exists(paste(z@wd, "/", z@name, ".fasta", sep=""))
         || !file.exists(paste(z@wd, "/", z@name, ".fasta.nin", sep=""))
         || !file.exists(paste(z@wd, "/", z@name, ".fasta.nsq", sep=""))
      ){
        # build the database and index for ref genome
        makeblastdb(file=(file.path(z@wd, paste( z@name, ".fasta", sep=""), fsep=.Platform$file.sep)), dbtype = "nucl")
        blast(file.path(z@wd, paste(z@name, ".fasta", sep=""), fsep=.Platform$file.sep))
      }
      else{
        #files already exist so no need to build the database again
        blast(file.path(z@wd, paste(z@name, ".fasta", sep=""), fsep=.Platform$file.sep))
      }
    )
    
    #return CTdb and GAdb for further calculations in server.R as a list of CTdb and GAdb
    dbList <- list("genomeDB"= genomeDB)
    print("finished computation of non-bisulfite reference genome")
    return (dbList)
  }
  
})

setGeneric(
  "getAssemblyName",
  function(x){
    standardGeneric("getAssemblyName")
  }
)

setMethod("getAssemblyName", "ReferenceGenome", function(x){
  if (x@name == "BSgenome.Hsapiens.UCSC.hg19"){
    return ("hg19")
  }
  if (x@name == "BSgenome.Hsapiens.UCSC.hg38"){
    return ("hg38")
  }
  if (x@name == "BSgenome.Mmusculus.UCSC.mm10"){
    return ("mm10")
  }
  if (x@name == "BSgenome.Mmusculus.UCSC.mm9"){
    return ("mm9")
  }  
  return ("The genme assembly you choose is invalid!")  
})