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
  function(z){
    standardGeneric("getBlastDB")
  }
)

setMethod("getBlastDB", "ReferenceGenome", function(z){

  print("checking ref genome wd")
  print(getwd())
  #create a folder, if there is not yet one for this genome - this line throws a warning, if the dir exists already
  dir.create(file.path(z@wd, fsep=.Platform$file.sep), recursive = TRUE)
  #get the fasta file for the corresponding genome
  fastafile <- export(getReferenceGenome(z), (file.path(z@wd, paste(z@name, ".fasta", sep=""), fsep=.Platform$file.sep)), compress="no", compression_level=NA, verbose=TRUE)
  #now check if makeblastdb is needed for this genome
  CTdb <- (
    if(!file.exists(paste(z@wd, z@name, ".fa.nhr", sep=""))
       || !file.exists(paste(z@wd, z@name, ".fa.nin", sep=""))
       || !file.exists(paste(z@wd, z@name, ".fa.nsq", sep=""))
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

  GAdb <- (
    if(!file.exists(paste(z@wd, z@name, ".fa.nhr", sep=""))
       || !file.exists(paste(z@wd, z@name, ".fa.nin", sep=""))
       || !file.exists(paste(z@wd, z@name, ".fa.nsq", sep=""))
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
  dbList <- list("CTdb"= CTdb, "GAdb"= GAdb)
  return (dbList)

})
