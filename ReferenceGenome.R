library(BSgenome)
library(methods)
library(devtools)
library(rBLAST)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)


setClass("referencegenome",
         #contains="BSgenome", #--> use this, if inheritance of BSgenome is indeed needed
         slots = list(genome="BSgenome",
                      name="character",
                      wd="character"))

#setting generic is not needed here, we are overwriting an already existing function
setMethod("show", signature(object="referencegenome"), function(object){
  return (object@name)
})

setGeneric(
  "setReferenceGenome",
  function(x){
    standardGeneric("setReferenceGenome")
  }
)

setMethod("setReferenceGenome", "referencegenome", function(x){
  cat("set the referencegenome here. This is to be implemented. ")
})

setGeneric(
  "getReferenceGenome",
  function(y){
    standardGeneric("getReferenceGenome")
  }
)

setMethod("getReferenceGenome", "referencegenome", function(y){
  return (y@genome)
})

setGeneric(
  "getBlastDB",
  function(z, is_bisulfite){
    standardGeneric("getBlastDB")
  }
)

setMethod("getBlastDB", "referencegenome", function(z, is_bisulfite){

  # get database for bisulfite primers 
  if (is_bisulfite == TRUE){
    print("bisulfite primer quality control")

    #now check if makeblastdb is needed for this genome
    # get name of desired reference geomoe 
    name <- z@name
    
    # now get CT and GA converted reference genomes from the database folder
    # reference genomes have to be created by converting non-bisulfite converted versions of the reference genome
    # this can be done via the server: linux commands (can be entered via commandline) #name# is the corresponding file name, depending on the job run
    # 1.a) create C to T converted reference genome: "cat #path_to_fastafile_of_reference_genome_non_bis# | tr 'cC' 'tT' > #name_for_new_CT_converted_genome.fa#"
    # 1.b) to reverse conversion of the "chr" at the beginning of the fastafile to "thr", the following command needs to be entered: "sed -i 's/>t/>c/g'"
    # 2. create G to A converted reference genome: "cat #path_to_fastafile_of_reference_genome_non_bis# | tr 'gG' 'aA' > #name_for_new_GA_converted_genome.fa#"
    # alternatively, a program to achieve the same without command line could be written (this could be integrated into EpiPrimer to enable users to upload ANY reference genome of their own.)
    
    bis_name <- paste0("Bis_", gsub("UCSC.", "", gsub("BSgenome.", "", name)))
    bis_path_CT <- file.path(getwd(), "database", bis_name, "CTgenome.fa", fsep = .Platform$file.sep)
    bis_path_GA <- file.path(getwd(), "database", bis_name, "GAgenome.fa", fsep = .Platform$file.sep)
    bis_path_CT_nin <- file.path(getwd(), "database", bis_name, "CTgenome.nin", fsep = .Platform$file.sep)
    bis_path_GA_nin <- file.path(getwd(), "database", bis_name, "GAgenome.nin", fsep = .Platform$file.sep)
    bis_path_CT_nsq <- file.path(getwd(), "database", bis_name, "CTgenome.nsq", fsep = .Platform$file.sep)
    bis_path_GA_nsq <- file.path(getwd(), "database", bis_name, "GAgenome.nsqr", fsep = .Platform$file.sep)
    print(bis_path_CT)
    print(bis_path_GA)
    
    # check if a wrong database is used.
    if (!file.exists(bis_path_CT) || !file.exists(bis_path_GA)){
      stop("Invalid referencegenome - bisulfite reference genome not found!")
    }
    
    CTdb <- (
      if(!file.exists(bis_path_CT)
         || !file.exists(bis_path_CT_nin)
         || !file.exists(bis_path_CT_nsq)
      ){
        # build the database and index for ref genome
        print("making new database...")
        makeblastdb(file=bis_path_CT, dbtype = "nucl")
        blast(bis_path_CT)
      }
      else{
        #files already exist so no need to build the database again
        blast(bis_path_CT)
      }
    )
    
    GAdb <- (
      if(!file.exists(bis_path_GA)
         || !file.exists(bis_path_GA_nin)
         || !file.exists(bis_path_GA_nsq)
      ){
        # build the database and index for ref genome
        print("making new database...")
        makeblastdb(file=bis_path_GA, dbtype = "nucl")
        blast(bis_path_GA)
      }
      else{
        #files already exist so no need to build the database again
        blast(bis_path_GA)
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

setMethod("getAssemblyName", "referencegenome", function(x){
  if (x@name == "BSgenome.Hsapiens.UCSC.hg18"){
    return ("hg18")
  }
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
  return ("The genome assembly you choose is invalid!")  
})