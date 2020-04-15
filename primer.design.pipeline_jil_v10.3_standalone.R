#fixed creating graph folder when it's not chosen 

###################################################################################
#structure of input files.
#region.file: tab.del with header: 
#"chr","start","end","assembly",""sequenceID","target.start","target.end" ('target' cols are optional)
# "chr": - e.g. "chr1"
# "start"/"end": 12345 / 12349
# "assembly: "hg18", "hg19", "mm9","mm10"...
# 
#sequence.file: tab.del with header: "sequence","sequenceID","target.start","target.end" ('target' cols are optional)
#RE.file : tab.del with header: "RE.id","RE.seq","linker.seq"
####################################################################################################
#v1.1:  [improvement of NOME.primer.design.script_v1.0]
#       - html links are working
#       - option to check for primer SNPs
#       - results are restructured 
#       - faster calculation in amplicon size filtering section.
#v1.2:  - increased performance in oligo self alignment calculations 
#       - (see improvement of NOME primer design script-->now v1.1).
#v1.3:  - creates a primer toplist table in addition (good for large projects).
#v1.4:  - now you can use either .bed files (regions) or .txt files as input.
#v2.0:  - standard bisulfite primer design is now enabled.
#v2.4:  - SNP analysis is improved.
#v2.5:  - SNPs analysis is even faster now (new get.snp function is used).
#       - Report gives a more detailed output.
#       - Toplist with best primers for each amplicon is produced.
#       - Several Bugs are removed.
#v2.10: - Bug fixed for output of amplicon sequences.
#v2.11: - Bug fixed for output of amplicon sequences. (ONLY Checked for BisulfitePrimers)
#ATTENTION: Have to check again the primer design for "minus" strand...but "plus" works fine...
#v3.00: - Also compatible with hairpin Bisulfite PrimerDesign.
#       - Allows to force a specific sub-region to be present in PCR product (target region).
#v3.01: - A logging file is now enabled.
#v3.04: - Calculates GC contents for Primers and amplicons.
#v4.0:  - Support of Hairpin-Bisulfite and Hairpin-NOME PrimerDesign.
#v4.1:  - More output tables for hp analyses.
#v5.0p2:  - Enabled Parallelization of PrimerDesign for considerable Speed-up on multi-core CPUs.
#v6.0:  - Enabled Removal of low complexity primers (repetitive 2-mers and 3-mers).
#v6.1:  - Working on the SNP analysis part, not yet finished.
#v6.2:  - MAJOR BUG removed: chr column in output files were wrong due to unwanted factor conversion. CHECK AGAIN!
#v6.3:  - Completely fixed bug from previous version 6.2
#v6.4:  - New feature: to the designed primers, one can add NGS adaptors.
#v6.5:  - Remove Parallelization...too much RAM issues...
#v6.6:  - Bugfix: target regions can now be specified individually for each input sequence.
#v6.7:  - Bugfix: for input.type= regions, the absolute position of amplicon should be calculated in a column.
#       - New Feature: oin utput file a column with ready to use (for ucsc) adress of amplicon will be generated.
#v6.8:  - Work on enabling analysis of top/bottom strand. Also rename "plus" to "top" and "minus" to "bottom".
#v6.9:  - to the input.sequence matrix add one column for the converted version of input sequences.
#v6.10: - Stream lined code for more efficient work.
#v7.0:  - Major changes in subscripts. Run performance increased >2-fold by more efficient code for primer.design parts
#       - Major bugfix: some regions were removed in earlier versions without reason. Now ther should be more results under certain conditions.
#       - Major bugfix: reverse complement was not working in a correct way, so the primer 2 sequence was incorrect.
#       - For each amplicon a .fasta file will be generated (compatible with BiQ)
#       - remove "species" as required field in input.file.
#       - redesign of some output files.
#       - introduce mode ="exact" or "fast" for full range primersearch or quick and dirty runs.
#       - bugfix: hairpin construction had an error, now fixed.
#       - bugfix: hairpin analysis reports also fragments without linke in it...not yet fixed...
#       - new feature: ucsc compatible track for all designed primers
#v7.1:  - more new tracks: .bed file tracks and one more ucsc track using only the toplist_primers
#       - the other new track files can be used as direct input for BiSearch tool. 
#       - PLease be aware that BiSearch only features one assembly per species,
#       - therefore exact positions might differ from the primer tool.
#       - more REs and linker sequences for hairpin analyses.
#v7.2   - debugging.
#v7.2standalone - all helper functions are attached.
#               - added option to install missing Rpackages.
#v7.3standalone - fixed bug in instalation of missing packages.
#               - NCBI2R is not available by CRAN anymore. User has to manually install it from archive.
#v7.4standalone:  - debugging
#v7.5standalone:  - bugfix in hairpinizer module
#v7.6standalone:  - added DdeI Restriction enzyme to hairpin analysis.
#v8.0standalone:  - added SNP analysis (beta).
#v8.1standalone:  - added SNP analysis (beta), added SNP.ids. improved logging file.
#v8.2standalone:  - added filter options based on SNP analysis (beta).
#v8.3standalone:  - optional graphical outputs.
#v8.4standalone:  - disabled installation of missing packs.
#                 - maximum number of input regions [hardcoded].
#                 - maximum size (bp) of input regions [hardcoded].
#v8.5standalone:  - adjusted optional graphical outputs for hairpin jobs.
#v9.0standalone:  - include genomic primer design (beta).
#                 - lots of bugfixes.
#v9.2standalone:  - include repeat analysis (beta).
#v9.3standalone:  - include gene annotation analysis (beta).
#                - include CpG Island annotation analysis (beta).
#                - bugfixes (NOME graphical output no GpCs).
#v9.4standalone:  - bugfix - Exported .FASTA reference sequences are exported in genomic version (not bisulfite).
#v9.5standalone:  - bugfix if nor results for cpgi annotation
#v9.6standalone:  - bugfix: removed error if no primers are left after filtering for primer dimers.
#v9.7standalone:  - bugfix: graphics issue when some info tracks are longer than input regions.
#v9.8standalone:  - bugfix: some primerpairs are reported twice...fixed.
#v10.0standalone:  - added CLEVER primer design.
#v10.1standalone:  - checks input file if 'chr' column is in the correct format.
#v10.2standalone:  - fixed error when fetching annotaion for hg38 from ucsc
#v10.3standalone:  - fixed error when fetching annotaion for hg38 from ucsc,
#                  - added short info/explanation/comments for each argument of the main function.
###################################################################################################
#
# to do: - GUI
#        - local fetching of sequences for popular genomes (human/mouse, hg19,hg38,mm9,mm10)
#        - hairpin primers have issues: amplicons with LINKER at the end are also reported.
#
###################################################################################################
###################################################################################################

primer.design.pipeline<-function(table.in,#filename.in = NULL, # direct path to the file with i) regions (.bed format) or ii) sequences (filepath) (character)
                                 input.type="regions",#"regions" (for a file with regions in .bed format) or "sequences" (for a ".txt" file with sequences) (fixed vocabulary) (character)
                                 path.out = NULL, # directory where results will be stored (directory) (character)
                                 analysis.id="PrimerAnalysis",# Do not use'SonderZeichen'
                                 primer.type="bisulfite",#"NOME" or "bisulfite" or "hp_bisulfite" or "hp_NOME" or "genomic" or "hp_genomic" or "CLEVER" or "hp_CLEVER" (fixed vocabulary) (character) or "CrispRCas9PCR"
                                 mode="fast", #"exact" or "fast" (fixed vocabulary) (character)
                                 strand="top",#"both","top","bottom" (fixed vocabulary) (character)
                                 min.length.primer=23, #minimum length of primers (numeric)
                                 max.length.primer=34, #maximum length of primers (numeric)
                                 min.Tm.primer=48,     #minimum Tm of primers (numeric)
                                 max.Tm.primer=60,     #maximum Tm of primers (numeric)
                                 max.Tm.difference.primer=2, #maximum Tm difference of primerpairs (numeric)
                                 low.complexity.primer.removal=TRUE, #Removes low complex primers e.g. TATATATATATATA (TRUE/FALSE)
                                 max.bins.low.complexity=7,#max length of monomers (e.g. T stretches) (numeric)
                                 remove.primers.with.n=TRUE, #remove primers that contain 'N'bases
                                 min.length.amplicon=150, #minimum length of amplicon (numeric)
                                 max.length.amplicon=500, #maximum length of amplicon (numeric)
                                 min.number.gc.amplicon=0, #minimum number of GpCs in amplicon (set to 0 to disable) (numeric)
                                 min.number.cg.amplicon=5, #minimum number of CpGs in amplicon (set to 0 to disable) (numeric)
                                 use.target.columns=FALSE, #only select amplicons covering a preselected window in specified regions (beta) (TRUE/FALSE)
                                 primer.align.binsize=12,  #maximum size of allowed primer:primer self-alignments (numeric)
                                 min.C2T.primer1=3,   # minimum n of C to T conversions in primer1 (numeric)
                                 min.G2A.primer2=3,   # minimum n of G to A conversions in primer2  (numeric)
                                 check4snps=TRUE,     # check & annotate primers/amplicons for underlying SNPs (TRUE/FALSE)
                                 min.MAF.snp=0.01,  #max allowed minor allele frequency for SNPs
                                 min.snps.amplicon=0, #minimum n of allowed SNPs undelying amplicons (set to 0 to disable) (numeric)
                                 max.snps.amplicon=20,#maximum n of allowed SNPs undelying amplicons (set to 0 to disable) (numeric)
                                 min.snps.primer1=0, #minimum n of allowed SNPs undelying primer1 (set to 0 to disable) (numeric)
                                 max.snps.primer1=0, #maximum n of allowed SNPs undelying primer1 (set to 0 to disable) (numeric)
                                 min.snps.primer2=0, #minimum n of allowed SNPs undelying primer2 (set to 0 to disable) (numeric)
                                 max.snps.primer2=0, #maximum n of allowed SNPs undelying primer2 (set to 0 to disable) (numeric)
                                 check4repeats=FALSE,# check & annotate primers/amplicons for underlying repeats [uses Repeatmasker annotations] (TRUE/FALSE)
                                 allow.repeats.in.primers=FALSE, # allow repeats in primers (TRUE/FALSE)
                                 allow.repeats.in.amplicon=TRUE, # allow repeats in amplicons (TRUE/FALSE)
                                 annotate.genes=TRUE,# check & annotate primers/amplicons for underlying genes [uses RefSeq gene annotations] (TRUE/FALSE)
                                 annotate.cpg.islands=TRUE,# check & annotate primers/amplicons for underlying CpG islands [uses UCSC CGI annotations] (TRUE/FALSE)
                                 create.toplist=TRUE, # For each job with multiple primerdesigns a top performing primerpair is suggested (TRUE/FALSE)
                                 hp.filename=NA,              # required for hairpin analyses only (filename) (optional)
                                 hp.length.min=50,            # required for hairpin analyses only (minimum length of one arm in the hairpin molecule) (numeric)
                                 hp.length.max=300,           # required for hairpin analyses only (maximum length of one arm in the hairpin molecule) (numeric)
                                 #parallel=TRUE,               # CPU usage: TRUE= use multiple cores, FALSE= one core usage
                                 #n.unused.cpus=1,
                                 chop.size=30,#genomic primers only: input sequence slicing (numeric)
                                 add.ngs.adaptor.f="TCTTTCCCTACACGACGCTCTTCCGATCT", #adds this sequence 5' to primer1 to make it NGS compatible (optional) (character)
                                 add.ngs.adaptor.r="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT", #adds this sequence 5' to primer2 to make it NGS compatible (optional) (character)
                                 #install.missing.packages=FALSE,
                                 create.graphics=TRUE, #will generate a graphical overview for each job with succesfull primerdesign(s) (TRUE/FALSE)
                                 ...){          
 
  if(is.null(table.in) | is.null(path.out)){
  
    stop("input filename and export directory are missing.")  
    
  } 
  

scriptID<-"10.3standalone"
mydate<-paste(Sys.time())
start.wd<-getwd()
start.run<-mydate

maximum.input.regions<-100              # tool will not run if job has more regions
maximum.input.region.length.bp<-4000    # tool will not run if some regions are larger.

#snp.species<-c("Hsapiens","Mmusculus")

target.start.column.id="target.start"# specifiy start/end to force a 
target.end.column.id="target.end"# region of input seq to be in PCR product.
allow.primer.in.target.region=FALSE

#################################################################################
create.amplicon.barplots=FALSE
#################################################################################
#parallelization
#
#require(parallel)
#if(parallel){
#  n.cpus<-detectCores()
#  n.use.cpus<-n.cpus-n.unused.cpus
#  if(n.use.cpus<1){n.use.cpus<-1}
#  cl<-makeCluster(n.use.cpus, type="PSOCK")
#  }
#if(!parallel){
#  n.use.cpus<-1
#  cl<-makeCluster(n.use.cpus, type="PSOCK")
#  }
#################################################################################

external.packages<-c("R2HTML", "gsubfn", "rtf", "xtable", #"NCBI2R", 
                     "rtracklayer", "ggplot2")

#if(install.missing.packages){
#  missing.packs<-external.packages[!external.packages %in% rownames(installed.packages())]
#  
#  for (imp in 1:length(missing.packs)){
#    #temp.pack<-packs[imp]
#    try(install.packages(missing.packs[imp],quiet = TRUE),silent=TRUE)
#  }
#  
#  missing.packs<-external.packages[!external.packages %in% rownames(installed.packages())]
#  
#  if(length(missing.packs)>0){
#    stop(paste("Some packages could not be installed: ",missing.packs,sep=""))
#  }
#}

#################################################################################
require(R2HTML)

#################################################################################

#Set up dir
dir.create(path.out)
path.wd<-paste(path.out,analysis.id,"/",sep="/")
dir.create(path.wd)
setwd(path.wd)
#path.html<-paste(path.out,"html_",analysis.id,"/",sep="")
#print(paste("Working diectory set to: ",path.html,sep=""))
path.temp<-paste(path.wd,"temp/",sep="")
dir.create(path.temp)
path.warnings<-paste(path.wd,"warnings/",sep="")
dir.create(path.warnings)
path.sequences<-paste(path.wd,"sequences/",sep="")
dir.create(path.sequences)
path.tracks<-paste(path.wd,"tracks/",sep="")
dir.create(path.tracks)
if(create.graphics){
path.graphics<-paste(path.wd,"graphics/",sep="")
dir.create(path.graphics)}
#################################################################################
#################################################################################
#################################################################################
#set up logfile function
fn.log<-paste(path.wd,analysis.id,"#logfile.txt",sep="")
log<-function(text,add=TRUE){write.log(paste(text,sep=""),file=fn.log,add=add)}

log("Logfile initiated.",add=FALSE)
#################################################################################
#################################################################################
#################################################################################
if(min.length.primer<0){
  min.length.primer=13
  log(paste("auto re-assigned min.length.primer to: ",min.length.primer,sep=""))
}

if(max.length.primer<0){
  max.length.primer=min.length.primer+1
  log(paste("auto re-assigned max.length.primer to: ",max.length.primer,sep=""))
}

if(min.length.amplicon<0){
  min.length.amplicon=100
  log(paste("auto re-assigned min.length.amplicon to: ",min.length.amplicon,sep=""))
}

if(max.length.amplicon<0){
  max.length.amplicon=min.length.amplicon+1
  log(paste("auto re-assigned max.length.amplicon to: ",max.length.ampliconsep=""))
}

if(min.Tm.primer<0){
  min.Tm.primer=40
  log(paste("auto re-assigned min.Tm.primer to: ",min.Tm.primer,sep=""))
}

if(max.Tm.primer<0){
  max.Tm.primer=min.Tm.primer+1
  log(paste("auto re-assigned max.Tm.primer to: ",max.Tm.primer,sep=""))
}
if(min.number.gc.amplicon<0){
  min.number.gc.amplicon=0
  log(paste("auto re-assigned min.number.gc.amplicon to: ",min.number.gc.amplicon,sep=""))
} 

if(min.number.cg.amplicon<0){
  min.number.cg.amplicon=0
  log(paste("auto re-assigned min.number.cg.amplicon to: ",min.number.cg.amplicon,sep=""))
}

if(max.Tm.difference.primer<0){
  max.Tm.difference.primer=1
  log(paste("auto re-assigned max.Tm.difference.primer to: ",max.Tm.difference.primer,sep=""))
}

if(max.bins.low.complexity<0){
  max.bins.low.complexity=7
  log(paste("auto re-assigned max.bins.low.complexity to: ",max.bins.low.complexity,sep=""))
}

if(primer.align.binsize<0){
  primer.align.binsize=12
  log(paste("auto re-assigned primer.align.binsize to: ",primer.align.binsize,sep=""))
}

if(min.C2T.primer1<0){
  min.C2T.primer1=3
  log(paste("auto re-assigned min.C2T.primer1 to: ",min.C2T.primer1,sep=""))
}

if(min.G2A.primer2<0){
  min.G2A.primer2=3
  log(paste("auto re-assigned min.G2A.primer2 to: ",min.G2A.primer2,sep=""))
}

if(min.snps.amplicon<0){
  min.snps.amplicon=0
  log(paste("auto re-assigned min.snps.amplicon2 to: ",min.snps.amplicon,sep=""))
}

if(max.snps.amplicon<0){
  max.snps.amplicon=0
  log(paste("auto re-assigned max.snps.amplicon to: ",max.snps.amplicon,sep=""))
}

if(min.snps.primer1<0){
  min.snps.primer1=0
  log(paste("auto re-assigned min.snps.primer1 to: ",min.snps.primer1,sep=""))
}

if(max.snps.primer1<0){
  max.snps.primer1=0
  log(paste("auto re-assigned max.snps.primer1 to: ",max.snps.primer1,sep=""))
}

if(min.snps.primer2<0){
  min.snps.primer2=0
  log(paste("auto re-assigned min.snps.primer2 to: ",min.snps.primer2,sep=""))
}

if(max.snps.primer2<0){
  max.snps.primer2=0
  log(paste("auto re-assigned max.snps.primer2 to: ",max.snps.primer2,sep=""))
}

if(hp.length.min<0){
  hp.length.min=50
  log(paste("auto re-assigned hp.length.min to: ",hp.length.min,sep=""))
}

if(hp.length.max<0){
  hp.length.max=0
  log(paste("auto re-assigned hp.length.max to: ",hp.length.max,sep=""))
}

if(max.length.primer < min.length.primer){
  max.length.primer <- min.length.primer+1
    log(paste("auto re-assigned max.length.primer to: ",max.length.primer,sep=""))
}

if(min.length.primer > max.length.primer){
  min.length.primer <- max.length.primer-1
  log(paste("auto re-assigned min.length.primer to: ",min.length.primer,sep=""))
}

if(max.length.amplicon < min.length.amplicon){
  max.length.amplicon <- min.length.amplicon+1
  log(paste("auto re-assigned max.length.amplicon to: ",max.length.amplicon,sep=""))
}

if(min.length.amplicon > max.length.amplicon){
  min.length.amplicon <- max.length.amplicon-1
  log(paste("auto re-assigned min.length.amplicon to: ",min.length.amplicon,sep=""))
}

if(max.Tm.primer < min.Tm.primer){
  max.Tm.primer <- min.Tm.primer+1
  log(paste("auto re-assigned max.Tm.primer to: ",max.Tm.primer,sep=""))
}

if(min.Tm.primer > max.Tm.primer){
  min.Tm.primer <- max.Tm.primer-1
  log(paste("auto re-assigned min.Tm.primer to: ",min.Tm.primer,sep=""))
}

if(min.snps.amplicon > max.snps.amplicon){
  min.snps.amplicon <- max.snps.amplicon-1
  log(paste("auto re-assigned min.snps.amplicon to: ",min.snps.amplicon,sep=""))
}

if(max.snps.amplicon < min.snps.amplicon){
  max.snps.amplicon <- min.snps.amplicon+1
  log(paste("auto re-assigned max.snps.amplicon to: ",max.snps.amplicon,sep=""))
}

if(min.snps.primer1 > max.snps.primer1){
  min.snps.primer1 <- max.snps.primer1-1
  log(paste("auto re-assigned min.snps.primer1 to: ",min.snps.primer1,sep=""))
}

if(max.snps.primer1 < min.snps.primer1){
  max.snps.primer1 <- min.snps.primer1+1
  log(paste("auto re-assigned max.snps.primer1 to: ",max.snps.primer1,sep=""))
}

if(min.snps.primer2 > max.snps.primer2){
  min.snps.primer2 <- max.snps.primer2-1
  log(paste("auto re-assigned min.snps.primer2 to: ",min.snps.primer2,sep=""))
}

if(max.snps.primer2 < min.snps.primer2){
  max.snps.primer2 <- min.snps.primer2+1
  log(paste("auto re-assigned max.snps.primer2 to: ",max.snps.primer2,sep=""))
}

if(hp.length.min > hp.length.max){
  hp.length.min <- hp.length.max-1
  log(paste("auto re-assigned hp.length.min to: ",hp.length.min,sep=""))
}

if(hp.length.max < hp.length.min){
  hp.length.max <- mhp.length.min+1
  log(paste("auto re-assigned hp.length.max to: ",hp.length.max,sep=""))
}

######################################Setting the input type##########################################

if(c("chr") %in% colnames(table.in)){
  input.type = "regions"
  log(paste("identified input.type as regions!"))
}

if(c("sequence") %in% colnames(table.in)){
  input.type = "sequences"
  log(paste("identified input.type as sequences!"))
}

#################################################################################
#################################################################################
#################################################################################
if(primer.type=="genomic"){
  if (chop.size < min.length.primer){
    chop.size<-min.length.primer
    log(paste("auto re-assigned chop.size to: ",chop.size,sep=""))
    
  }
  
  if (chop.size < max.length.primer){
    chop.size<-max.length.primer
    log(paste("auto re-assigned chop.size to: ",chop.size,sep=""))
    
  }
}#if(primer.type=="genomic"){
#################################################################################
#################################################################################
#################################################################################
# html.report<-function(filename.out="html.report.html",filenames,txt.header=TRUE,txt.sep="\t"){
  
  # require(R2HTML)
  # require(tools)
  
  # htmlWidth<-600
  # htmlHeight<-600
  # GraphBorderSize<-5
  
  # graphics<-c("png","jpeg","svg","jpg")
  # txt<-c("txt","csv","tsv","fasta","bed")
  # graphics<-c(graphics,toupper(graphics))
  # txt<-c(txt,toupper(txt))
  
  # for (i in 1:length(filenames)){
    # temp<-file_path_as_absolute(filenames[i])
    # splits<-strsplit(temp,"/")
    # fn<-splits[[1]][length(splits[[1]])]
    # file.type<-file_ext(fn)
    # temp.name<-gsub(file.type,"",fn)
    
    # if(file.type %in% graphics){
      # HTMLInsertGraph(Caption=paste(temp.name),
                      # file=filename.out,append=TRUE,
                      # GraphFileName=temp,
                      # GraphBorder=GraphBorderSize,
                      # Width=htmlWidth,Height=htmlHeight)
      # HTMLhr(file=filename.out,append=TRUE)
    # }#if graphix
    
    # if(file.type %in% txt){
      # df<-read.table(filenames[i],header=txt.header,sep=txt.sep)
      # if(nrow(df)==0 | ncol(df)<2){
        # df<-readLines(filenames[i])
      # }
      # HTML(x=df,file=filename.out,
           # classcaption="captiondataframe",
           # caption=paste(temp.name),
           # Border=1,innerBorder=1,append=TRUE,align="center",
           # captionalign="top",row.names=FALSE)
      # HTMLhr(file=filename.out,append=TRUE)
    # }#if txt
    
  # }#for i
# }
#################################################################################
#################################################################################
#################################################################################
#################################################################################
##### use rtracklayer to get data (e.g. snps) from ucsc:

ucsc.info<-function(assembly="hg19",
                    chr=NULL,                   # c("chr1","chr4","chr5"...)
                    start=NULL,                #c(1234567,1234567,1234567)
                    end=NULL,                  #c(1235567,1235567,1235567)
                    strand=NULL,               #c(+","-","+")
                    seqinfo=NULL,
                    track="snp147Common",###'rmsk' for repeats
                    table="snp147Common",###'rmsk' for repeats
                    browser="UCSC"){
  
  if(is.null(chr) | is.null(start) | is.null(end)){
    stop("input of chromosome, start and end are mandatory")
  }
  
  if(length(track)>1 | length(table)>1 | length(browser)>1){
    stop("only 1 entry for track, table and browser")
  }
  
  if(length(assembly)>1){
    warning("only one assembly per session allowed...use only first one")
    
    assembly.chosen<-assembly[1]
    
    filter.by.assembly<-assembly==assembly.chosen
    assembly<-assembly[filter.by.assembly]
    chr<-chr[filter.by.assembly]
    start<-start[filter.by.assembly]
    end<-end[filter.by.assembly]
    
    if(!is.null(strand)){
      strand<-strand[filter.by.assembly]
    }
    
    if(!is.null(seqinfo)){
      seqinfo<-seqinfo[filter.by.assembly]
    }
  }
  
  user <- unname(Sys.info()["user"])
  
# if (user == "shiny") {
# 
#   # Set library locations
#   ds <- .libPaths(c("/home/users/amerg/R/x86_64-redhat-linux-gnu-library/3.3","/opt/Rlib/3.2","/usr/lib64/R/library","/usr/share/R/library"))
# print(ds)
# }
  require (rtracklayer)
  mySession = browserSession(browser)
  genome(mySession) <- assembly
  
  #GRanges(seqnames=NULL,ranges=NULL, strand=NULL,
  #..., seqlengths=NULL, seqinfo=NULL): Creates a GRanges object.
  #
  #see help GRanges
  #gr0 <- GRanges(seqnames = chr, IRanges(start = start,end = end))
#  
#  if(!is.null(strand)){strand(gr0)<-strand}
#  
#  if(!is.null(seqinfo)){
#    to.seq.info <- Seqinfo(seqinfo)
#    #seqinfo(gr0) <- merge(seqinfo(gr0), to.seq.info)
#    seqinfo(gr0)<-merge(seqinfo(gr0), to.seq.info)
#    #seqlevels(gr0) <- seqlevels(to.seq.info)
#  }
  
  myquer<-ucscTableQuery(mySession,track,GRangesForUCSCGenome(assembly,chr,IRanges(start,end,names=seqinfo)))
  
  tableName(myquer)<- paste(table)

  #genebrowser.info <- getTable(session = mySession, trackName = track, tableName = table,)
  genebrowser.info <- getTable(myquer)
  
  return(genebrowser.info)
  
}# end of ucsc function

##################################################################################################
##################################################################################################
##################################################################################################
# #create main .html file
# log("Create HTML scaffolds...")

# if(primer.type=="hp_bisulfite" | primer.type=="hp_NOME" | primer.type=="hp_genomic" | primer.type=="hp_CLEVER"){
  # links<-c(paste("html_",analysis.id,"/settings_",analysis.id,".html",sep=""),
          # paste("html_",analysis.id,"/bedfile_",analysis.id,".html",sep=""),
          # paste("html_",analysis.id,"/RestrictionEnzymesAndLinkerInfo_",analysis.id,".html",sep=""),
          # paste("html_",analysis.id,"/Hairpinizer_Results_",analysis.id,".html",sep=""),
          # paste("html_",analysis.id,"/summary_",analysis.id,".html",sep=""),
          # paste("html_",analysis.id,"/primer_",analysis.id,".html",sep=""),
          # paste("html_",analysis.id,"/primer_",analysis.id,"_toplist.html",sep=""))

  # labs<-c("Settings","InputFile","Restriction Enzymes & Linkers","Hairpinizer Results","Summary","Primer","Primer_Toplist")
# }
# if(primer.type=="bisulfite" | primer.type=="NOME" | primer.type=="genomic" | primer.type=="CLEVER" |primer.type=="CrispRCas9PCR"){
  # links<-c(paste("html_",analysis.id,"/settings_",analysis.id,".html",sep=""),
           # paste("html_",analysis.id,"/bedfile_",analysis.id,".html",sep=""),
           # paste("html_",analysis.id,"/summary_",analysis.id,".html",sep=""),
           # paste("html_",analysis.id,"/primer_",analysis.id,".html",sep=""),
           # paste("html_",analysis.id,"/primer_",analysis.id,"_toplist.html",sep=""))
  
  # labs<-c("Settings","InputFile","Summary","Primer","Primer_Toplist")
# }
# html.main.page(fn.links = links,labels = labs,filename = paste(path.out,"PrimerDesign_",analysis.id,".html",sep=""))

# ###############################################################################

# #create overview .html file
settings<-data.frame(PARAMETERS=c("Date",
                                  "Script_Version",
                                  #"Name of Inputfile",
                                  "Type of Inputfile",
                                  "Output Path",
                                  "Analysis ID",
                                  "Primer Design Type",
                                  "Primer Design Mode",
                                  "Analyzed DNA Strands",
                                  "Minimum Length of Primers",
                                  "Maximum Length of Primers",
                                  "Minimum Tm of Primers",
                                  "Maximum Tm of Primers",
                                  "Maximum Tm difference",
                                  "Remove Low Complexity Primers",
                                  "Maximum Number of Low Complex Bins (1mers/2mers/3mers)",
                                  "Remove 'N'containing primers",
                                  "Minimum Length of Amplicon",
                                  "Maximum Length of Amplicon",
                                  "Use target columns",
                                  "Minimum Number of GpCs in Amplicon",
                                  "Minimum Number of CpGs in Amplicon",
                                  paste("Minimum 'C' to 'T' conversions in primer 1",sep=""),
                                  paste("Minimum 'G' to 'A' conversions in primer 2",sep=""),
                                  "Primer Dimer BinSize",
                                  "Check for Primer SNPs",
                                  "Minimum MAF of SNPs",
                                  "Minimum SNPs in amplicon",
                                  "Maximum SNPs in amplicon",
                                  "Minimum SNPs in primer1",
                                  "Maximum SNPs in primer1",
                                  "Minimum SNPs in primer2",
                                  "Maximum SNPs in primer2",
                                  "Check for Repeats",
                                  "Allow Repeats in Primers",
                                  "Allow Repeats in Amplicon",
                                  "Annotate Genes",
                                  "Annotate CpG Islands",
                                  "Create Toplist Table",
                                  "Create graphics"),
                SETTINGS=c(mydate,
                           paste(scriptID),
                           #paste(filename.in),
                           paste(input.type),
                           paste(path.out),
                           paste(analysis.id),
                           paste(primer.type),
                           paste(mode),
                           paste(strand),
                           paste(min.length.primer),
                           paste(max.length.primer),
                           paste(min.Tm.primer),
                           paste(max.Tm.primer),
                           paste(max.Tm.difference.primer),
                           paste(low.complexity.primer.removal),
                           paste(max.bins.low.complexity),
                           paste(remove.primers.with.n),
                           paste(min.length.amplicon),
                           paste(max.length.amplicon),
                           paste(use.target.columns),
                           paste(min.number.gc.amplicon),
                           paste(min.number.cg.amplicon),
                           paste(min.C2T.primer1),
                           paste(min.G2A.primer2),
                           paste(primer.align.binsize),
                           paste(check4snps),
                           paste(min.MAF.snp),
                           paste(min.snps.amplicon),
                           paste(max.snps.amplicon),
                           paste(min.snps.primer1),
                           paste(max.snps.primer1),
                           paste(min.snps.primer2),
                           paste(max.snps.primer2),
                           paste(check4repeats),
                           paste(allow.repeats.in.primers),
                           paste(allow.repeats.in.amplicon),
                           paste(annotate.genes),
                           paste(annotate.cpg.islands),
                           paste(create.toplist),
                           paste(create.graphics)))

if(primer.type=="hp_bisulfite" | primer.type=="hp_NOME" | primer.type=="hp_genomic" | primer.type=="hp_CLEVER"){
  
  settings<-rbind(settings,
                  data.frame(PARAMETERS=c("hp.filename",
                                          "hp.length.min",
                                          "hp.length.max"),
                             SETTINGS=c(hp.filename,
                                        hp.length.min,
                                        hp.length.max)))
}

if(primer.type=="hp_genomic" | primer.type=="genomic"){
  
  settings<-rbind(settings,
                  data.frame(PARAMETERS="input.chop.size",
                             SETTINGS=chop.size))
}

 write.table(settings,file=paste(path.wd,"settings_",analysis.id,".txt",sep=""),
             col.names=T,row.names=F,sep="\t",dec=".",quote=F)
# # html.report(filenames =  paste(path.html,"settings_",analysis.id,".txt",sep=""),
            # # filename.out = paste(path.html,"settings_",analysis.id,".html",sep=""),txt.header = TRUE,txt.sep = "\t")
log("Done.")
#################################################################################

#start data import
log("Module: Data Import.")

if (input.type=="regions"){
  log("Fetch sequences from UCSC...")
  
  #iregs <- filetable.in #iregs<-read.table(filename.in,header=T,sep="\t",dec=".")
  iregs<- table.in
  print(iregs)
  
  def.cols<-c("chr","start","end","assembly","sequenceID","sequence.length","sequence.adress","sequence")
  print(def.cols)
  nregs<-nrow(iregs)#number of input regions
  print(nregs)
  
  
  if("sequenceID" %in% colnames(iregs)){
    
    iregs = iregs[order(iregs$sequenceID),]
    reg.rle = rle(as.character(iregs$sequenceID))
    
    if(any(reg.rle$lengths > 1)){
      
      ids.new <- paste0(rep(reg.rle$values, times = reg.rle$lengths), "_",unlist(lapply(reg.rle$lengths, seq_len)))      
      iregs$sequenceID = ids.new 
      log("Found duplicated sequence IDs...unify IDs") 
      
    }
  } 
  
  if(nregs>maximum.input.regions){
    
    stop(paste("Too many input regions [max: ",maximum.input.regions,"]",sep=""))
    
  }
  
  if(!all(grepl("chr",iregs$chr))){
    
    stop(paste("'chr' column in input file is not in the correct format ('chr'[chr])",sep=""))
    
  }
  
  iregchr<-as.character(iregs$chr)
  print(iregchr)
  bed<-data.frame(chr=rep(NA,nregs),start=rep(NA,nregs),end=rep(NA,nregs),
                  assembly=rep(NA,nregs),
                  sequenceID=rep(NA,nregs),sequence.length=rep(NA,nregs),
                  sequence.adress=rep(NA,nregs),sequence=rep(NA,nregs))
				  print(bed)
  bed$chr<-as.character(bed$chr)
  bed$chr<-iregchr
  print(bed$chr)
  bed$start<-as.numeric(bed$start)
  print(bed$start)
  bed$end<-as.numeric(bed$end)
  print(bed$end)
  undef.cols<-colnames(iregs)[!colnames(iregs) %in% def.cols]
  print(undef.cols)
  bed[,undef.cols]<-NA
  print(bed[,undef.cols])
  #for each entry in input.file retrieve dna seq
  for (ir in 1:nregs){
    dna<-getDNA(iregs[ir,"chr"],iregs[ir,"start"],iregs[ir,"end"],iregs[ir,"assembly"])

    sid<-iregs[ir,"sequenceID"]
    if(is.na(sid)|sid==""){
      sid<-paste(c(iregs[ir,"assembly"],"_",iregs[ir,"chr"],"_",iregs[ir,"start"],"_",iregs[ir,"end"]),sep="")
    }
    adr<-paste(c(as.character(iregs[ir,"chr"]),"_",iregs[ir,"start"],"_",iregs[ir,"end"]),sep="")
    bed[ir,c("start",
             "end",
             "assembly",
             "sequenceID",
             "sequence.length",
             "sequence.adress",
             "sequence")]<-c(iregs[ir,"start"],iregs[ir,"end"],paste(iregs[ir,"assembly"]),
                             paste(iregs[ir,"species"]),paste(sid),nchar(dna),paste(adr,collapse=""),
                                                                                                       dna)   
    bed[ir,undef.cols]<-iregs[ir,undef.cols]
  }
   print(bed)
  if(any(nchar(bed$sequence)) > maximum.input.region.length.bp){
    stop(paste(" one or more input regions are too large [max: ",maximum.input.region.length.bp,"]",sep=""))
  }
  
#######################################################################################

#remove entries with no sequence retrieved
  bad.seq.f<-bed$sequence.length<1  
  bad.out<-bed[bad.seq.f,]
  
  if(nrow(bad.out)>0){
    write.table(bad.out,file=paste(path.sequences,analysis.id,"EntriesWithNoSequences.txt",sep=""),sep="\t",dec=".",
                col.names=T,row.names=F)  
  }
  
  bed<-bed[!bad.seq.f,]

###################################################################################

#perform input sequence characterization
  log("Start characterization of input sequences...")
  
  seq.ana.all<-list()
  for (iiii in 1:nrow(bed)){
    
    seqanai<-sequence.characterization(paste(bed[iiii,"sequence"]))
    write.table(seqanai,file=paste(path.sequences,"sequence.characterization_",bed[iiii,"sequenceID"],".txt",sep=""),
                sep="\t",dec=".",col.names=T,row.names = FALSE,quote=F)
    
    seq.ana.all[[paste(bed[iiii,"sequenceID"])]]<-seqanai
  }

  log("Done.")
  
#######################################################################################
#######################################################################################
#######################################################################################

#retrieve snp information for input sequences  
  if (check4snps){
    
    supported.assemblies.snp<-c("hg18","hg19","hg38","mm9","mm10")
    
    log("Start SNP analysis...")
    
    if(input.type=="sequences"){
      #should check biomart....
      log(paste("WARNING: To check for SNPs input type must be regions."))
      fileConn<-file(paste(path.warnings,"Warning.SNPcheck.data.error.txt",sep=""))
      writeLines(c("WARNING: To check for SNPs input type must be regions."), fileConn)
      close(fileConn)
      hp.initial.input.type <- input.type
    }
    
    
    if(input.type=="regions"){
      fetched_snp = mapply(rtl.get.snp.info.by.region,
                           paste0(bed$assembly),
                           paste0(bed$chr),
                           as.numeric(as.character(bed$start)),
                           as.numeric(as.character(bed$end)))
        
      for(isnp in 1:ncol(fetched_snp)){
        if(isnp == 1){
            all_my_snps = as.data.frame(fetched_snp[,isnp])
        }
      
        if(isnp > 1){
          all_my_snps = rbind(all_my_snps,as.data.frame(fetched_snp[,isnp]))
        }
      }
      
      colnames(all_my_snps) = gsub("chrom","chr",
                                   gsub("chromStart","start",
                                        gsub("chromEnd","end",
                                             gsub("submitters","source",
                                                  gsub("name","rs_id",colnames(all_my_snps))))))

      write.table(all_my_snps,file=paste(path.tracks,"#SNP.info.input.regions.txt",sep=""),
                  sep="\t",dec=".",col.names=T,row.names=F,quote=F)
      
      
  }#if regions
}#check4snps

  
#######################################################################################
#######################################################################################
#######################################################################################

#retrieve repeat information for input sequences  
if (check4repeats){
  
  supported.assemblies.snp<-c("hg18","hg19","hg38","mm9","mm10")
  
  log("Start Repeat analysis...")
  
  if(input.type=="sequences"){
    
    #should check biomart....
    log(paste("WARNING: To check for repeats input type must be regions."))
    fileConn<-file(paste(path.warnings,"Warning.Repeatcheck.data.error.txt",sep=""))
    writeLines(c("WARNING: To check for repeats input type must be regions."), fileConn)
    close(fileConn)
  }
  
  
  if(input.type=="regions"){
    
    assems<-paste(unique(bed[,"assembly"]))
    
    for (iasem in 1:length(assems)){
      
      assem<-assems[iasem]
      bedia<-bed[as.character(bed$assembly)==assem,]
      
      if(!assem %in% supported.assemblies.snp){
        log(paste("...assembly [",assem, "] not supported for repeat analysis",sep=""))   
        next()
      }
      
      if(assem %in% supported.assemblies.snp){
        
        chrom<-paste(bedia[,"chr"])
        allstarts<-as.numeric(as.character(bedia[,"start"]))
        allends<-as.numeric(as.character(bedia[,"end"]))
        
        my_repeats<-fetch.repeat.info.rest(assembly = assem,
                           chr = chrom,
                           start = allstarts,
                           end = allends)
        
        if(nrow(my_repeats)>0){
        
        my_repeats$assembly<-assem
        
        }
        
        if(nrow(my_repeats)>0){
         
          m_repeats<-data.frame(assembly=assem,chr=chrom,start=NA,end=NA,track=NA,table=NA)
           
        }
        
        if(iasem==1){
          all_my_repeats<-my_repeats
        }#if assem
        
        if(iasem>1){
          all_my_repeats<-rbind(all_my_repeats,my_repeats)
        }#if >1
        
      }#if assem
    }#for assem
    write.table(all_my_repeats,file=paste(path.tracks,"#repeat.info.input.regions.txt",sep=""),
                sep="\t",dec=".",col.names=T,row.names=F,quote=F)
  }#if regions
}#check4repeats

#######################################################################################
#######################################################################################
#######################################################################################

#retrieve gene information for input sequences  
if (annotate.genes){
  
  supported.assemblies.snp<-c("hg18","hg19","hg38","mm9","mm10")
  
  log("Start Gene annotation analysis...")
  
  if(input.type=="sequences"){
    
    #should check biomart....
    log(paste("WARNING: To check for gene annotations input type must be regions."))
    fileConn<-file(paste(path.warnings,"Warning.GeneAnnotationcheck.data.error.txt",sep=""))
    writeLines(c("WARNING: To check for Gene Annotation input type must be regions."), fileConn)
    close(fileConn)
  }
  
  
  if(input.type=="regions"){
    
    assems<-paste(unique(bed[,"assembly"]))
    
    for (iasem in 1:length(assems)){
      
      assem<-assems[iasem]
      bedia<-bed[as.character(bed$assembly)==assem,]
      
      if(!assem %in% supported.assemblies.snp){
        log(paste("...assembly [",assem, "] not supported for gene annotation analysis",sep=""))   
        next()
      }
      
      if(assem %in% supported.assemblies.snp){
        
        chrom<-paste(bedia[,"chr"])
        allstarts<-as.numeric(as.character(bedia[,"start"]))
        allends<-as.numeric(as.character(bedia[,"end"]))
        
        my_genes<-fetch.gene.info.rest(assembly = assem,
                              chr = chrom,
                              start = allstarts,
                              end = allends)
        
        my_genes$assembly<-assem
        
        if(iasem==1){
          all_my_genes<-my_genes
        }#if assem
        
        if(iasem>1){
          all_my_genes<-rbind(all_my_genes,my_genes)
        }#if >1
        
      }#if assem
    }#for assem
    
    write.table(all_my_genes,file=paste(path.tracks,"#gene.info.input.regions.txt",sep=""),
                sep="\t",dec=".",col.names=T,row.names=F,quote=F)
  }#if regions
}#annotate.genes

#######################################################################################
#######################################################################################
#######################################################################################

#retrieve cpg island information for input sequences  
if (annotate.genes){
  
  supported.assemblies.snp<-c("hg18","hg19","hg38","mm9","mm10")
  
  log("Start CpG Island annotation analysis...")
  
  if(input.type=="sequences"){
    
    #should check biomart....
    log(paste("WARNING: To check for CpG island annotations input type must be regions."))
    fileConn<-file(paste(path.warnings,"Warning.CpGIcheck.data.error.txt",sep=""))
    writeLines(c("WARNING: To check for CpG Islands input type must be regions."), fileConn)
    close(fileConn)
  }
  
  
  if(input.type=="regions"){
    
    assems<-paste(unique(bed[,"assembly"]))
    
    for (iasem in 1:length(assems)){
      
      assem<-assems[iasem]
      bedia<-bed[as.character(bed$assembly)==assem,]
      
      if(!assem %in% supported.assemblies.snp){
        log(paste("...assembly [",assem, "] not supported for cpg island annotation analysis",sep=""))   
        next()
      }
      
      if(assem %in% supported.assemblies.snp){
        
        # cpgitrack<-switch(assem, 
        #                   "hg18" = "cpgIslandExt",
        #                   "hg19" = "cpgIslandExt",
        #                   "hg38" = "cpgIslandExt",
        #                   "mm9"  = "cpgIslandExt",
        #                   "mm10" = "cpgIslandExt")
        # 
        # cpgitable<-switch(assem, 
        #                   "hg18" = "cpgIslandExt",
        #                   "hg19" = "cpgIslandExt",
        #                   "hg38" = "cpgIslandExt",
        #                   "mm9"  = "cpgIslandExt",
        #                   "mm10" = "cpgIslandExt")
        
        chrom<-paste(bedia[,"chr"])
        allstarts<-as.numeric(as.character(bedia[,"start"]))
        allends<-as.numeric(as.character(bedia[,"end"]))
        
        my_cpgi<-fetch.gene.info.rest(assembly = assem,
                            chr = chrom,
                            start = allstarts,
                            end = allends)
        
        if(!is.null(nrow(my_cpgi)) && nrow(my_cpgi)>0){
        
        my_cpgi$assembly<-assem
        my_cpgi$cpgi.db<-cpgitable
        
        }
        
        if(iasem==1){
          all_my_cpgis<-my_cpgi
        }#if assem
        
        if(iasem>1){
          all_my_cpgis<-rbind(all_my_cpgis,my_cpgi)
        }#if >1
        
      }#if assem
    }#for assem
    
    write.table(all_my_cpgis,file=paste(path.tracks,"#cpgi.info.input.regions.txt",sep=""),
                sep="\t",dec=".",col.names=T,row.names=F,quote=F)
  }#if regions
}#annotate.cpgi
  
#######################################################################################
#######################################################################################
#######################################################################################
    
  if(primer.type == "bisulfite" | primer.type == "hp_bisulfite" | primer.type == "CrispRCas9PCR"){
   bed$sequence.converted<-sapply(bed$sequence,bisulfite.conversion,strand=strand)
   }
  
  if(primer.type == "NOME" | primer.type == "hp_NOME"){
    bed$sequence.converted<-sapply(bed$sequence,gcCon,strand=strand)   
  }
  
  if(primer.type == "genomic" | primer.type == "hp_genomic" | primer.type == "CLEVER" | primer.type == "hp_CLEVER" ){
    bed$sequence.converted<-bed$sequence
  }
  
  hp.initial.input.type<-input.type
  
  bed$Index<-1:nrow(bed)
  write.table(bed,file=paste(path.sequences,analysis.id,"_PrimerDesign_InputSequences.txt",sep=""),sep="\t",dec=".",
              col.names=TRUE,row.names=FALSE,quote=F)
  #html.report(filenames =  paste(path.sequences,analysis.id,"_PrimerDesign_InputSequences.txt",sep=""),
              #filename.out = paste(path.html,"bedfile_",analysis.id,".html",sep=""),txt.header = TRUE,txt.sep = "\t")
  log("Done.")
}

if (input.type=="sequences"){
  
  log("Import Sequences from .txt File...")
  bed<-table.in #read.table(file=filename.in,sep="\t",dec=".",header=TRUE)
  #html.report(filenames =  filename.in,filename.out = paste(path.html,"bedfile_",analysis.id,".html",sep=""),txt.header = TRUE,txt.sep = "\t")
  log("Done.")
  
  log("Start characterization of input sequences...")

  seq.ana.all<-list()

  for (iiii in 1:nrow(bed)){
    seqanai<-sequence.characterization(paste(bed[iiii,"sequence"]))
    write.table(seqanai,file=paste(path.sequences,"sequence.characterization_",bed[iiii,"sequenceID"],".txt",sep=""),
                sep="\t",dec=".",col.names=T,row.names = FALSE,quote=F)
    seq.ana.all[[paste(bed[iiii,"sequenceID"])]]<-seqanai
  }
  log("Done.")
  
  if("sequenceID" %in% colnames(bed)){
    
    bed = bed[order(bed$sequenceID),]
    reg.rle = rle(as.character(bed$sequenceID))
    
    if(any(reg.rle$lengths > 1)){
      
      ids.new <- paste0(rep(reg.rle$values, times = reg.rle$lengths), "_",unlist(lapply(reg.rle$lengths, seq_len)))      
      bed$sequenceID = ids.new 
      
      log("Found duplicated sequence IDs...unify IDs") 
    }
  } 

  }
################################################################################

##hairpin specific section: calculate restriction digests
if(primer.type=="hp_bisulfite" | primer.type=="hp_NOME" | primer.type=="hp_genomic" | primer.type=="hp_CLEVER"){
  
  if(!is.na(hp.filename)){
    hp.info<-read.table(file=paste(hp.filename),header=TRUE,sep="\t",dec=".")
    write.table(hp.info,file=paste(path.sequences,"RestrictionEnzymesAndLinkerInfo_",analysis.id,".txt",sep=""),
                sep="\t",dec=".",col.names=TRUE,row.names=FALSE)     
  }
  
  if(is.na(hp.filename)){
    log("Required hairpin infofile is missing...use default hairpin settings...")
    hp.info<-data.frame(RE.id=c("MspI",
                                "AluI",
                                "TaqI",
                                "Nla3",
                                "PstI",
                                "BamHI",
                                "DdeI",
                                "BsaWI",
                                "Eco47I",
                                "PpuMI",
                                "BlpI",
                                "Bpu10I",
                                "Bsu36I"),
                        RE.seq=c("CCGG",
                                 "AGCT",
                                 "TCGA",
                                 "CATG",
                                 "CTGCAG",
                                 "GGATCC",
                                 "CTNAG",
                                 "WCCGGW",
                                 "GGWCC",
                                 "RGGWCCY",
                                 "GCTNAGC",
                                 "CCTNAGC",
                                 "CCTNAGG"),
                        linker.sequence=c(toupper("CGGGGCCCATddddddddATGGGCCC"),
                                          toupper("GGGCCCATddddddddATGGGCCC"),
                                          toupper("CGGGGCCCATddddddddATGGGCCC"),
                                          toupper("GGGCCTAATATAGTATAGGCCCCATG"),
                                          toupper("GGGCCTAATATAGTATAGGCCCTGCA"),
                                          toupper("GATCGGGCCCATddddddddATGGGCCC"),
                                          toupper("tnagggSccatddddddddatgggScc"),
                                          toupper("CCGGCGG5G6GATDDDDDDDDATCG7G8CG"),
                                          toupper("GWCGGG5G6GATDDDDDDDDATCG7G8CC"),
                                          toupper("GWCGGG5G6GATDDDDDDDDATCG7G8CC"),
                                          toupper("TNAGCG5G6GATDDDDDDDDATCG7G8GC"),
                                          toupper("TNAGCG5G6GATDDDDDDDDATCG7G8GC"),
                                          toupper("TNAGCG5G6GATDDDDDDDDATCG7G8GC")))
    write.table(hp.info,file=paste(path.sequences,"RestrictionEnzymesAndLinkerInfo_",analysis.id,".txt",sep=""),sep="\t",dec=".",col.names=TRUE,row.names=FALSE)
  }

################################################################################

#add webpage to hmtl report
    # html.report(filenames = paste(path.sequences,"RestrictionEnzymesAndLinkerInfo_",analysis.id,".txt",sep=""),
                # filename.out = paste(path.html,"RestrictionEnzymesAndLinkerInfo_",analysis.id,".html",sep=""),
                # txt.header = TRUE,txt.sep = "\t")
    
################################################################################

    hp.res.ids<-paste(hp.info[,"RE.id"])
    hp.res.seqs<-paste(hp.info[,"RE.seq"])
    hp.linker.seqs<-paste(hp.info[,"linker.sequence"])
  
###################################################################################### 
  
#Module: Hairpinizer
    log("Module: hairpinizer...")
    
    hp.res<-list()
    for (iib in 1:nrow(bed)){
      ihp.id<-paste(bed[iib,"sequenceID"])
      ihp.seq<-paste(bed[iib,"sequence"])
      log(paste("Module: hairpinizer...analyze sequence ",ihp.id,sep=""))
      
      for (ihp in 1:nrow(hp.info)){
        ihp.res.id<-hp.res.ids[ihp]
        ihp.res.seq<-hp.res.seqs[ihp]
        ihp.linker.seq<-hp.linker.seqs[ihp]        
        
        log(paste("Analyze RE ",ihp.res.id,"...",sep=""))
        hp.res[[paste(ihp.id,"_",ihp.res.id,sep="")]]<-hairpinizer(ihp.seq,
                                                                   inputID=ihp.id,
                                                                   REseq=ihp.res.seq,
                                                                   RE.id=ihp.res.id,
                                                                   linker.seq=ihp.linker.seq,
                                                                   length.hairpin.region=hp.length.max)
        log(paste("Finished RE ",ihp.res.id,sep=""))
      }
      log(paste("Module: hairpinizer...finished sequence ",ihp.id,sep=""))
    } 
    log("Finished hairpinizer.")

    log("Export results...")
    
    for (iib in 1:nrow(bed)){
      ihp.id<-paste(bed[iib,"sequenceID"])
      hps<-do.call(what = rbind,hp.res[grepl(ihp.id,names(hp.res))])
      hps$id<-rownames(hps)
      write.table(hps,file=paste(path.wd,analysis.id,"_Hairpinizer_ResultsBySequence_",ihp.id,".txt",sep=""),
                  sep="\t",dec=".",col.names=TRUE,row.names=FALSE,quote=F)
    }

    for (ihp in 1:nrow(hp.info)){
      ihp.res.id<-hp.res.ids[ihp]
      hpr<-do.call(what = rbind,hp.res[grepl(ihp.res.id,names(hp.res))])
      hpr$id<-rownames(hpr)#[grepl(ihp.res.id,names(hp.res))])
      write.table(hpr,file=paste(path.wd,analysis.id,"_Hairpinizer_ResultsByEnzyme_",ihp.res.id,".txt",sep=""),
                  sep="\t",dec=".",col.names=TRUE,row.names=FALSE,quote=F)
    }
    
    log("Done.")
################################################################################

#complete results
    hp<-do.call(rbind,hp.res)
    hp$id<-rownames(hp)
    
    hp$chr.absolute<-NA
    hp$start.absolute<-NA
    hp$end.absolute<-NA
    hp$chr.hpreg.absolute<-NA
    hp$start.hpreg.absolute<-NA
    hp$end.hpreg.absolute<-NA
    
    if(!(exists("hp.initial.input.type"))){
      log("Hairpinprimer input type must be a regions file!")
      stop("Hairpinprimer input type must be a regions file!")
    }
    
    if(hp.initial.input.type == "regions"){
      
      for (iicor in 1:nrow(hp)){
        
        iicor.rowid<-paste(hp[iicor,"sequenceID"])
        iicor.input.id<-paste(hp[iicor,"inputID"])
        iicor.absolute.start<-as.numeric(as.character(bed[bed$sequenceID == iicor.input.id,"start"]))
        iicor.absolute.end<-as.numeric(as.character(bed[bed$sequenceID == iicor.input.id,"end"]))
        
        hp[iicor,"chr.absolute"]<-paste(bed[bed$sequenceID == iicor.input.id,"chr"])
        hp[iicor,"start.absolute"] <- iicor.absolute.start
        hp[iicor,"end.absolute"] <- iicor.absolute.end
        
        hp[iicor,"chr.hpreg.absolute"]<-paste(bed[bed$sequenceID == iicor.input.id,"chr"])
        hp[iicor,"start.hpreg.absolute"] <- as.numeric(iicor.absolute.start + as.numeric(as.character(hp[iicor,"hp_Region_Start_Relative"])) - 1)
        hp[iicor,"end.hpreg.absolute"] <- as.numeric(iicor.absolute.start + as.numeric(as.character(hp[iicor,"hp_Region_End_Relative"])) - 1)
        
        
      }#for iicor
    }#if(hp.initial.input.type=="regions")
    
    write.table(hp,file=paste(path.wd,analysis.id,"_Hairpinizer_Results_all.txt",sep=""),
                sep="\t",dec=".",col.names=TRUE,row.names=FALSE,quote=F)

################################################################################

#pre filter for primer design.
    hpf<-hp[!is.na(hp$hp_Region_Length),]
    hpf<-hpf[hpf$hp_Region_Length>=hp.length.min & hpf$hp_Region_Length<=hp.length.max,]
    hpf$hp_Sequence<-as.character(hpf$hp_Sequence)
    hpf$linker_Sequence<-as.character(hpf$linker_Sequence)
    
    hpf$LinkerStartInMolecule<-nchar(hpf$hp_Sequence)+1
    hpf$LinkerEndInMolecule<-hpf$LinkerStartInMolecule+nchar(hpf$linker_Sequence)-1
 
################################################################################

#export temporal results
    log("Export temporal results.")
    write.table(hpf,file=paste(path.wd,analysis.id,"_Hairpinizer_Results_PreFiltered.txt",sep=""),
                sep="\t",dec=".",col.names=TRUE,row.names=FALSE,quote=F)
 
################################################################################

#Filter for target regions...
    hpf2<-hpf

    #if(!is.na(target.start.column.id) & !is.na(target.end.column.id)){
    #  log("Filter for fragments that cover target region.")
    #  
    #  for (iib in 1:nrow(bed)){
    #    iib.bed.id<-paste(bed[iib,"sequenceID"])
    #    iib.start<-bed[iib,target.start.column.id]
    #    iib.end<-bed[iib,target.end.column.id]
    #    ida<-hpf2[hpf2$inputID==iib.bed.id,]
    #    hpf2<-hpf2[hpf2$inputID!=iib.bed.id,]    
    #    idaf<-ida[ida$hp_Region_Start_Relative<=iib.start & ida$hp_Region_End_Relative>=iib.end,]
    #    
    #    if(nrow(idaf)>0){
    #      hpf2<-rbind(hpf2,idaf)
    #    }#if                      
    #  }# for iib
    #}# if !is.na target.start.column.id

################################################################################

#if necessary re-assign target start/end column ids.
target.start.column.id<-"LinkerStartInMolecule"
target.end.column.id<-"LinkerEndInMolecule"
use.target.columns<-TRUE

################################################################################

#export of results
    log("Export temporal results.")
    write.table(hpf2,file=paste(path.wd,analysis.id,"_Hairpinizer_Results_Final.txt",sep=""),
                sep="\t",dec=".",col.names=TRUE,row.names=FALSE,quote=F)

    #html.report(filenames =  paste(path.html,analysis.id,"_Hairpinizer_Results_Final.txt",sep=""),
                #filename.out = paste(path.html,"Hairpinizer_Results_",analysis.id,".html",sep=""),txt.header = TRUE,txt.sep = "\t")

    log("Done.")

################################################################################

    log("Bedfile swapping...")
    write.table(bed,file=paste(path.sequences,analysis.id,"_PrimerDesign_Initial.InputSequences.txt",sep=""),sep="\t",dec=".",
                col.names=TRUE,row.names=FALSE,quote=F)
    #html.report(filenames =  paste(path.sequences,analysis.id,"PrimerDesign_Initial.InputSequences.txt",sep=""),
    #            filename.out = paste(path.html,"bedfile_",analysis.id,".html",sep=""),txt.header = TRUE,txt.sep = "\t")
    log("Done.")
    
    bed.ini<-bed    
    
    bed<-hpf2
    
    write.table(bed,file=paste(path.sequences,analysis.id,"InputSequences_AfterBedFileSwap.txt",sep=""),sep="\t",dec=".",
                col.names=TRUE,row.names=FALSE,quote=F)
    log("Done.")
  
    if(nrow(bed)<1){
      log(paste("After filtering no hairpin sequences left for primer design."))
      stop()
    }
    
################################################################################    
    
    #perform input sequence characterization on hairpin regions
    log("Start characterization of input sequences...")
    
    for (iiii in 1:nrow(bed)){
      
      seqanai<-sequence.characterization(paste(bed[iiii,"sequence"]))
      write.table(seqanai,file=paste(path.sequences,"sequence.characterization_",bed[iiii,"sequenceID"],".txt",sep=""),
                  sep="\t",dec=".",col.names=T,row.names = FALSE,quote=F)
      
      seq.ana.all[[paste(bed[iiii,"sequenceID"])]]<-seqanai
    }
    
    log("Done.")
    
  }#if hp_bisulfite or hp_NOMe

################################################################################

#strand selection
if (strand=="both"){
  strand<-c("top","bottom")
}

# disable strand selection for hairpin analyses
if(primer.type=="hp_bisulfite" | primer.type=="hp_NOME" | primer.type=="hp_genomic"){
    strand<-"top"
    log("Enforced strand-setting to top-strand since it is a hairpin analysis.")
    }

###################################################################################

#primer.design section
index<-1
#np<-data.frame()
for (istrand in strand){ #top and bottom strand
  
  for (iseq in 1:nrow(bed)){ #for each input seq
    
    sequence<-as.character(bed[iseq,"sequence"])   
    sequence.id<-as.character(bed[iseq,"sequenceID"])
    
    if(primer.type=="bisulfite" | primer.type=="NOME" | primer.type=="genomic" | primer.type=="CLEVER" | primer.type == "CrispRCas9PCR"){
      
    if(input.type=="regions"){
      i.amp.chr<-as.numeric(paste(gsub("chr","",bed[iseq,"chr"])))
      i.amp.start<-as.numeric(paste(bed[iseq,"start"]))
      i.amp.end<-as.numeric(paste(bed[iseq,"end"]))
    }
    if(input.type=="sequences"){
      i.amp.chr<-NA
      i.amp.start<-NA
      i.amp.end<-NA
    }
    }#if(primer.type=="hp_bisulfite" | primer.type=="hp_NOME" | primer.type=="hp_genomic"){
    
    if(primer.type=="hp_bisulfite" | primer.type=="hp_NOME" | primer.type=="hp_genomic" | primer.type=="hp_CLEVER"){
      
      if(hp.initial.input.type  =="sequences"){
        
        i.amp.chr<-NA
        i.amp.start<-NA
        i.amp.end<-NA

        log(paste("Enforced Input.type setting to 'sequence' since it is a hairpin analysis."))
      }
      
      if(hp.initial.input.type=="regions"){
        
        input.type<-"sequences"
        i.amp.chr<-as.numeric(paste(gsub("chr","",bed[iseq,"chr.hpreg.absolute"])))
        i.amp.start<-as.numeric(paste(bed[iseq,"start.hpreg.absolute"]))
        i.amp.end<-as.numeric(paste(bed[iseq,"end.hpreg.absolute"]))
        
        log(paste("Enforced Input.type setting to 'sequence' since it is a hairpin analysis."))
        
      }
    }

    log(paste(index,": Start primer design for amplicon ",sequence.id," [",istrand," Strand]...",sep=""))
    
    #####################################Calling the various primer designs according to primer type########################################### 
    if(primer.type=="NOME"){
    np<-nome.primer.design(sequence=sequence,
                           sequence.id=sequence.id,
                           min.Tm.primer=min.Tm.primer,
                           max.Tm.primer=max.Tm.primer,
                           min.number.gc.amplicon=min.number.gc.amplicon,
                           min.number.cg.amplicon=min.number.cg.amplicon,
                           max.Tm.difference.primer=max.Tm.difference.primer,
                           primer.align.binsize=primer.align.binsize,
                           min.length.primer=min.length.primer,
                           max.length.primer=max.length.primer, 
                           low.complexity.primer.removal=low.complexity.primer.removal,
                           max.bins.low.complexity=max.bins.low.complexity,
                           remove.primers.with.n=remove.primers.with.n,
                           min.C2T.primer1=min.C2T.primer1,
                           min.G2A.primer2=min.G2A.primer2,
                           min.length.amplicon=min.length.amplicon,
                           max.length.amplicon=max.length.amplicon,
                           strand=istrand,
                           mode=mode)
    }
    
    if(primer.type=="bisulfite"){
      np<-bisulfite.primer.design(sequence=sequence,
                             sequence.id=sequence.id,
                             min.Tm.primer=min.Tm.primer,
                             max.Tm.primer=max.Tm.primer,
                             min.number.gc.amplicon=min.number.gc.amplicon,
                             min.number.cg.amplicon=min.number.cg.amplicon,
                             max.Tm.difference.primer=max.Tm.difference.primer,
                             primer.align.binsize=primer.align.binsize,
                             min.length.primer=min.length.primer,
                             max.length.primer=max.length.primer, 
                             low.complexity.primer.removal=low.complexity.primer.removal,
                             max.bins.low.complexity=max.bins.low.complexity,
                             remove.primers.with.n=remove.primers.with.n,
                             min.C2T.primer1=min.C2T.primer1,
                             min.G2A.primer2=min.G2A.primer2,
                             min.length.amplicon=min.length.amplicon,
                             max.length.amplicon=max.length.amplicon,
                             strand=istrand,
                             mode=mode)
      colnames(np)<-gsub("NOME","bisulfite",colnames(np))
    }
    
    if(primer.type=="genomic"){
      np<-genomic.primer.design(sequence=sequence,
                                  sequence.id=sequence.id,
                                  min.Tm.primer=min.Tm.primer,
                                  max.Tm.primer=max.Tm.primer,
                                  min.number.gc.amplicon=min.number.gc.amplicon,
                                  min.number.cg.amplicon=min.number.cg.amplicon,
                                  max.Tm.difference.primer=max.Tm.difference.primer,
                                  primer.align.binsize=primer.align.binsize,
                                  min.length.primer=min.length.primer,
                                  max.length.primer=max.length.primer, 
                                  low.complexity.primer.removal=low.complexity.primer.removal,
                                  max.bins.low.complexity=max.bins.low.complexity,
                                  remove.primers.with.n=remove.primers.with.n,
                                  min.C2T.primer1=0,
                                  min.G2A.primer2=0,
                                  min.length.amplicon=min.length.amplicon,
                                  max.length.amplicon=max.length.amplicon,
                                  strand=istrand,
                                  mode=mode,
                                  chop.size=chop.size)
      colnames(np)<-gsub("NOME","genomic",colnames(np))
    }
    
    if(primer.type=="CLEVER"){
      np<-CLEVER.primer.design(sequence=sequence,
                                sequence.id=sequence.id,
                                min.Tm.primer=min.Tm.primer,
                                max.Tm.primer=max.Tm.primer,
                                min.number.gc.amplicon=min.number.gc.amplicon,
                                min.number.cg.amplicon=min.number.cg.amplicon,
                                max.Tm.difference.primer=max.Tm.difference.primer,
                                primer.align.binsize=primer.align.binsize,
                                min.length.primer=min.length.primer,
                                max.length.primer=max.length.primer, 
                                low.complexity.primer.removal=low.complexity.primer.removal,
                                max.bins.low.complexity=max.bins.low.complexity,
                                remove.primers.with.n=remove.primers.with.n,
                                min.C2T.primer1=0,
                                min.G2A.primer2=0,
                                min.length.amplicon=min.length.amplicon,
                                max.length.amplicon=max.length.amplicon,
                                strand=istrand,
                                mode=mode)
      colnames(np)<-gsub("NOME","CLEVER",colnames(np))
    }
    
    #
    
    if(primer.type=="CrispRCas9PCR"){
      np<-crispR.Cas9.amp.primer.design(sequence=sequence,
                                        sequence.id=sequence.id,
                                        min.Tm.primer=min.Tm.primer,
                                        max.Tm.primer=max.Tm.primer,
                                        min.number.gc.amplicon=min.number.gc.amplicon,
                                        min.number.cg.amplicon=min.number.cg.amplicon,
                                        max.Tm.difference.primer=max.Tm.difference.primer,
                                        primer.align.binsize=primer.align.binsize,
                                        min.length.primer=min.length.primer,
                                        max.length.primer=max.length.primer, 
                                        low.complexity.primer.removal=low.complexity.primer.removal,
                                        max.bins.low.complexity=max.bins.low.complexity,
                                        remove.primers.with.n=remove.primers.with.n,
                                        min.C2T.primer1=0,
                                        min.G2A.primer2=0,
                                        min.length.amplicon=min.length.amplicon,
                                        max.length.amplicon=max.length.amplicon,
                                        strand=istrand,
                                        mode=mode)
      colnames(np)<-gsub("NOME","CrispRCas9PCR",colnames(np))
    }
#
      if(primer.type=="hp_bisulfite"){
      np<-bisulfite.primer.design(sequence=sequence,
                                  sequence.id=sequence.id,
                                  min.Tm.primer=min.Tm.primer,
                                  max.Tm.primer=max.Tm.primer,
                                  min.number.gc.amplicon=min.number.gc.amplicon,
                                  min.number.cg.amplicon=min.number.cg.amplicon,
                                  max.Tm.difference.primer=max.Tm.difference.primer,
                                  primer.align.binsize=primer.align.binsize,
                                  min.length.primer=min.length.primer,
                                  max.length.primer=max.length.primer, 
                                  low.complexity.primer.removal=low.complexity.primer.removal,
                                  max.bins.low.complexity=max.bins.low.complexity,
                                  remove.primers.with.n=remove.primers.with.n,
                                  min.C2T.primer1=min.C2T.primer1,
                                  min.G2A.primer2=min.G2A.primer2,
                                  min.length.amplicon=min.length.amplicon,
                                  max.length.amplicon=max.length.amplicon,
                                  strand=istrand,
                                  mode=mode)
      colnames(np)<-gsub("NOME","bisulfite",colnames(np))
      }

      if(primer.type=="hp_NOME"){
      np<-nome.primer.design(sequence=sequence,
                             sequence.id=sequence.id,
                             min.Tm.primer=min.Tm.primer,
                             max.Tm.primer=max.Tm.primer,
                             min.number.gc.amplicon=min.number.gc.amplicon,
                             min.number.cg.amplicon=min.number.cg.amplicon,
                             max.Tm.difference.primer=max.Tm.difference.primer,
                             primer.align.binsize=primer.align.binsize,
                             min.length.primer=min.length.primer,
                             max.length.primer=max.length.primer, 
                             low.complexity.primer.removal=low.complexity.primer.removal,
                             max.bins.low.complexity=max.bins.low.complexity,
                             remove.primers.with.n=remove.primers.with.n,
                             min.C2T.primer1=min.C2T.primer1,
                             min.G2A.primer2=min.G2A.primer2,
                             min.length.amplicon=min.length.amplicon,
                             max.length.amplicon=max.length.amplicon,
                             strand=istrand,
                             mode=mode)  
      }
    
    if(primer.type=="hp_genomic"){
      np<-genomic.primer.design(sequence=sequence,
                                  sequence.id=sequence.id,
                                  min.Tm.primer=min.Tm.primer,
                                  max.Tm.primer=max.Tm.primer,
                                  min.number.gc.amplicon=min.number.gc.amplicon,
                                  min.number.cg.amplicon=min.number.cg.amplicon,
                                  max.Tm.difference.primer=max.Tm.difference.primer,
                                  primer.align.binsize=primer.align.binsize,
                                  min.length.primer=min.length.primer,
                                  max.length.primer=max.length.primer, 
                                  low.complexity.primer.removal=low.complexity.primer.removal,
                                  max.bins.low.complexity=max.bins.low.complexity,
                                  remove.primers.with.n=remove.primers.with.n,
                                  min.C2T.primer1=0,
                                  min.G2A.primer2=0,
                                  min.length.amplicon=min.length.amplicon,
                                  max.length.amplicon=max.length.amplicon,
                                  strand=istrand,
                                  mode=mode,
                                  chop.size=chop.size)
      colnames(np)<-gsub("NOME","genomic",colnames(np))
    }
    
    if(primer.type=="hp_CLEVER"){
      np<-CLEVER.primer.design(sequence=sequence,
                                  sequence.id=sequence.id,
                                  min.Tm.primer=min.Tm.primer,
                                  max.Tm.primer=max.Tm.primer,
                                  min.number.gc.amplicon=min.number.gc.amplicon,
                                  min.number.cg.amplicon=min.number.cg.amplicon,
                                  max.Tm.difference.primer=max.Tm.difference.primer,
                                  primer.align.binsize=primer.align.binsize,
                                  min.length.primer=min.length.primer,
                                  max.length.primer=max.length.primer, 
                                  low.complexity.primer.removal=low.complexity.primer.removal,
                                  max.bins.low.complexity=max.bins.low.complexity,
                                  remove.primers.with.n=remove.primers.with.n,
                                  min.C2T.primer1=0,
                                  min.G2A.primer2=0,
                                  min.length.amplicon=min.length.amplicon,
                                  max.length.amplicon=max.length.amplicon,
                                  strand=istrand,
                                  mode=mode)
      colnames(np)<-gsub("NOME","CLEVER",colnames(np))
    }
    
    print(paste("Index: ",index,"   nrow.np: ",nrow(np)))
    log("Primer design tool finished.")

if(is.na(np[1,1])){
  
  log(paste("No primer pair found for: ", sequence.id,sep=""))
  next()
  
}
    
############################################################################
  
# primer.design results post-processing
    log(paste("Process primer design intermediate result (",nrow(np)," primer pairs)...",sep=""))
    log(paste("Current amplicon located on  ",i.amp.chr,sep=""))
    
if(nrow(np)>0){
  
  np$amplicon.chr<-i.amp.chr
  
    if(primer.type=="hp_bisulfite" | primer.type=="hp_NOME" | primer.type=="hp_genomic" | primer.type=="hp_CLEVER"){
      
      if(hp.initial.input.type=="sequences"){
        np$amplicon.start<-NA
        np$amplicon.end<-NA      
      }
      
      if(hp.initial.input.type=="regions"){
        np$amplicon.start<-as.numeric(i.amp.start)+as.numeric(np$amplicon.start.relative)-1
        np$amplicon.end<-as.numeric(np$amplicon.start)+as.numeric(np$amplicon.length)-1  
      }
      
    }#if(primer.type=="hp_bisulfite" | primer.type=="hp_NOME" | primer.type=="hp_genomic" | primer.type=="hp_CLEVER"){
  
    log("Done.")

#if(!is.na(np[1,"amplicon.length"])){
  if(!exists("results") & !is.na(np[1,1])){
    
    results<-np
    log("Export temporal results...")
    write.table(results,file=paste(path.temp,analysis.id,"_TempResults.",sequence.id,".txt",sep=""),
    col.names=T,row.names=T,sep="\t",dec=".",quote=F)
    log("Done.")
  }# if exists results
 
  if(exists("results") & !is.na(np[1,1])){
    write.table(np,file=paste(path.temp,analysis.id,"_TempResults.",sequence.id,".txt",sep=""),
                col.names=T,row.names=T,sep="\t",dec=".",quote=F)
    
    results<-rbind(results,np[,colnames(np) %in% colnames(results)])
    
    log("Export all temporal results...")
    write.table(results,file=paste(path.temp,analysis.id,"_TempResults.all.txt",sep=""),
                col.names=T,row.names=T,sep="\t",dec=".",quote=F)
    log("Done.")
    
  }#else
  }# if nrow np>0
    
index<-index+1

  }# for isec
}#istrand

###################################################################################

if(exists(x="results")){

#if there are results, they will be processed
log("Final processing of complete results...")
  
results2<-as.data.frame(results)

#if using regions mode then absolute positions will be added
if(input.type=="regions"){
  
  log("Add absolute primer positions to primer table...")

  results2$amplicon.chr<-iregs[match(results2$sequence.id,iregs$sequenceID),"chr"]
  results2$amplicon.start<-as.numeric(iregs[match(results2$sequence.id,iregs$sequenceID),"start"])+results2$amplicon.start.relative-1
  results2$amplicon.end<-results2$amplicon.start+as.numeric(results2$amplicon.length)-1
  results2$assembly<-iregs[match(results2$sequence.id,iregs$sequenceID),"assembly"]
  results2$primer1.start<-as.character(results2$amplicon.start)
  results2$primer1.end<-as.character(as.numeric(results2$primer1.start)+as.numeric(results2$primer1.length)-1)
  results2$primer2.start<-as.character(as.numeric(results2$amplicon.end)-as.numeric(results2$primer2.length)+1)
  results2$primer2.end<-as.character(results2$amplicon.end)
  results2$amplicon.adress<-paste(results2$amplicon.chr,":",results2$amplicon.start,"-",results2$amplicon.end,sep="")
}

if(primer.type=="hp_bisulfite" | primer.type=="hp_NOME" | primer.type=="hp_genomic" | primer.type=="hp_CLEVER"){
#if using regions mode then absolute positions will be added
if(hp.initial.input.type=="regions"){
  
  log("Add absolute primer positions to primer table...")
  
  results2$amplicon.chr<-paste(bed[match(results2$sequence.id,bed$sequenceID),"chr.hpreg.absolute"])
  results2$amplicon.start<-as.numeric(bed[match(results2$sequence.id,bed$sequenceID),"start.hpreg.absolute"])+
    results2$amplicon.start.relative-1
  results2$amplicon.end<-as.numeric(bed[match(results2$sequence.id,bed$sequenceID),"start.hpreg.absolute"])+
    nchar(bed[match(results2$sequence.id,bed$sequenceID),"sequence"])-
    nchar(bed[match(results2$sequence.id,bed$sequenceID),"linker_Sequence"])-1
  results2$assembly<-bed[match(results2$sequence.id,bed$sequenceID),"assembly"]
  results2$primer1.start<-as.character(results2$amplicon.start)
  results2$primer1.end<-as.character(as.numeric(results2$primer1.start)+as.numeric(results2$primer1.length)-1)
  results2$primer2.start<-as.character(as.numeric(results2$amplicon.end)-as.numeric(results2$primer2.length)+1)
  results2$primer2.end<-as.character(results2$amplicon.end)
  results2$amplicon.adress<-paste(results2$amplicon.chr,":",results2$amplicon.start,"-",results2$amplicon.end,sep="")
  
}
}#if(primer.type=="hp_bisulfite" | primer.type=="hp_NOME" | primer.type=="hp_genomic"){

results2$primer.type<-primer.type
#results2$mode<-mode

##################################################################################

#export
log(paste("Export final temporal results..."))
write.table(results2,file=paste(path.temp,analysis.id,"_TempResults.final.txt",sep=""),
            col.names=T,row.names=T,sep="\t",dec=".",quote=F)
log("Done.")

####################################################################################

#Filter na long amplicons.
results2<-results2[!is.na(results2$amplicon.length),]
if(nrow(results2)>0){
  
###################################################################################

if (check4snps){
  
  log("Assign SNP information to designed amplicons...")
  
  if(exists("x") && x == "all_my_snps"){
    
    results2$SNP.db<-NA
    results2$amplicon.nSNPs<-NA
    results2$primer1.nSNPs<-NA
    results2$primer2.nSNPs<-NA
    results2$amplicon.SNP.ids<-NA
    results2$primer1.SNP.ids<-NA
    results2$primer2.SNP.ids<-NA
    
    for (iamplicons in 1:nrow(results2)){
    
      iseq.id<-paste(results2[iamplicons,"sequence.id"])
      iseq.chr<-paste(results2[iamplicons,"amplicon.chr"])
      iseq.amp.start<-results2[iamplicons,"amplicon.start"]
      iseq.amp.end<-results2[iamplicons,"amplicon.end"]
      iseq.p1.start<-results2[iamplicons,"primer1.start"]
      iseq.p1.end<-results2[iamplicons,"primer1.end"]
      iseq.p2.start<-results2[iamplicons,"primer2.start"]
      iseq.p2.end<-results2[iamplicons,"primer2.end"]
      
      if(exists("all_my_snps") && nrow(all_my_snps) != 0 && length(levels(all_my_snps) != 0) ){
        results2[iamplicons,"SNP.db"] <- levels(all_my_snps[gsub("chr","",as.character(all_my_snps$chr))==as.character(gsub("chr","",iseq.chr)) & 
                                                            as.numeric(as.character(all_my_snps$start))>=iseq.amp.start & 
                                                            as.numeric(as.character(all_my_snps$start))<=iseq.amp.end ,"source"])[1]
      
        results2[iamplicons,"amplicon.nSNPs"] <- nrow(all_my_snps[gsub("chr","",as.character(all_my_snps$chr))==gsub("chr", "", iseq.chr) & 
                                                                  as.numeric(as.character(all_my_snps$start))>=iseq.amp.start & 
                                                                  as.numeric(as.character(all_my_snps$start))<=iseq.amp.end ,]) 
      
        results2[iamplicons,"amplicon.SNP.ids"]<-paste(all_my_snps[(gsub("chr","",as.character(all_my_snps$chr))==as.character(gsub("chr", "", iseq.chr))) & 
                                                                   (as.numeric(as.character(all_my_snps$start))>=iseq.amp.start) & 
                                                                   (as.numeric(as.character(all_my_snps$start))<=iseq.amp.end) ,"rs_id"],collapse=",")
      
      
        results2[iamplicons,"primer1.nSNPs"] <- nrow(all_my_snps[gsub("chr","",as.character(all_my_snps$chr))==as.character(gsub("chr", "", iseq.chr)) & 
                                                                 as.numeric(as.character(all_my_snps$start))>=iseq.p1.start & 
                                                                 as.numeric(as.character(all_my_snps$start))<=iseq.p1.end ,])
      
        results2[iamplicons,"primer1.SNP.ids"]<-paste(all_my_snps[(gsub("chr","",as.character(all_my_snps$chr))==as.character(gsub("chr", "", iseq.chr))) & 
                                                                  (as.numeric(as.character(all_my_snps$start))>=iseq.p1.start) & 
                                                                  (as.numeric(as.character(all_my_snps$start))<=iseq.p1.end) ,"rs_id"],collapse=",")
      
        results2[iamplicons,"primer2.nSNPs"] <- nrow(all_my_snps[gsub("chr","",as.character(all_my_snps$chr))==as.character(gsub("chr", "", iseq.chr)) & 
                                                                 as.numeric(as.character(all_my_snps$start))>=iseq.p2.start & 
                                                                 as.numeric(as.character(all_my_snps$start))<=iseq.p2.end,])
      
        results2[iamplicons,"primer2.SNP.ids"]<-paste(all_my_snps[(gsub("chr","",as.character(all_my_snps$chr))==as.character(gsub("chr", "", iseq.chr))) & 
                                                                  (as.numeric(as.character(all_my_snps$start))>=iseq.p2.start) & 
                                                                  (as.numeric(as.character(all_my_snps$start))<=iseq.p2.end) ,"rs_id"],collapse=",")
        }
      
      }# iamplicons
      
      log("Done.")
      
      ##################################################

    #export
    log(paste("Export final temporal results..."))
    write.table(results2,file=paste(path.temp,analysis.id,"_TempResults.final_SNP_all.txt",sep=""),
                col.names=T,row.names=T,sep="\t",dec=".",quote=F)
    log("Done.")

  #####FILTERSTEPS
  
  log("Filter amplicons/primers by SNP info...")
  results2 <- results2[results2$amplicon.nSNPs>=min.snps.amplicon & results2$amplicon.nSNPs<=max.snps.amplicon,]
  results2 <- results2[results2$primer1.nSNPs>=min.snps.primer1 & results2$primer1.nSNPs<=max.snps.primer1,]
  results2 <- results2[results2$primer2.nSNPs>=min.snps.primer2 & results2$primer2.nSNPs<=max.snps.primer2,]
  log(paste("done: ",nrow(results2)," amplicons survived.",sep=""))
  
  #export
  log(paste("Export final temporal results..."))
  write.table(results2,file=paste(path.temp,analysis.id,"_TempResults.final_SNP_filtered.txt",sep=""),
              col.names=T,row.names=T,sep="\t",dec=".")
  log("Done.")
  log(paste("Finished SNP module..."))
  
  }#if exists
  
  if(!exists("all_my_snps")){
    
  log(paste("WARNING: No SNP info found for SNP analysis...skipped"))  
    
    results2$SNPdb<-NA
    results2$amplicon.nSNPs<-NA
    results2$primer1.nSNPs<-NA
    results2$primer2.nSNPs<-NA
    results2$amplicon.SNP.ids<-NA
    results2$primer1.SNP.ids<-NA
    results2$primer2.SNP.ids<-NA
    
  }#if!exists
  
}#if check SNPs

###################################################################################
###################################################################################

if (check4repeats){
  
  log("Assign repeat information to designed amplicons...")
  
  if(exists(x = "all_my_repeats")){
    
    #results2$repeat.db<-"Repeatmasker"
    #results2$amplicon.n.repeats<-NA
    #results2$primer1.n.repeats<-NA
    #results2$primer2.n.repeats<-NA
    #results2$amplicon.repeat.ids<-NA
    #results2$primer1.repeat.ids<-NA
    #results2$primer2.repeat.ids<-NA
    
    if(nrow(results2) >= 1){}
    for (iamplicons in 1:nrow(results2)){
      
      iseqid<-paste(results2[iamplicons,"sequence.id"])
      iseq.chr<-paste(results2[iamplicons,"amplicon.chr"])
      iseq.amp.start<-results2[iamplicons,"amplicon.start"]
      iseq.amp.end<-results2[iamplicons,"amplicon.end"]
      iseq.p1.start<-results2[iamplicons,"primer1.start"]
      iseq.p1.end<-results2[iamplicons,"primer1.end"]
      iseq.p2.start<-results2[iamplicons,"primer2.start"]
      iseq.p2.end<-results2[iamplicons,"primer2.end"]
      
      
      p1.reps<-all_my_repeats[as.character(all_my_repeats$chr)==as.character(gsub("chr","",iseq.chr)) &
                          (as.numeric(as.character(all_my_repeats$start))>=iseq.p1.start & 
                             as.numeric(as.character(all_my_repeats$start))<=iseq.p1.end) |
                          (as.numeric(as.character(all_my_repeats$start))<=iseq.p1.start & 
                             as.numeric(as.character(all_my_repeats$end))>=iseq.p1.start) |
                          (as.numeric(as.character(all_my_repeats$start))>=iseq.p1.start & as.numeric(as.character(all_my_repeats$start))<=iseq.p1.end),]
      
      p2.reps<-all_my_repeats[as.character(all_my_repeats$chr)==as.character(gsub("chr","",iseq.chr)) &
                                (as.numeric(as.character(all_my_repeats$start))>=iseq.p2.start & as.numeric(as.character(all_my_repeats$start))<=iseq.p2.end) |
                                (as.numeric(as.character(all_my_repeats$start))<=iseq.p2.start & as.numeric(as.character(all_my_repeats$end))>=iseq.p2.start) |
                                (as.numeric(as.character(all_my_repeats$start))>=iseq.p2.start & as.numeric(as.character(all_my_repeats$start))<=iseq.p2.end),]
      
      amp.reps<-all_my_repeats[as.character(all_my_repeats$chr)==as.character(gsub("chr","",iseq.chr)) &
                                (as.numeric(as.character(all_my_repeats$start))>=iseq.amp.start & as.numeric(as.character(all_my_repeats$start))<=iseq.amp.end) |
                                (as.numeric(as.character(all_my_repeats$start))<=iseq.amp.start & as.numeric(as.character(all_my_repeats$end))>=iseq.amp.start) |
                                (as.numeric(as.character(all_my_repeats$start))>=iseq.amp.start & as.numeric(as.character(all_my_repeats$start))<=iseq.amp.end),]
      
      results2[iamplicons,"amplicon.n.repeats"] <- nrow(amp.reps)
      results2[iamplicons,"amplicon.repeat.ids"] <- paste(amp.reps[,"repClass"],collapse=",")
      results2[iamplicons,"primer1.n.repeats"] <- nrow(p1.reps)
      results2[iamplicons,"primer1.repeat.ids"] <- paste(p1.reps[,"repClass"],collapse=",")
      results2[iamplicons,"primer2.n.repeats"] <- nrow(p2.reps)
      results2[iamplicons,"primer2.repeat.ids"] <- paste(p2.reps[,"repClass"],collapse=",")
      
    }# iamplicons
    }
    
    log("Done.")
    
    ##################################################
    
    #export
    log(paste("Export final temporal results..."))
    write.table(results2,file=paste(path.temp,analysis.id,"_TempResults.final_Repeats_all.txt",sep=""),
                col.names=T,row.names=T,sep="\t",dec=".",quote=F)
    log("Done.")
    
    #####FILTERSTEPS
    
    if(allow.repeats.in.primers==FALSE){
      
      log("Filter amplicons/primers by repeat info on primers...")
      results2 <- results2[results2$primer1.n.repeats == 0 & results2$primer2.n.repeats == 0,]
      log("Done.")
      
    }
    
    if(allow.repeats.in.amplicon==FALSE){
      
      log("Filter amplicons/primers by repeat info on amplicon...")
      results2 <- results2[results2$amplicon.n.repeats == 0,]
      log("Done.")
    
    }
    
    log(paste("done: ",nrow(results2)," amplicons survived.",sep=""))
    
    #export
    log(paste("Export final temporal results..."))
    write.table(results2,file=paste(path.temp,analysis.id,"_TempResults.final_Repeats_filtered.txt",sep=""),
                col.names=T,row.names=T,sep="\t",dec=".")
    log("Done.")
    log(paste("Finished repeat module..."))
    
  }#if exists
  
  if(!exists("all_my_repeats")){
    
    log(paste("WARNING: No info found for repeat analysis...skipped"))  
    
    results2$repeat.db<-NA
    results2$amplicon.n.repeats<-NA
    results2$primer1.n.repeats<-NA
    results2$primer2.n.repeats<-NA
    results2$amplicon.repeat.ids<-NA
    results2$primer1.repeat.ids<-NA
    results2$primer2.repeat.ids<-NA
    
  }#if!exists
  
}#if check repeats

  
###################################################################################
###################################################################################
  
  if (annotate.genes){
    
    log("Assign gene information to designed amplicons...")
    
    if(exists(x = "all_my_genes") && is.data.frame(results2) && nrow(results2) > 0){
      
      results2$gene.db<-NA
      results2$amplicon.n.genes<-NA
      results2$amplicon.gene.ids<-NA

      for (iamplicons in 1:nrow(results2)){
        
        iseqid<-paste(results2[iamplicons,"sequence.id"])
        iseq.chr<-paste(results2[iamplicons,"amplicon.chr"])
        iseq.amp.start<-results2[iamplicons,"amplicon.start"]
        iseq.amp.end<-results2[iamplicons,"amplicon.end"]
        iseq.p1.start<-results2[iamplicons,"primer1.start"]
        iseq.p1.end<-results2[iamplicons,"primer1.end"]
        iseq.p2.start<-results2[iamplicons,"primer2.start"]
        iseq.p2.end<-results2[iamplicons,"primer2.end"]
        
        amp.genes<-all_my_genes[all_my_genes$chrom==iseq.chr &
                                 (all_my_genes$txStart>=iseq.amp.start & all_my_genes$txStart<=iseq.amp.end) |
                                 (all_my_genes$txStart<=iseq.amp.start & all_my_genes$txEnd>=iseq.amp.start) |
                                 (all_my_genes$txStart>=iseq.amp.start & all_my_genes$txStart<=iseq.amp.end),]
        
        results2[iamplicons,"amplicon.n.genes"] <- nrow(amp.genes)
        results2[iamplicons,"amplicon.gene.ids"] <- paste(amp.genes[,"gene_name"],collapse=",")
        
      }# iamplicons
      
      log("Done.")
      
      ##################################################
      
      #export
      log(paste("Export final temporal results..."))
      write.table(results2,file=paste(path.temp,analysis.id,"_TempResults.final_Genes_all.txt",sep=""),
                  col.names=T,row.names=T,sep="\t",dec=".",quote=F)
      log("Done.")
      
      log(paste("Finished Gene annotation module..."))
      
    }#if exists
    
    if(!exists("all_my_genes")){
      
      log(paste("WARNING: No info found for Gene annotation...skipped"))  
      
      results2$gene.db<-NA
      results2$amplicon.n.genes<-NA
      results2$amplicon.gene.ids<-NA
      
    }#if!exists
    
  }#if annotate genes

###################################################################################
###################################################################################

if (annotate.cpg.islands){
  
  log("Assign cpgi information to designed amplicons...")
  
  if(exists(x = "all_my_cpgis")){
    
    results2$cpgi.db<-"UCSC"
    results2$amplicon.n.cpgi<-NA
    results2$amplicon.cpgi.ids<-NA
    
    for (iamplicons in 1:nrow(results2)){
      
      iseqid<-paste(results2[iamplicons,"sequence.id"])
      iseq.chr<-paste(results2[iamplicons,"amplicon.chr"])
      iseq.amp.start<-results2[iamplicons,"amplicon.start"]
      iseq.amp.end<-results2[iamplicons,"amplicon.end"]
      iseq.p1.start<-results2[iamplicons,"primer1.start"]
      iseq.p1.end<-results2[iamplicons,"primer1.end"]
      iseq.p2.start<-results2[iamplicons,"primer2.start"]
      iseq.p2.end<-results2[iamplicons,"primer2.end"]
      
      amp.cpgis<-all_my_cpgis[all_my_cpgis$chrom==iseq.chr &
                                (all_my_cpgis$chromStart>=iseq.amp.start & all_my_cpgis$chromStart<=iseq.amp.end) |
                                (all_my_cpgis$chromStart<=iseq.amp.start & all_my_cpgis$chromEnd>=iseq.amp.start) |
                                (all_my_cpgis$chromStart>=iseq.amp.start & all_my_cpgis$chromStart<=iseq.amp.end),]
      
      results2[iamplicons,"amplicon.n.cpgi"] <- nrow(amp.cpgis)
      results2[iamplicons,"amplicon.cpgi.ids"] <- paste(amp.cpgis[,"name"],collapse=",")
      
    }# iamplicons
    
    log("Done.")
    
    ##################################################
    
    #export
    log(paste("Export final temporal results..."))
    write.table(results2,file=paste(path.temp,analysis.id,"_TempResults.final_cpgi_all.txt",sep=""),
                col.names=T,row.names=T,sep="\t",dec=".",quote=F)
    log("Done.")
    
    log(paste("Finished CpG Island annotation module..."))
    
  }#if exists
  
  if(!exists("all_my_cpgis")){
    
    log(paste("WARNING: No info found for CpG Island annotation...skipped"))  
    
    results2$cpgi.db<-NA
    results2$amplicon.n.cpgi<-NA
    results2$amplicon.cpgi.ids<-NA
    
  }#if!exists
  
}#if annotate cpgis

###################################################################################  
###################################################################################
  
#Remove results that do not cover target region.
if(use.target.columns){
  remcount<-1
  
  for (iseq in 1:nrow(bed)){
    log("Remove amplicons that not cover specified target region...")
    iamptarem<-bed[iseq,"sequenceID"]
    
    if(allow.primer.in.target.region == FALSE){
      target.start<-as.numeric(bed[iseq,target.start.column.id])-max.length.primer
      target.end<-as.numeric(bed[iseq,target.end.column.id])+max.length.primer
    }
    
    if(allow.primer.in.target.region){
      target.start<-as.numeric(bed[iseq,target.start.column.id])
      target.end<-as.numeric(bed[iseq,target.end.column.id])
    }
    
    if(!is.na(target.start) & !is.na(target.end)){
    tempres<-results2[results2$sequence.id==iamptarem,]
    fi.tg<-tempres$primer1.end.relative <= target.start & tempres$primer2.start.relative >= target.end
    fi.notg<-tempres$primer1.end.relative > target.start & tempres$primer2.start.relative < target.end
    survivors<-length(fi.tg[fi.tg==TRUE])
    deaths<-length(fi.tg[fi.tg==FALSE]) 
    log(paste("For ",iamptarem," ",deaths," amplicons were removed and ",survivors," were retained.",sep=""))
           
        if(remcount==1){
          ret<-tempres[fi.tg,]
          rem<-tempres[fi.notg,]
        }
        if(remcount>1){
          ret<-rbind(ret,tempres[fi.tg,])
          rem<-rbind(rem,tempres[fi.notg,])
        }
        remcount<-remcount+1
      }
  }
  
  write.table(rem,file=paste(path.temp,analysis.id,"_Primers_offtargetregion.txt",sep=""),
              col.names=T,row.names=T,sep="\t",dec=".")
  results2<-ret   
log("Done.")
}

##########################################################################################

#fix v9.8: remove double results
  temp.p1p2.ids<-paste0(results2$amplicon.id,"_",results2$primer1.length,"_",results2$primer2.length)
  results2<-results2[!duplicated(temp.p1p2.ids),]
  write.table(results2,file=paste(path.temp,analysis.id,"_TempResults.final_after.dup.removal.txt",sep=""),
              col.names=T,row.names=T,sep="\t",dec=".",quote=F)

##########################################################################################
 
  #create columns with added NGS adaptors to primers. 
  
  if(!is.null(add.ngs.adaptor.f)){
  
    log("Add NGS adaptors for primer 1...")
    results2$primer1.ngs.added<-paste(add.ngs.adaptor.f,results2$primer1.sequence,sep="")
    log("Done.")
  }
  
  if(is.null(add.ngs.adaptor.f)){
    
    log("Add NGS adaptors for primer 1...")
    results2$primer1.ngs.added<-NA
    log("Done.")
  }
  
  if(!is.null(add.ngs.adaptor.r)){
    
    log("Add NGS adaptors for primer 2...")
    results2$primer2.ngs.added<-paste(add.ngs.adaptor.r,results2$primer2.sequence,sep="")
    log("Done.")
  
  }
  
  if(is.null(add.ngs.adaptor.r)){
    
    log("Add NGS adaptors for primer 2...")
    results2$primer2.ngs.added<-NA
    log("Done.")
    
  }
  
###################################################################################

  #export results
log("Export results...")

if(file.exists(paste(path.wd,"primer_",analysis.id,"_wholelist",".txt",sep=""))){
  primer.old<-read.table(paste(path.wd,"primer_",analysis.id,"_wholelist",".txt",sep=""))
  results2<-rbind(primer.old,results2)
  write.table(results2,file=paste(path.wd,"primer_",analysis.id,"_wholelist",".txt",sep=""),
              col.names=TRUE,row.names=FALSE,sep="\t",dec=".") 
}
if(!file.exists(paste(path.wd,"primer_",analysis.id,"_wholelist",".txt",sep=""))){
  write.table(results2,file=paste(path.wd,"primer_",analysis.id,"_wholelist",".txt",sep=""),
              col.names=TRUE,row.names=FALSE,sep="\t",dec=".") 
}
log("Done.")

log("Create HTML...")
n.total.primer<-nrow(results2[!is.na(results2$amplicon.id),])
# html.report(filenames =  paste(path.html,"primer_",analysis.id,".txt",sep=""),
            # filename.out = paste(path.html,"primer_",analysis.id,".html",sep=""),txt.header = TRUE,txt.sep = "\t")
log("Done.")

###################################################################################

#for all amplicons generate .fasta file
afn<-paste(path.sequences,results2$amplicon.id,".FASTA",sep="")
mapply(FUN = write.fasta,results2$amplicon.sequence.genomic,results2$amplicon.id,afn)
results3<-results2

###################################################################################

#create toplist tables
if (create.toplist){
  log("Calculate toplist tables...")
  tl<-results3[!is.na(results3$amplicon.id),]
  al<-levels(factor(tl$sequence.id))
  
  for (ial in 1:length(al)){
    ta<-al[ial]
    pool<-tl[tl$sequence.id==ta,]
    if(nrow(pool)==1){   
      top.cand<-pool[1,]
      if(ial==1){
        toplist<-top.cand
      }
      else{
        toplist<-rbind(toplist,top.cand)
      }  
    }
    else{
      if(is.infinite(min(pool$primer1.nSNPs))){
        pool2<-pool
      }
      else{
        pool2<-pool[which.min(pool$primer1.nSNPs),]
      }
      if(is.infinite(min(pool2$primer2.nSNPs))){
        pool3<-pool2
      }
      else{
        pool3<-pool2[which.min(pool2$primer2.nSNPs),]
      }
      if(primer.type=="NOME" | primer.type=="hp_NOME"){
      if(is.infinite(max(pool3$nGCs))){
        pool4<-pool3
      }
      else{
        pool4<-pool3[which.max(pool3$nGCs),]
      }
      }
      
      if(primer.type=="bisulfite" | primer.type=="hp_bisulfite" | primer.type=="CLEVER" | primer.type=="hp_CLEVER" | primer.type=="CrispRCas9PCR"){
        if(is.infinite(max(pool3$nCGs))){
          pool4<-pool3
        }
        else{
          pool4<-pool3[which.max(pool3$nCGs),]
        }
      }
      
      if(primer.type=="genomic" | primer.type=="hp_genomic"){
        if(is.infinite(max(pool3$nCGs))){
          pool4<-pool3
        }
        else{
          pool4<-pool3[which.max(pool3$nCGs),]
        }
      }
      
      if(is.infinite(min(pool4$amplicon.length))){
        pool5<-pool4
      }
      else{
        pool5<-pool4[which.min(pool4$amplicon.length),]
      }
      
    lv.str<-levels(factor(pool5$dna.strand))
    if ("top" %in% lv.str){
      pool6<-pool5[pool5$dna.strand=="top",]
    }
    else{
      pool6<-pool5[pool5$dna.strand=="bottom"]
    }
    pool7<-pool6[which.min(pool6$primer1.length),]
    pool8<-pool7[which.min(pool7$primer2.length),]
    
    top.cand<-pool8[1,]
    
    if(ial==1){
      toplist<-top.cand
    }
    else{
      toplist<-rbind(toplist,top.cand)
    }
  }
  }
  toplist<-as.data.frame(toplist)
log("Done.")

###################################################################################

if(nrow(toplist)>0){
  log("Export toplist tables...")
  
  if(file.exists(paste(path.wd,"primer_",analysis.id,"_toplist.txt",sep=""))){
    top.old<-read.table(paste(path.wd,"primer_",analysis.id,"_toplist.txt",sep=""),header=TRUE,sep="\t",dec=".")
    toplist<-rbind(top.old,toplist)
    write.table(toplist,file=paste(path.wd,"primer_",analysis.id,"_toplist.txt",sep=""),
                col.names=TRUE,row.names=FALSE,sep="\t",dec=".") 
  }
  
  if(!file.exists(paste(path.wd,"primer_",analysis.id,"_toplist.txt",sep=""))){
    write.table(toplist,file=paste(path.wd,"primer_",analysis.id,"_toplist.txt",sep=""),
                col.names=TRUE,row.names=FALSE,sep="\t",dec=".") 
  }
  
  #html.report(filenames =  paste(path.html,"primer_",analysis.id,"_toplist.txt",sep=""),filename.out = paste(path.html,"primer_",analysis.id,"_toplist.html",sep=""),txt.header = TRUE,txt.sep = "\t")
  log("Done.")
  }

}

###################################################################################

#generate ucsc browser custom tracks for each sequence
if(input.type=="regions"){  
  
  log("Create UCSC custom tracks...")
  
  #for complete primer list.
  results2$ucsc.color<-700
  results2$nl<-"\n"
  
  fa<-levels(factor(results2$amplicon.id))[1]
  faf<-results2$amplicon.id==fa
  broch<-paste(results2[faf,"amplicon.chr"][1])
  bromi<-min(results2[faf,"amplicon.start"])-1000
  broma<-max(results2[faf,"amplicon.end"])+1000
  ffn<-paste(path.tracks,analysis.id,"_ucsc_tracks_all.txt",sep="")
  
  write.text(paste('browser position ',broch,':',bromi,'-',broma,'\n',
                   'browser hide all','\n',
                   'track name="Amplicons" description="Primerpairs" visibility=2 color=0,128,0 useScore=1','\n',sep=''),
             file=ffn) 
  
  write.text(apply(X = results2[,c("amplicon.chr","amplicon.start","amplicon.end","primer.pair.id","ucsc.color","nl")],1,paste,sep="\t"),
             file=ffn,add=TRUE)
  
  
  if (create.toplist && !(is.data.frame(toplist) && nrow(toplist)==0)){ 
  
    #for primer toplist.
    toplist$ucsc.color<-700
    toplist$nl<-"\n"
  
    fa<-levels(factor(toplist$amplicon.id))[1]
    faf<-toplist$amplicon.id==fa
    broch<-paste(toplist[faf,"amplicon.chr"][1])
    bromi<-min(toplist[faf,"amplicon.start"])-1000
    broma<-max(toplist[faf,"amplicon.end"])+1000
    ffn<-paste(path.tracks,analysis.id,"_ucsc_tracks_toplist.txt",sep="")
  
    write.text(paste('browser position ',broch,':',bromi,'-',broma,'\n',
                     'browser hide all','\n',
                     'track name="Amplicons" description="Primerpairs" visibility=2 color=0,128,0 useScore=1','\n',sep=''),
               file=ffn) 
    
    write.text(apply(X = toplist[,c("amplicon.chr","amplicon.start","amplicon.end","primer.pair.id","ucsc.color","nl")],1,paste,sep="\t"),
               file=ffn,add=TRUE)
  }# if (create.toplist)
  
  print("Done.")
}
###################################################################################

#generate custom tracks for BiSearch tool.
log("Create BiSearch custom tracks...")
  
  write.table(results2[,c("primer.pair.id","primer1.sequence","primer2.sequence")], 
              file=paste(path.tracks,analysis.id,"_BiSearch_tracks_all.txt",sep=""),
              sep=" ",quote = FALSE,col.names=F,row.names=F)
  
  if(create.toplist){
    write.table(toplist[,c("primer.pair.id","primer1.sequence","primer2.sequence")], 
                file=paste(path.tracks,analysis.id,"_BiSearch_tracks_toplist.txt",sep=""),
                sep=" ",quote = FALSE,col.names=F,row.names=F)
  }

print("Done.")

###################################################################################

#generate .bed files for amplicon and/or primers.
if(input.type=="regions"){ 
  
  log("Create .amplicon and primer .bed files...")

  write.table(results2[,c("amplicon.chr","primer1.start","primer1.end", "primer.pair.id")], 
              file=paste(path.tracks,analysis.id,"_bed_primer1_all.bed",sep=""),sep="\t",quote = FALSE,col.names=F,row.names=F)
  write.table(results2[,c("amplicon.chr","primer2.start","primer2.end", "primer.pair.id")], 
              file=paste(path.tracks,analysis.id,"_bed_primer2_all.bed",sep=""),sep="\t",quote = FALSE,col.names=F,row.names=F) 
  write.table(results2[,c("amplicon.chr","amplicon.start","amplicon.end", "primer.pair.id")], 
              file=paste(path.tracks,analysis.id,"_bed_amplicon_all.bed",sep=""),sep="\t",quote = FALSE,col.names=F,row.names=F)
  
  if(create.toplist){
  write.table(toplist[,c("amplicon.chr","primer1.start","primer1.end", "primer.pair.id")], 
              file=paste(path.tracks,analysis.id,"_bed_primer1_toplist.bed",sep=""),sep="\t",quote = FALSE,col.names=F,row.names=F)
  write.table(toplist[,c("amplicon.chr","primer2.start","primer2.end", "primer.pair.id")], 
              file=paste(path.tracks,analysis.id,"_bed_primer2_toplist.bed",sep=""),sep="\t",quote = FALSE,col.names=F,row.names=F) 
  write.table(toplist[,c("amplicon.chr","amplicon.start","amplicon.end", "primer.pair.id")], 
              file=paste(path.tracks,analysis.id,"_bed_amplicon_toplist.bed",sep=""),sep="\t",quote = FALSE,col.names=F,row.names=F)
    
   }
  
  print("Done.")
}#if(input.type=="regions"){

###################################################################################

log("Generate summary...")
end.run<-paste(Sys.time())
smry<-data.frame(PARAMETERS=c("Analysis Start","Analysis End","Number of Input Sequences",
            "Number of PrimerPairs [Total]","Number of Sequences with Primers","Number of Sequences without Primers"),
            SETTINGS=c(start.run,end.run, paste(nrow(bed)),paste(n.total.primer),
                       paste(nrow(toplist)),paste(nrow(bed)-nrow(toplist))))
write.table(smry,file=paste(path.wd,"summary_",analysis.id,".txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",dec=".")
#html.report(filenames =  paste(path.html,"summary_",analysis.id,".txt",sep=""),filename.out = paste(path.html,"summary_",analysis.id,".html",sep=""),txt.header = TRUE,txt.sep = "\t")

# primer designs by input sequence
if(!(is.na(results2$sequence.id))){
  tbl<-data.frame(table(results2$sequence.id))
  colnames(tbl)<-c("sequence.id","amplicons[n]")
  write.table(tbl,file=paste(path.wd,"primerdesigns.by.sequence_",analysis.id,".txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",dec=".")
  #html.report(filenames =  paste(path.html,"primerdesigns.by.sequence_",analysis.id,".txt",sep=""),filename.out = paste(path.html,"summary_",analysis.id,".html",sep=""),txt.header = TRUE,txt.sep = "\t")
}

log("Done.")

###################################################################################

log("Write whitelist...")
white<-bed[paste(bed$sequenceID) %in% paste(toplist[,1]),]

  if(!file.exists(paste(path.wd,"primer_", analysis.id, "_whitelist.txt", sep=""))){
    white.old<-read.table(paste(path.wd, "primer_", analysis.id, "_whitelist.txt", sep=""), header=TRUE, sep="\t",dec=".")
    
    
    white<-rbind(white.old, white) 
    write.table(white,file=paste(path.wd,"primer_",analysis.id,"_whitelist.txt",sep=""),
                 col.names=TRUE,row.names=FALSE,sep="\t",dec=".") 
  } else {
    write.table(white,file=paste(path.wd,"primer_",analysis.id,"_whitelist.txt",sep=""),
                col.names=TRUE,row.names=FALSE,sep="\t",dec=".") 
  }

log("Done.")

log("Write blacklist...")
black<-bed[!paste(bed$sequenceID) %in% paste(white$sequenceID),]
write.table(black,file=paste(path.wd,"primer_",analysis.id,"_blacklist.txt",sep=""),
            col.names=TRUE,row.names=FALSE,sep="\t",dec=".") 
log("Done.")

###################################################################################
###################################################################################
if(create.graphics){
  
  log("Create graphics...")
  require(ggplot2)
  
  if(length(levels(factor(results2$sequence.id))) >= 1){
  
  for(ibp in 1:length(levels(factor(results2$sequence.id)))){
    
    ibps<-paste(levels(factor(results2$sequence.id)))[ibp]
    plot.dat<-results2[as.character(results2$sequence.id)==ibps,]
    
    if(create.amplicon.barplots){
    
    png(filename=paste(path.graphics,"bar.amplicon.length_",ibps,".png",sep=""),height=1000,width=1200,pointsize=24)
    par(mar=c(8,4,4,2)+0.1)
     barplot(plot.dat[,"amplicon.length"],
            names.arg = paste(plot.dat$amplicon.id),las=3,cex.names = 0.5,
            ylab="amplicon length",
            main=paste(ibps,": size of designed amplicons [n=",nrow(plot.dat),"]",sep=""),
            cex.main=0.8)
    dev.off()
    
    png(filename=paste(path.graphics,"bar.number.cgs_",ibps,".png",sep=""),height=1000,width=1200,pointsize=24)
    par(mar=c(8,4,4,2)+0.1)
    barplot(plot.dat[,"nCGs"],
            names.arg = paste(plot.dat$amplicon.id),las=3,cex.names = 0.5,
            ylab="number of CpGs",
            main=paste(ibps,": occurence of CpG sites in designed amplicons [n=",nrow(plot.dat),"]",sep=""),
            cex.main=0.8)
    dev.off()
    
    png(filename=paste(path.graphics,"bar.number.gcs_",ibps,".png",sep=""),height=1000,width=1200,pointsize=24)
    par(mar=c(8,4,4,2)+0.1)
    barplot(plot.dat[,"nGCs"],
            names.arg = paste(plot.dat$amplicon.id),las=3,cex.names = 0.5,
            ylab="number of GpCs",
            main=paste(ibps,": occurence of GpC sites in designed amplicons [n=",nrow(plot.dat),"]",sep=""),
            cex.main=0.8)
    dev.off()
    }
    
 ###############################################################################################################
    
    if(primer.type == "bisulfite" | primer.type == "NOME" | primer.type == "genomic" | primer.type == "CLEVER" | primer.type=="CrispRCas9PCR"){
    
    #per sequence id primer design overview plot
    
    #bed[ir,c("start","end","assembly","sequenceID","sequence.length","sequence.adress","sequence")]
    ipt<-NA
    bed.length<-nchar(as.character(bed[as.character(bed$sequenceID)==ibps,"sequence"]))
    sa<-seq.ana.all[[ibps]]#read.table(file="E:\\Science_Jil\\Projects\\SB_PrimerDesign\\#LiveDemo\\html_bis_withSNPanalysis\\sequences\\sequence.characterization_sequence1.txt",header=T,sep="\t",dec=".")
    sels<-results2[as.character(results2$sequence.id)==ibps,]
    
    }# if(not hairpin)
    
    #hp.initial.input.type <- input.type
    
    if((primer.type == "hp_bisulfite" | primer.type == "hp_NOME" | primer.type == "hp_genomic" | primer.type == "hp_CLEVER") &&
       hp.initial.input.type == "regions"){
         
         #per sequence id primer design overview plot
         
         #bed[ir,c("start","end","assembly","sequenceID","sequence.length","sequence.adress","sequence")]
         ipt<-paste(hpf2[hpf2$sequenceID==ibps,"inputID"])
         bed.length<-nchar(as.character(hpf2[as.character(hpf2$sequenceID)==ibps,"sequence"]))
         sa<-seq.ana.all[[ibps]]#read.table(file="E:\\Science_Jil\\Projects\\SB_PrimerDesign\\#LiveDemo\\html_bis_withSNPanalysis\\sequences\\sequence.characterization_sequence1.txt",header=T,sep="\t",dec=".")
         sels<-results2[as.character(results2$sequence.id)==ibps,]
    }#if hp
    
    lol<-data.frame(relative.position=rep(c(1,nrow(sa)),6),
                    cg.pos=NA,
                    gc.pos=NA,
                    snp_start=NA,
                    snp_end=NA,
                    p1_start=NA,
                    p1_end=NA,
                    p2_start=NA,
                    p2_end=NA,
                    amp.start=NA,
                    amp.end=NA,
                    linker.start=NA,
                    linker.end=NA,
                    repeat.start=NA,
                    repeat.end=NA,
                    gene.start=NA,
                    gene.end=NA,
                    cpgi.start=NA,
                    cpgi.end=NA,
                    sequence.id=ibps,
                    amplicon.id=c(rep("input.sequence [CpG]",2),
                                  rep("input.sequence [GpC]",2),
                                  rep("input.sequence [SNP]",2),
                                  rep("input.sequence [Repeats]",2),
                                  rep("input.sequence [CpG Island]",2),
                                  rep("input.sequence [Genes]",2)),
                    feature=NA)
    
    if((primer.type == "hp_bisulfite" | primer.type == "hp_NOME" | primer.type == "hp_genomic" | primer.type == "hp_CLEVER") &&
       hp.initial.input.type == "regions"){
      
      lol<-rbind(lol,data.frame(relative.position=rep(c(1,nrow(sa)),2),
                                cg.pos=NA,
                                gc.pos=NA,
                                snp_start=NA,
                                snp_end=NA,
                                p1_start=NA,
                                p1_end=NA,
                                p2_start=NA,
                                p2_end=NA,
                                amp.start=NA,
                                amp.end=NA,
                                linker.start=NA,
                                linker.end=NA,
                                repeat.start=NA,
                                repeat.end=NA,
                                gene.start=NA,
                                gene.end=NA,
                                cpgi.start=NA,
                                cpgi.end=NA,
                                sequence.id=ibps,
                                amplicon.id=rep("input.sequence [linker]",2),
                                feature=NA))
    }# linker line
    
    lol<-rbind(lol,data.frame(relative.position=NA,
                              cg.pos=sa[sa$cg==TRUE,"basecount"],
                              gc.pos=NA,
                              snp_start=NA,
                              snp_end=NA,
                              p1_start=NA,
                              p1_end=NA,
                              p2_start=NA,
                              p2_end=NA,
                              amp.start=NA,
                              amp.end=NA,
                              linker.start=NA,
                              linker.end=NA,
                              repeat.start=NA,
                              repeat.end=NA,
                              gene.start=NA,
                              gene.end=NA,
                              cpgi.start=NA,
                              cpgi.end=NA,
                              sequence.id=ibps,
                              amplicon.id="input.sequence [CpG]",
                              feature="CpG"))
    
    lol<-rbind(lol,data.frame(relative.position=NA,
                              cg.pos=NA,
                              gc.pos=sa[sa$gc==TRUE,"basecount"],
                              snp_start=NA,
                              snp_end=NA,
                              p1_start=NA,
                              p1_end=NA,
                              p2_start=NA,
                              p2_end=NA,
                              amp.start=NA,
                              amp.end=NA,
                              linker.start=NA,
                              linker.end=NA,
                              repeat.start=NA,
                              repeat.end=NA,
                              gene.start=NA,
                              gene.end=NA,
                              cpgi.start=NA,
                              cpgi.end=NA,
                              sequence.id=ibps,
                              amplicon.id="input.sequence [GpC]",
                              feature="GpC"))
    
    if((primer.type == "hp_bisulfite" | primer.type == "hp_NOME" | primer.type == "hp_genomic" | primer.type == "hp_CLEVER") &&
      hp.initial.input.type == "regions"){
      
    lol<-rbind(lol,data.frame(relative.position=NA,
                              cg.pos=NA,
                              gc.pos=NA,
                              snp_start=NA,
                              snp_end=NA,
                              p1_start=NA,
                              p1_end=NA,
                              p2_start=NA,
                              p2_end=NA,
                              amp.start=NA,
                              amp.end=NA,
                              linker.start=as.numeric(as.character(hpf2[hpf2$sequenceID==ibps,"LinkerStartInMolecule"])),
                              linker.end=as.numeric(as.character(hpf2[hpf2$sequenceID==ibps,"LinkerEndInMolecule"])),
                              repeat.start=NA,
                              repeat.end=NA,
                              gene.start=NA,
                              gene.end=NA,
                              cpgi.start=NA,
                              cpgi.end=NA,
                              sequence.id=ibps,
                              amplicon.id="input.sequence [linker]",
                              feature="linker"))
    }#if hp
    
    lol<-rbind(lol,data.frame(relative.position=NA,
                              cg.pos=NA,
                              gc.pos=NA,
                              snp_start=NA,
                              snp_end=NA,
                              p1_start=as.numeric(as.character(sels$primer1.start.relative)),
                              p1_end=as.numeric(as.character(sels$primer1.end.relative)),
                              p2_start=as.numeric(as.character(sels$primer2.start.relative)),
                              p2_end=as.numeric(as.character(sels$primer2.end.relative)),
                              amp.start=NA,
                              amp.end=NA,
                              linker.start=NA,
                              linker.end=NA,
                              repeat.start=NA,
                              repeat.end=NA,
                              gene.start=NA,
                              gene.end=NA,
                              cpgi.start=NA,
                              cpgi.end=NA,
                              sequence.id=sels$sequence.id,
                              amplicon.id=sels$amplicon.id,
                              feature="primer"))
    
    lol<-rbind(lol,data.frame(relative.position=NA,
                              cg.pos=NA,
                              gc.pos=NA,
                              snp_start=NA,
                              snp_end=NA,
                              p1_start=NA,
                              p1_end=NA,
                              p2_start=NA,
                              p2_end=NA,
                              amp.start=as.numeric(as.character(sels$amplicon.start.relative)),
                              amp.end=as.numeric(as.character(sels$amplicon.end.relative)),
                              linker.start=NA,
                              linker.end=NA,
                              repeat.start=NA,
                              repeat.end=NA,
                              gene.start=NA,
                              gene.end=NA,
                              cpgi.start=NA,
                              cpgi.end=NA,
                              sequence.id=sels$sequence.id,
                              amplicon.id=sels$amplicon.id,
                              feature=NA))
    
    if(check4snps & exists("all_my_snps")){
      
      if((primer.type == "bisulfite" | primer.type == "NOME" | primer.type == "genomic" | primer.type == "CLEVER" | primer.type=="CrispRCas9PCR") &&
         input.type == "regions"){
        
      selchr<-paste(unique(sels[as.character(sels$sequence.id)==ibps,"amplicon.chr"]))  
      bedstart<-as.numeric(bed[bed$sequenceID==ibps,"start"])
      bedend<-as.numeric(bed[bed$sequenceID==ibps,"end"])
      
      }#if not hp
      
      
      if((primer.type == "hp_bisulfite" | primer.type == "hp_NOME" | primer.type == "hp_genomic" | primer.type == "hp_CLEVER") &&
         hp.initial.input.type == "regions"){
        
      selchr<-paste(unique(hpf2[as.character(hpf2$sequenceID)==ibps,"chr.absolute"]))  
      bedstart<-as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"start.absolute"]))
      bedend<-as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"start.absolute"])) + 
        as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"LinkerStartInMolecule"])) - 1
      lnk.eim<-as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"LinkerEndInMolecule"]))
      
      }#if hp
      
      
      tolo<-all_my_snps[gsub("chr","",all_my_snps$chr)==gsub("chr","",selchr) & 
                          as.numeric(as.character(all_my_snps$start)) >= bedstart & 
                          as.numeric(as.character(all_my_snps$end)) <= bedend,]
      
      tolo$start.relative=as.numeric(as.character(tolo$start))-bedstart+1
      tolo$end.relative=as.numeric(as.character(tolo$end))-bedstart+1
      
      if((primer.type == "hp_bisulfite" | primer.type == "hp_NOME" | primer.type == "hp_genomic" | primer.type == "hp_CLEVER") &
         hp.initial.input.type == "regions"){
        
      #trick to also show SNPs on the second part of the hairpin
      tolo.temp<-tolo
      tolo.temp$end.relative <- lnk.eim + as.numeric(as.character(tolo$start.relative)) - 1 
      tolo.temp$start.relative <- lnk.eim + as.numeric(as.character(tolo$end.relative)) - 1
      tolo<-rbind(tolo,tolo.temp)
      
      }#if hp
      
      #write.table(tolo,file=paste(path.tracks,"tolo.object.for.plots_temp_",ibps,".txt",sep=""),
      #            sep="\t",dec=".",col.names=T,row.names=F,quote=F)
      
      srs<-as.numeric(as.character(tolo$start.relative))
      sre<-as.numeric(as.character(tolo$end.relative))
      
      if(length(srs)==0){srs<-NA}
      if(length(sre)==0){sre<-NA}
      
      lol<-rbind(lol,data.frame(relative.position=NA,
                                cg.pos=NA,
                                gc.pos=NA,
                                snp_start=srs,
                                snp_end=sre,
                                p1_start=NA,
                                p1_end=NA,
                                p2_start=NA,
                                p2_end=NA,
                                amp.start=NA,
                                amp.end=NA,
                                linker.start=NA,
                                linker.end=NA,
                                repeat.start=NA,
                                repeat.end=NA,
                                gene.start=NA,
                                gene.end=NA,
                                cpgi.start=NA,
                                cpgi.end=NA,
                                sequence.id=ibps,
                                amplicon.id="input.sequence [SNP]",
                                feature="SNP"))
    }# check4snps

#############################################################################################################
    
    if(check4repeats & exists("all_my_repeats")){
      
      if((primer.type == "bisulfite" | primer.type == "NOME" | primer.type == "genomic" | primer.type == "CLEVER" | primer.type=="CrispRCas9PCR") &&
         input.type == "regions"){
        
        selchr<-paste(unique(sels[as.character(sels$sequence.id)==ibps,"amplicon.chr"]))  
        bedstart<-as.numeric(bed[bed$sequenceID==ibps,"start"])
        bedend<-as.numeric(bed[bed$sequenceID==ibps,"end"])
        
      }#if not hp
      
      
      if((primer.type == "hp_bisulfite" | primer.type == "hp_NOME" | primer.type == "hp_genomic" | primer.type == "hp_CLEVER") &&
         hp.initial.input.type == "regions"){
        
        selchr<-paste(unique(hpf2[as.character(hpf2$sequenceID)==ibps,"chr.absolute"]))  
        bedstart<-as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"start.absolute"]))
        bedend<-as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"start.absolute"])) + 
          as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"LinkerStartInMolecule"])) - 1
        
        lnk.eim<-as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"LinkerEndInMolecule"]))
        
      }#if hp
      
      tolo<-all_my_repeats[as.character(all_my_repeats$chr)==gsub("chr","",selchr),] #& 
                          #all_my_repeats$genoStart >= bedstart & 
                          #all_my_repeats$genoEnd <= bedend
      
      tolo$start.relative=as.numeric(as.character(tolo$start))-bedstart+1
      tolo$end.relative=as.numeric(as.character(tolo$end))-bedstart+1
      
      #tolo<-tolo[(tolo$start.relative>=0 & tolo$end.relative<=bed.length) |
      #             ((tolo$start.relative<0 & tolo$end.relative<=bed.length) & tolo$end.relative>=0) |
      #             (tolo$start.relative>=0 & tolo$start.relative<=bed.length) |
      #             (tolo$end.relative>=0 & tolo$end.relative<=bed.length),]
      
      
      tolo<-tolo[(!as.numeric(as.character(tolo$end.relative))<0) &
                  (!as.numeric(as.character(tolo$start.relative))>=bed.length),]
      
        if((primer.type == "hp_bisulfite" | primer.type == "hp_NOME" | primer.type == "hp_genomic" | primer.type == "hp_CLEVER") &&
         hp.initial.input.type == "regions"){
        
        #trick to also show SNPs on the second part of the hairpin
        tolo.temp<-tolo
        tolo.temp$end.relative <- lnk.eim + as.numeric(as.character(tolo$start.relative)) - 1 
        tolo.temp$start.relative <- lnk.eim + as.numeric(as.character(tolo$end.relative)) - 1
        tolo<-rbind(tolo,tolo.temp)
        
      }#if hp
      
      if(nrow(tolo[as.numeric(as.character(tolo$start.relative)) < 0,])>0){
        tolo[as.numeric(as.character(tolo$start.relative)) < 0,"start.relative"]<-0
      }
      
      if(nrow(tolo[as.numeric(as.character(tolo$end.relative)) < 0,])>0){
        tolo[as.numeric(as.character(tolo$end.relative)) > bed.length,"end.relative"]<-bed.length
      }
      
      #write.table(tolo,file=paste(path.tracks,"tolo.object.for.plots_temp_",ibps,".txt",sep=""),
      #            sep="\t",dec=".",col.names=T,row.names=F,quote=F)
      
      srs<-as.numeric(as.character(tolo$start.relative))
      sre<-as.numeric(as.character(tolo$end.relative))
      
      if(length(srs)==0){srs<-NA}
      if(length(sre)==0){sre<-NA}
      
      #lol$repeat.start<-NA
      #lol$repeat.end<-NA
      
      lol<-rbind(lol,data.frame(relative.position=NA,
                                cg.pos=NA,
                                gc.pos=NA,
                                snp_start=NA,
                                snp_end=NA,
                                p1_start=NA,
                                p1_end=NA,
                                p2_start=NA,
                                p2_end=NA,
                                amp.start=NA,
                                amp.end=NA,
                                linker.start=NA,
                                linker.end=NA,
                                repeat.start=sre,
                                repeat.end=srs,
                                gene.start=NA,
                                gene.end=NA,
                                cpgi.start=NA,
                                cpgi.end=NA,
                                sequence.id=ibps,
                                amplicon.id="input.sequence [Repeats]",
                                feature=NA))
    }# check4repeats

########################################################################################################    

    if(annotate.genes & exists("all_my_genes")){
      
      if((primer.type == "bisulfite" | primer.type == "NOME" | primer.type == "genomic" | primer.type == "CLEVER" | primer.type=="CrispRCas9PCR") &
         input.type == "regions"){
        
        selchr<-paste(unique(sels[as.character(sels$sequence.id)==ibps,"amplicon.chr"]))  
        bedstart<-as.numeric(bed[bed$sequenceID==ibps,"start"])
        bedend<-as.numeric(bed[bed$sequenceID==ibps,"end"])
        
      }#if not hp
      
      
      if((primer.type == "hp_bisulfite" | primer.type == "hp_NOME" | primer.type == "hp_genomic" | primer.type == "hp_CLEVER") &&
         hp.initial.input.type == "regions"){
        
        selchr<-paste(unique(hpf2[as.character(hpf2$sequenceID)==ibps,"chr.absolute"]))  
        bedstart<-as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"start.absolute"]))
        bedend<-as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"start.absolute"])) + 
          as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"LinkerStartInMolecule"])) - 1
        lnk.eim<-as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"LinkerEndInMolecule"]))
        
      }#if hp
      
      tolo<-all_my_genes[all_my_genes$chrom==selchr ,]
      
      gns<-tolo
      all.exons.start<-as.numeric(unlist(strsplit(as.character(gns$exonStarts),",")))
      all.exons.end<-as.numeric(unlist(strsplit(as.character(gns$exonEnds),",")))
      
      if(length(all.exons.start) == length(all.exons.end)){
        
        add.n.rows<-length(all.exons.start)-nrow(tolo)
        adNA<-as.data.frame(matrix(NA,nrow = add.n.rows,ncol=ncol(tolo)))
        colnames(adNA)<-colnames(tolo)
        
        tolo<-rbind(tolo,adNA)
        tolo$exon.start=all.exons.start
        tolo$exon.end=all.exons.end
        
        tolo$start.relative=as.numeric(as.character(tolo$exon.start))-bedstart+1
        tolo$end.relative=as.numeric(as.character(tolo$exon.end))-bedstart+1
        
       # tolo<-tolo[(tolo$start.relative>=0 & tolo$end.relative<=bed.length) |
      #             (tolo$start.relative<0 & tolo$end.relative<=bed.length & tolo$end.relative>=0) |
      #             (tolo$start.relative>=0 & tolo$start.relative<=bed.length) |
      #             (tolo$end.relative>=0 & tolo$end.relative<=bed.length),]
        
        tolo<-tolo[(!as.numeric(as.character(tolo$end.relative))<0) &
                     (!as.numeric(as.character(tolo$start.relative))>=bed.length),]
        
      }#if(length(all.exons.start) == length(all.exons.end)){
      
      if(length(all.exons.start) != length(all.exons.end)){
        
      tolo$start.relative=tolo$txStart-bedstart+1
      tolo$end.relative=tolo$txEnd-bedstart+1
      
      }#if(length(all.exons.start) != length(all.exons.end)){
      
      
      if((primer.type == "hp_bisulfite" | primer.type == "hp_NOME" | primer.type == "hp_genomic" | primer.type == "hp_CLEVER") &&
         hp.initial.input.type == "regions"){
        
        #trick to also show SNPs on the second part of the hairpin
        tolo.temp<-tolo
        tolo.temp$end.relative <- lnk.eim + tolo$start.relative - 1 
        tolo.temp$start.relative <- lnk.eim + tolo$end.relative - 1
        tolo<-rbind(tolo,tolo.temp)
        
      }#if hp
      
      if(nrow(tolo[as.numeric(as.character(tolo$start.relative)) < 0,])>0 && !is.na(tolo[as.numeric(as.character(tolo$start.relative)) <0,"start.relative"])){
        tolo[as.numeric(as.character(tolo$start.relative)) < 0,"start.relative"]<-0
      }
      
      if(nrow(tolo[as.numeric(as.character(tolo$end.relative)) < 0,])>0 && !is.na(tolo[as.numeric(as.character(tolo$end.relative)) > bed.length,"end.relative"])){
        tolo[as.numeric(as.character(tolo$end.relative)) > bed.length,"end.relative"]<-bed.length
      }
      
      #write.table(tolo,file=paste(path.tracks,"tolo.object.for.plots_temp_",ibps,".txt",sep=""),
      #            sep="\t",dec=".",col.names=T,row.names=F,quote=F)
      
      srs<-as.numeric(as.character(tolo$start.relative))
      sre<-as.numeric(as.character(tolo$end.relative))
      
      if(length(srs)==0){srs<-NA}
      if(length(sre)==0){sre<-NA}
      
      lol<-rbind(lol,data.frame(relative.position=NA,
                                cg.pos=NA,
                                gc.pos=NA,
                                snp_start=NA,
                                snp_end=NA,
                                p1_start=NA,
                                p1_end=NA,
                                p2_start=NA,
                                p2_end=NA,
                                amp.start=NA,
                                amp.end=NA,
                                linker.start=NA,
                                linker.end=NA,
                                repeat.start=NA,
                                repeat.end=NA,
                                gene.start=srs,
                                gene.end=sre,
                                cpgi.start=NA,
                                cpgi.end=NA,
                                sequence.id=ibps,
                                amplicon.id="input.sequence [Genes]",
                                feature=NA))
    }# annotate.genes

##############################################################################################################
    
    if(annotate.cpg.islands & exists("all_my_cpgis")){
      
      if((primer.type == "bisulfite" | primer.type == "NOME" | primer.type == "genomic" | primer.type == "CLEVER" | primer.type=="CrispRCas9PCR") &
         input.type == "regions"){
        
        selchr<-paste(unique(sels[as.character(sels$sequence.id)==ibps,"amplicon.chr"]))  
        bedstart<-as.numeric(bed[bed$sequenceID==ibps,"start"])
        bedend<-as.numeric(bed[bed$sequenceID==ibps,"end"])
        
      }#if not hp
      
      
      if((primer.type == "hp_bisulfite" | primer.type == "hp_NOME" | primer.type == "hp_genomic" | primer.type == "hp_CLEVER") &&
         hp.initial.input.type == "regions"){
        
        selchr<-paste(unique(hpf2[as.character(hpf2$sequenceID)==ibps,"chr.absolute"]))  
        bedstart<-as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"start.absolute"]))
        bedend<-as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"start.absolute"])) + 
          as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"LinkerStartInMolecule"])) - 1
        lnk.eim<-as.numeric(unique(hpf2[hpf2$sequenceID==ibps,"LinkerEndInMolecule"]))
        
      }#if hp
      
      tolo<-all_my_cpgis[all_my_cpgis$chrom==selchr,]# & 
                           #all_my_cpgis$chromStart >= bedstart & 
                           #all_my_cpgis$chromEnd <= bedend,]
      
      tolo$start.relative=tolo$chromStart-bedstart+1
      tolo$end.relative=tolo$chromEnd-bedstart+1
      
     # tolo<-tolo[(tolo$start.relative>=0 & tolo$end.relative<=bed.length) |
    #               (tolo$start.relative<0 & tolo$end.relative<=bed.length & tolo$end.relative>=0) |
    #               (tolo$start.relative>=0 & tolo$start.relative<=bed.length) |
    #               (tolo$end.relative>=0 & tolo$end.relative<=bed.length),]
      
      tolo<-tolo[(!as.numeric(as.character(tolo$end.relative))<0) &
                   (!as.numeric(as.character(tolo$start.relative))>=bed.length),]
      
      if((primer.type == "hp_bisulfite" | primer.type == "hp_NOME" | primer.type == "hp_genomic" | primer.type == "hp_CLEVER") &&
         hp.initial.input.type == "regions"){
        
        #trick to also show SNPs on the second part of the hairpin
        tolo.temp<-tolo
        tolo.temp$end.relative <- lnk.eim + tolo$start.relative - 1 
        tolo.temp$start.relative <- lnk.eim + tolo$end.relative - 1
        tolo<-rbind(tolo,tolo.temp)
        
      }#if hp
      
      if(nrow(tolo[as.numeric(as.character(tolo$start.relative)) < 0,])>0){
        tolo[as.numeric(as.character(tolo$start.relative)) < 0,"start.relative"]<-0
      }
      
      if(nrow(tolo[as.numeric(as.character(tolo$end.relative)) < 0,])>0){
        tolo[as.numeric(as.character(tolo$end.relative)) > bed.length,"end.relative"]<-bed.length
      }
      
      
      #write.table(tolo,file=paste(path.tracks,"tolo.object.for.plots_temp_",ibps,".txt",sep=""),
      #            sep="\t",dec=".",col.names=T,row.names=F,quote=F)
      
      srs<-as.numeric(as.character(tolo$start.relative))
      sre<-as.numeric(as.character(tolo$end.relative))
      
      if(length(srs)==0){srs<-NA}
      if(length(sre)==0){sre<-NA}
      
      lol<-rbind(lol,data.frame(relative.position=NA,
                                cg.pos=NA,
                                gc.pos=NA,
                                snp_start=NA,
                                snp_end=NA,
                                p1_start=NA,
                                p1_end=NA,
                                p2_start=NA,
                                p2_end=NA,
                                amp.start=NA,
                                amp.end=NA,
                                linker.start=NA,
                                linker.end=NA,
                                repeat.start=NA,
                                repeat.end=NA,
                                gene.start=NA,
                                gene.end=NA,
                                cpgi.start=srs,
                                cpgi.end=sre,
                                sequence.id=ibps,
                                amplicon.id="input.sequence [CpG Island]",
                                feature=NA))
    }# annotate.cpgis
    
####################################################################################################################################        
 #export(lol object)
    lol$feature<-as.character(lol$feature)
    
    #v9.7 bugfix
    lolmax<-max(lol$relative.position,na.rm=T)
    lol[lol>lolmax]<-lolmax

    write.table(lol,file=paste(path.temp,"lol.object.for.plots_",ibps,".txt",sep=""),
                sep="\t",dec=".",col.names=T,row.names=F,quote=F)
            
##############PREPARE Data for the plot
    
    if(primer.type!="NOME" & primer.type!="hp_NOME"){
      lol<-lol[lol$amplicon.id != "input.sequence [GpC]",]
    }
    
    if(check4snps == FALSE){
      lol<-lol[lol$amplicon.id != "input.sequence [SNP]",]
    }
    
    if(check4repeats == FALSE){
      lol<-lol[lol$amplicon.id != "input.sequence [Repeats]",]
    } 
    
    if(annotate.genes == FALSE){
      lol<-lol[lol$amplicon.id != "input.sequence [Genes]",]
    }
    
    if(annotate.cpg.islands ==FALSE){
      lol<-lol[lol$amplicon.id != "input.sequence [CpG Island]",]
    }

    if(!primer.type %in% c("hp_bisulfite","hp_NOME","hp_genomic","hp_CLEVER")){
      lol<-lol[lol$amplicon.id != "input.sequence [linker]",]
    }
 
##############CREATE THE PLOT 

    lopl<-ggplot(lol,aes(relative.position,amplicon.id))+
             geom_line(colour="black",size=0.7)
    
    if(check4repeats){
       if(nrow(lol[!is.na(lol$repeat.start) & !is.na(lol$repeat.end),])>0){
        lopl<-lopl+geom_segment(colour="black",size=3,
                                aes(x=repeat.start,xend=repeat.end,y=amplicon.id,yend=amplicon.id))#repeats
      }
    }
    
    if(annotate.genes){
      if(nrow(lol[!is.na(lol$gene.start) & !is.na(lol$gene.end),])>0){
        lopl<-lopl+geom_segment(colour="orange",size=3,
                                aes(x=gene.start,xend=gene.end,y=amplicon.id,yend=amplicon.id))#genes
      }
    }
    
    if(annotate.cpg.islands){
      if(nrow(lol[!is.na(lol$cpgi.start) & !is.na(lol$cpgi.end),])>0){
        lopl<-lopl+geom_segment(colour="darkgreen",size=3,
                                aes(x=cpgi.start,xend=cpgi.end,y=amplicon.id,yend=amplicon.id))#cpg islands
      }
    }
           
   if(nrow(lol[!is.na(lol$amp.start) & !is.na(lol$amp.end),])>0){
   lopl<-lopl+geom_segment(colour="darkblue",size=2,
                aes(x=amp.start,xend=amp.end,y=amplicon.id,yend=amplicon.id))   #amplicon
   }
 
   if(nrow(lol[!is.na(lol$p1_start) & !is.na(lol$p1_end),])>0){
     lopl<-lopl+geom_segment(colour="darkred",size=2.5,
                aes(x=p1_start,xend=p1_end,y=amplicon.id,yend=amplicon.id))   #primer1
   }
 
   if(nrow(lol[!is.na(lol$p2_start) & !is.na(lol$p2_end),])>0){
     lopl<-lopl+geom_segment(colour="darkred",size=2.5,
                aes(x=p2_start,xend=p2_end,y=amplicon.id,yend=amplicon.id))   #primer2
   }
     
    if(primer.type =="hp_bisulfite" | primer.type == "hp_NOME" | primer.type == "hp_genomic" | primer.type == "hp_CLEVER"){
      
           if(nrow(lol[!is.na(lol$linker.start) & !is.na(lol$linker.end),])>0){
             lopl<-lopl+geom_segment(colour="darkgrey",size=2.5,
                        aes(x=linker.start,xend=linker.end,y=amplicon.id,yend=amplicon.id))#linker
           }
    }
           
    if(nrow(lol[!is.na(lol$cg.pos),])>0){
      lopl<-lopl+geom_point(shape = 23, colour = "black", fill = "white", size = 4, stroke = 0.5,
                      aes(x = cg.pos,y=amplicon.id))#CpGs
    } 
   
   if(primer.type=="NOME" | primer.type=="hp_NOME"){
   
     if(nrow(lol[!is.na(lol$gc.pos),])>0){
       lopl<-lopl+geom_point(shape = 23, colour = "white", fill = "black", size = 4, stroke = 0.5,
                       aes(x = gc.pos,y=amplicon.id))#GpCs
     }
     
    }
   
  if(check4snps & exists("all_my_snps")){  
   if(nrow(lol[!is.na(lol$snp_start) & !is.na(lol$snp_end),])>0){
     lopl<-lopl+geom_segment(colour="black",size=6,
                             aes(x=snp_start,xend=snp_end,y=amplicon.id,yend=amplicon.id))#SNPs
   }
  }

      lopl<-lopl+facet_wrap(~sequence.id)+#separate plots by ampliconID  
           
      theme(axis.text.x = element_text(colour = "black",angle=0,size=10,hjust=0))+
      theme(axis.text.y = element_text(colour = "darkblue",angle=0,size=7,hjust=0))+
      theme(panel.grid.minor=element_blank())+
      theme(panel.background=element_rect(fill="grey90"))+
      theme(panel.background=element_rect(linetype="solid"))+
      labs(title=paste(analysis.id))
     
     ggsave(filename=paste("amplicon.design_",ibps,".png",sep=""),
            plot = lopl,path = path.graphics,width = 4,height = 4)
    }#for 
    
  }  
  log("Done.")
##########add amplicon.design plots to html report
  # log("Add graphics to html report...")
  
    # fnhm<-list.files(path.graphics)[grepl("amplicon.design",list.files(path.graphics))]
    # fnhm<-fnhm[order(fnhm)]
    
    # #for (ala in AllAmps){
    # fnir<-fnhm#[grepl(ala,fnhm)]
    
    # html.slice.fnir<-paste("<img src='",path.graphics,fnir,"' border=",2,
                           # " width=",500," height=",500,">",collapse="")
    # html.slice.fnir<-gsub("/ ","/",html.slice.fnir)
    # html.caption.fnir<-paste(fnir,collapse = "     ")
    
    # write.text(text = paste("<p align=center>",html.slice.fnir,
                            # "<br><i class=caption> ",html.caption.fnir,
                            # " </i>",collapse=""),
               # file = paste(path.html,"summary_",analysis.id,".html",sep=""),add = TRUE) 
    # write.text(text = paste("<hr  width=100% size=2>"),
               # file = paste(path.html,"summary_",analysis.id,".html",sep=""),add = TRUE) 
    
    
    # #add barplots plots to html report
    
    # fnhm<-list.files(path.graphics)[grepl("bar.",list.files(path.graphics))]
    # fnhm<-fnhm[order(fnhm)]
    
    # #for (ala in AllAmps){
    # fnir<-fnhm#[grepl(ala,fnhm)]
    
    # html.slice.fnir<-paste("<img src='",path.graphics,fnir,"' border=",2,
                           # " width=",500," height=",500,">",collapse="")
    # html.slice.fnir<-gsub("/ ","/",html.slice.fnir)
    # html.caption.fnir<-paste(fnir,collapse = "     ")
    
    # write.text(text = paste("<p align=center>",html.slice.fnir,
                            # "<br><i class=caption> ",html.caption.fnir,
                            # " </i>",collapse=""),
               # file = paste(path.html,"summary_",analysis.id,".html",sep=""),add = TRUE) 
    # write.text(text = paste("<hr  width=100% size=2>"),
               # file = paste(path.html,"summary_",analysis.id,".html",sep=""),add = TRUE)    
    # log("Done.")
}# if create.graphics
###################################################################################
###################################################################################
#return(black)
log("Primer design pipeline finished.")
}#if nrow(results2)>0

  if(nrow(results2)<1){
    black<-bed
    write.table(black,file=paste(path.wd,"primer_",analysis.id,"_blacklist.txt",sep=""),
                col.names=TRUE,row.names=FALSE,sep="\t",dec=".",quote=F)     
  }
  if(nrow(results2)>0){
    if(nrow(results3)<1){
    black<-bed
    write.table(black,file=paste(path.wd,"primer_",analysis.id,"_blacklist.txt",sep=""),
                col.names=TRUE,row.names=FALSE,sep="\t",dec=".",quote=F)  
    }
  }
#

log("reset working directory...")
setwd(start.wd)
log("done.")

log("end.")

return(black)#output for primerwheel

}#if(exists("results))

}#end

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#add functions

#######################################################################################
#function to retrieve genomic sequences from ucsc
getDNA<-function(chr,start,end,assembly="hg19",sleep=0.5){
  
  #url example for UCSC sequence retireval: http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr1:100000,200000
  chr<-gsub("chr","",chr)
  url.full<-paste("http://genome.ucsc.edu/cgi-bin/das/",assembly,"/dna?segment=chr",chr,":",formatC(start,format="f",digits=0),",",formatC(end,format="f",digits=0),sep="")
  rl<-readLines(url.full)
  rlf<-paste(rl[!grepl("<",rl)],collapse = "")
  Sys.sleep(time = sleep)
  return(rlf)
  
}

###########################################################################################################
###########################################################################################################
###########################################################################################################

##########function to do bisulfite conversion ("C" to "T", except "CG"s)
bisulfite.conversion<-function(sequence,strand="top"){#strand="top" or "bottom"
  require(gsubfn)
  #if("pathScripts" %in% ls(name=".GlobalEnv")){
  #  source(paste(pathScripts,"_Functions\\reverse.complement_jil_v1.0.r",sep=""))
  #}
  
  x<-as.character(toupper(sequence))
  
  if(strand=="bottom"){
    x<-reverse.complement(x)
  }
  
  #bis conversion
  cg2xy<-gsubfn("CG","XY",x,ignore.case=T)
  c2t<-gsubfn("C","T",cg2xy,ignore.case=T)
  out<-gsubfn("XY","CG",c2t,ignore.case=T)  
  
  out<-tolower(out)
  return(out)
}
#####Example
##mySeq<-"AAGGCCTTCGCGCGTTGCGCGCTCCCGCGGCCTCGCCTCCAACGTTCGAATTCCCGCGGCCCGGCGCGGCC"
#myCon<-bisulfite.conversion(mySeq,strand = "bottom")
#print(paste("original: ",mySeq))
#print(paste("converted:",myCon))

###############################################################################################
#####in silico conversion of genomic DNA to NOMe DNA
################################################################################
#install and load required packages.....
#install.packages("gsubfn")
#install.packages("rtf")
###########################################################################################################
###########################################################################################################
###########################################################################################################

##########function to convert "C" to "T", except those in a "CG" or "GC"context
##########V1.1 now takes care for "hidden" CGs or GCs in line-breaks...
##########V1.2 runs either on "top" or "bottom" strand...
gcCon<-function(sequence,strand="top"){
  #require(gplots)
  require(rtf)
  require(gsubfn)
  
  #if("pathScripts" %in% ls(name=".GlobalEnv")){
  #  source(paste(pathScripts,"_Functions\\reverse.complement_jil_v1.0.r",sep=""))
  #}
  
  x<-as.character(toupper(sequence))
  
  #check for number of lines and concatenate
  if (length(x)>1){
    x<-paste(x,collapse="")
  }
  
  if(strand=="bottom"){
    x<-reverse.complement(x)
  }
  
  #do in silico conversion   
  cg2xy<-gsubfn("CG","XY",x,ignore.case=T)
  gc2yz<-gsubfn("GC","YZ",cg2xy,ignore.case=T)
  yc2yz<-gsubfn("YC","YZ",gc2yz,ignore.case=T)
  out<-gsubfn("C","T",yc2yz,ignore.case=T)
  out<-gsubfn("X","C",out,ignore.case=T)
  out<-gsubfn("Z","C",out,ignore.case=T)
  out.final<-gsubfn("Y","G",out,ignore.case=T)
  
  return(out.final)
}

######Example
#mySeq<-"AAGGCCTTCGCGCGTTGCGCGCTCCCGCGGCCTCGCCTCCAACGTTCGAATTCCCGCGGCCCGGCGCGGCC"
#myCon<-gcCon(mySeq)
#print("Example...")
#print(paste("original: ",mySeq))
#print(paste("converted:",myCon))

###########################################################################################################
###########################################################################################################
###########################################################################################################


#this function calculates primer melting temperatures in °C.

calculate.Tm<-function(sequence){
  #formula used: Tm= 62.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC) #if seq longer than 13 nucleotides
  #formula used: (wA+xT) * 2 + (yG+zC) * 4  #if seq shorter than 14 nucleotides
  seq<-toupper(sequence)
  
  nA<-nchar(x = gsub("[C|G|T]","",seq))
  nC<-nchar(x = gsub("[A|G|T]","",seq))
  nG<-nchar(x = gsub("[C|A|T]","",seq))
  nT<-nchar(x = gsub("[C|A|G]","",seq))
  
  if (nchar(seq)>13){
    tm<-62.9 + 41*(nG+nC-16.4) / (nA+nT+nG+nC)
  }
  else{
    tm<-(nA+nT) * 2 + (nG+nC) * 4 
  }
  return(tm)  
}

###########################################################################################################
###########################################################################################################
###########################################################################################################

#reverse.complement a dna sequence

strReverse <- function(x){
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}

reverse.complement<-function(sequence,reverse=TRUE,complement=TRUE){
  temp<-as.character(tolower(sequence))
  if(reverse==TRUE){
    #temp<-paste(temp[length(temp):1]) #this was a major bug!!!!
    temp<-strReverse(temp)
  }
  if(complement==TRUE){
    temp<-gsub("a","1",temp)
    temp<-gsub("c","2",temp)
    temp<-gsub("g","c",temp)
    temp<-gsub("t","a",temp)
    
    temp<-gsub("1","t",temp)
    temp<-gsub("2","g",temp)
  }
  out<-temp
  return(out)
}

#example
#reverse.complement("agattcggggc")
###########################################################################################################
###########################################################################################################
###########################################################################################################


#write a sequence in .fasta format
write.fasta<-function(sequence="Example",sequence.id="my.fasta",filename=paste(getwd(),"/my.fasta.FASTA",sep="")){
  cat(paste('>',sequence.id, '\n', sequence,sep=""),file=filename)
}

#example
#write.fasta("AGATGT","my.first.fasta",paste(getwd(),"/my.fasta.FASTA",sep=""))

###########################################################################################################
###########################################################################################################
###########################################################################################################


#read.fasta file
read.my.fasta<-function(filename){
  require("seqinr")
  fasta<-read.fasta(file = filename, 
                    seqtype = "DNA", as.string = TRUE, forceDNAtolower = TRUE,
                    set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = FALSE)
  return(fasta)  
}

###########################################################################################################
###########################################################################################################
###########################################################################################################

#Export some text to a file...
#adds also date info to the specified text, useful for logging information...
#Warning: add=FALSE: If file already exists it will be overwritten!
#Optional: add=TRUE, opens an existing file and adds the text without overwriting original contents.

write.log<-function(text="BLAH!",file=paste(getwd(),"logfile.txt",sep=""),add=TRUE){
  fileConn<-file(file)
  if(!add){
    writeLines(paste(date(),".....",text,sep=""),fileConn)
  }
  
  if(add){
    temp<-readLines(fileConn)
    writeLines(c(temp,text),fileConn)
  }
  close(fileConn)
  print(paste(date(),".....",text,sep=""))
}

###########################################################################################################
###########################################################################################################
###########################################################################################################


#Export some text to a file...
#Warning: If file already exists it will be overwritten!
#Optional: add=TRUE, opens an existing file and adds the text without overwriting original contents.

write.text<-function(text="BLAH!",file="text.txt",add=FALSE){
  cat(text,file=file,append=add)
}

###########################################################################################################
###########################################################################################################
###########################################################################################################


#exports a data.frame as an HTML table
#####INPUT:
# - data . the data.frame that will be exported.
# - filename: a filename for the created table.
# - caption: a html caption for the file.
# - label: a label for the html table.
# - include.rownames/colnames: TRUE/FALSE if row-and colnames will be included, respectively.
#####OUTPUT:
# a html file.

write.html<-function(data,filename,caption="HTML File",label="HTML File",include.rownames=TRUE,include.colnames=TRUE,digits=2){
  require(xtable)
  xHTML<-xtable(data, caption=caption, label=label, type="html",digits=digits)
  print.xtable(xHTML,type="html",file=paste(filename,".html",sep=""),
               append=FALSE, include.rownames=include.rownames, include.colnames=include.colnames) 
}

###########################################################################################################
###########################################################################################################
###########################################################################################################


# #creates a html main page that links to other html pages. Labels arguments holds the description of the html that  is linked to...
# html.main.page<-function(fn.links,
                         # labels,
                         # filename=paste(getwd(),"MainPage.html",sep="")){
  
  # #nested function html.report
  
  # html.report<-function(filename.out="html.report.html",filenames,txt.header=TRUE,txt.sep="\t"){
    
    # require(R2HTML)
  # require(tools)
    
    # htmlWidth<-600
    # htmlHeight<-600
    # GraphBorderSize<-5
    
    # graphics<-c("png","jpeg","svg","jpg")
    # txt<-c("txt","csv","tsv","fasta","bed")
    # graphics<-c(graphics,toupper(graphics))
    # txt<-c(txt,toupper(txt))
    
    # for (i in 1:length(filenames)){
      # temp<-file_path_as_absolute(filenames[i])
      # splits<-strsplit(temp,"/")
      # fn<-splits[[1]][length(splits[[1]])]
      # file.type<-file_ext(fn)
      # temp.name<-gsub(file.type,"",fn)
      
      # if(file.type %in% graphics){
        # HTMLInsertGraph(Caption=paste(temp.name),
                        # file=filename.out,append=TRUE,
                        # GraphFileName=temp,
                        # GraphBorder=GraphBorderSize,
                        # Width=htmlWidth,Height=htmlHeight)
        # HTMLhr(file=filename.out,append=TRUE)
      # }#if graphix
      
      # if(file.type %in% txt){
        # df<-read.table(filenames[i],header=txt.header,sep=txt.sep)
        # if(nrow(df)==0 | ncol(df)<2){
          # df<-readLines(filenames[i])
        # }
        # HTML(x=df,file=filename.out,
             # classcaption="captiondataframe",
             # caption=paste(temp.name),
             # Border=1,innerBorder=1,append=TRUE,align="center",
             # captionalign="top",row.names=FALSE)
        # HTMLhr(file=filename.out,append=TRUE)
      # }#if txt
      
    # }#for i
  # }
  
  # #nested function html.report
  
  # out<-data.frame(Contents=labels,
                  # Links=print(paste("<a href= ",fn.links," target='_blank'",'>',"link","</a>",sep="")))
  
  # write.table(out,file=paste(filename,".txt",sep=""),sep="\t",dec=".",col.names=TRUE,row.names=FALSE)
  
  # html.report(filename.out = filename,filenames = paste(filename,".txt",sep=""))
# }

###########################################################################################################
###########################################################################################################
###########################################################################################################

##############################################################################################
# this function expects a .bed file with the following set-up:
# - tab.delimited & Header
# "chr":  - columna specifying the chromosome
#"start"&"end" - columns specifying start/end of the sequence to be retrieved
#"assembly" - the genome assembly ("hg18" / "hg19" / "mm9" /"mm10")
#"species" - "the species to be used (currently only "Hsapiens" and "Mmusculus").

#####Output:
# - returns a data.frame containing the requested information (dna sequence & length) 
# - writes the data into a .txt delimited file,

get.sequence.info<-function(filename,filename.out){
  #if("pathScripts" %in% ls(name=".GlobalEnv")){
  #  source(paste(pathScripts,"_Functions/get.dna_jil_v1.0.r",sep=""))
  #}
  
  bed<-read.table(file=filename,header=TRUE,sep="\t",dec=".")
  bed$sequenceID<-as.character(bed$sequenceID)
  id.na<-is.na(bed$sequenceID)|bed$sequenceID==""
  bed[id.na,"sequenceID"]<-paste(bed[id.na,"chr"],"_",bed[id.na,"start"],"_",bed[id.na,"end"],sep="")
  bed$sequence.length<-bed$end-bed$start+1
  bed$sequence.adress<-paste(bed[,"chr"],"_",bed[,"start"],"_",bed[,"end"],sep="")
  bed$sequence<-NA
  
  for (i1 in 1:nrow(bed)){
    i1spec<-as.character(bed[i1,"species"])
    i1ass<-as.character(bed[i1,"assembly"])
    i1chr<-as.character(bed[i1,"chr"])
    i1start<-bed[i1,"start"]
    i1end<-bed[i1,"end"]
    
    if(i1spec=="Hsapiens" & i1ass=="hg19"){
      bed[i1,"sequence"]<-get.dna.human.hg19(chr=i1chr,start=i1start,end=i1end)
      next}
    if(i1spec=="Hsapiens" & i1ass=="hg18"){
      bed[i1,"sequence"]<-get.dna.human.hg18(chr=i1chr,start=i1start,end=i1end)
      next
    }
    if(i1spec=="Mmusculus" & i1ass=="mm9"){
      bed[i1,"sequence"]<-get.dna.mouse.mm9(chr=i1chr,start=i1start,end=i1end)
      next
    }
    if(i1spec=="Mmusculus" & i1ass=="mm10"){
      bed[i1,"sequence"]<-get.dna.mouse.mm10(chr=i1chr,start=i1start,end=i1end)
      next
    }
    else{
      print("This genome assembly is currently not supported.")
      bed[i1,"sequence"]<-NA     
      next  
    }    
  } 
  
  write.table(bed,file=filename.out,sep="\t",dec=".",col.names=TRUE,row.names=FALSE)
  return(bed)
}

###########################################################################################################
###########################################################################################################
###########################################################################################################
sequence.characterization<-function(sequence,
                                    feature1="cg",
                                    feature1.subpos=1,#1st position = 1, 2nd position =2,...
                                    feature2="gc",
                                    feature2.subpos=2){#1st position = 1, 2nd position =2,...
  myseq<-tolower(sequence)
  feature1<-tolower(feature1)
  feature2<-tolower(feature2)
  
  ##########################################################################################################
  
  #all pos splits
  s0<-strsplit(x = myseq,split = "")
  
  
  ###########
  #all feature1
  s1<-gregexpr(pattern = feature1,text = myseq)
  
  if(length(s1[[1]])>1 || length(s1[[1]])==1 & s1[[1]][1]>0 ){
    pos1<-as.numeric(s1[[1]])
    subpos1<-pos1+feature1.subpos-1
    n.f1<-length(subpos1)
    names(subpos1)<-paste(feature1,"_",1:n.f1,sep="")
  }
  
  if(length(s1[[1]])==1 & s1[[1]][1]<0 ){
    pos1<-0
    subpos1<-0
    n.f1<-0
  }
  
  ###########
  #all feature2
  s2<-gregexpr(pattern = feature2,text = myseq)
  if(length(s2[[1]])>1 || length(s2[[1]])==1 & s2[[1]][1]>0 ){
    pos2<-as.numeric(s2[[1]])
    subpos2<-pos2+feature2.subpos-1
    n.f2<-length(subpos2)
    
    names(subpos2)<-paste(feature2,"_",1:n.f2,sep="")}
  if(length(s2[[1]])==1 & s2[[1]][1]<0 ){
    pos2<-0
    subpos2<-0
    n.f2<-0
  }
  
  ###########
  #feature.overlapps
  ol1 <- NA
  ol2 <- NA
  
  if(n.f1>0){
    if(n.f2>0){
      ol1<-subpos1[subpos1 %in% subpos2]
      ol2<-subpos2[subpos2 %in% subpos1] 
    }
    
    if(n.f2==0){
      ol1<-NA
      ol2<-NA
    }
  }
  
  if(n.f2>0){
    if(n.f1==0){
      ol1<-NA
      ol2<-NA
    }
  }
  
  ###########
  #results as data.frame
  df<-data.frame(basecount=1:nchar(myseq),base=unlist(s0))
  
  ua<-unlist(s0)=="a"
  uc<-unlist(s0)=="c"
  ug<-unlist(s0)=="g"
  ut<-unlist(s0)=="t"
  un<-unlist(s0)=="n"
  
  df$a<-ua
  df$c<-uc
  df$g<-ug
  df$t<-ut
  df$n<-un
  
  al<-length(ua[ua==TRUE])
  cl<-length(uc[uc==TRUE])
  gl<-length(ug[ug==TRUE])
  tl<-length(ut[ut==TRUE])
  nl<-length(un[un==TRUE])
  
  df$a_count<-NA
  df$c_count<-NA
  df$g_count<-NA
  df$t_count<-NA
  df$n_count<-NA
  
  df$a_count[df$a==TRUE]<-paste("a_",1:al,sep="")
  df$c_count[df$c==TRUE]<-paste("c_",1:cl,sep="")
  df$g_count[df$g==TRUE]<-paste("g_",1:gl,sep="")
  df$t_count[df$t==TRUE]<-paste("t_",1:tl,sep="")
  df$n_count[df$n==TRUE]<-paste("n_",1:nl,sep="")
  
  df[,paste(feature1)]<-FALSE
  df[subpos1,paste(feature1)]<-TRUE
  df[,paste(feature2)]<-FALSE
  df[subpos2,paste(feature2)]<-TRUE
  df[,paste(feature1,"_overlapps_",feature2,sep="")]<-FALSE
  df[,paste(feature2,"_overlapps_",feature1,sep="")]<-FALSE
  
  ###########
  if(!is.na(ol1[1])){
    df[ol1,paste(feature1,"_overlapps_",feature2,sep="")]<-TRUE
  }
  if(!is.na(ol2[1])){
    df[ol2,paste(feature2,"_overlapps_",feature1,sep="")]<-TRUE
  }
  
  ###########
  df[,paste(feature1,"_count",sep="")]<-NA
  df[subpos1,paste(feature1,"_count",sep="")]<-names(subpos1)
  df[,paste(feature2,"_count",sep="")]<-NA
  df[subpos2,paste(feature2,"_count",sep="")]<-names(subpos2)
  
  df[is.na(df)]<-""
  
  return(df)
}

##################################################################################################################
##################################################################################################################
##################################################################################################################
#hairpinizer: code for processing sequence suitable for hairpin analysis
#v2.0 now support of degenerated restriction sites
hairpinizer<-function(seq,inputID="sequence1",REseq="CCGG",RE.id="MspI",
                      linker.seq="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",length.hairpin.region=300){
  #if("pathScripts" %in% ls(name=".GlobalEnv")){
  #  source(paste(pathScripts,"_Functions/reverse.complement_jil_v1.0.r",sep=""))
  #}
  
  seq<-tolower(seq)
  seq.length<-nchar(seq)
  REseq<-tolower(REseq)
  REseq.length<-nchar(REseq)
  #Analyze REseq for degeneration...
  REpatterns<-gsub("n","[a|c|g|t]",gsub("r","[a|g]",gsub("y","[c|t]",gsub("w","[a|t]",gsub("s","[c|g]",
                                                                                           gsub("m","[a|c]",gsub("k","[t|g]",gsub("b","[g|c|t]",gsub("h","[a|c|t]",gsub("d","[a|g|t]",gsub("v","[a|g|c]",REseq)))))))))))
  
  out<-data.frame(inputID=inputID,sequenceID=NA,RE=RE.id,REseq=REseq,REpatterns=REpatterns,hp_Region_Length=NA,
                  sequence=NA,hp_Region_Start_Relative=NA,
                  hp_Region_End_Relative=NA,hp_Sequence=NA,
                  hp_Sequence_ReverseComplement=NA,linker_Sequence=NA,Comment=NA)
  
  
  
  if(length.hairpin.region>seq.length){
    print(paste(date(),".....Hairpin Region is larger than input sequence",sep=""))
  }
  
  else{
    print(paste(date(),".....Check for presence of restriction sites...",sep=""))
    re.sites<-gregexpr(REpatterns, seq)#[[1]]
    
    if(re.sites[[1]][1]== -1){
      print(paste(date(),".....No restriction sites '",REseq,"' present in sequence",sep=""))
      out[1,"Comment"]<-paste("No restriction sites '",REseq,"' present in sequence",sep="")
    }
    
    else{
      re.hits<-as.numeric(re.sites[[1]])
      n.re.hits<-length(re.hits)
      print(paste(date(),".....Identified ",n.re.hits," Restriction sites for '",REseq,"' in sequence",sep=""))
      
      for (i in 1:n.re.hits){
        print(paste(date(),".....Analyze site ",i,"...",sep=""))
        re.site.start<-re.hits[i]
        hairpin.start<-re.site.start-length.hairpin.region
        re.site.end<-re.site.start+REseq.length-1
        if(i==1){
          if(hairpin.start<1){
            hairpin.start<-1
            out[i,"Comment"]<-"RE site is close to sequence start...trimmed."
          }
        }
        if(i>1){
          pre.site.end<-re.hits[i-1]+REseq.length-1
          dist.pre.site<-re.site.start-pre.site.end
          if(dist.pre.site<length.hairpin.region){
            hairpin.start<-pre.site.end+1
            out[i,"Comment"]<-"Upstream RE site is closer than desired hp length...trimmed."
          }
        }
        out[i,"hp_Region_Start_Relative"]<-hairpin.start
        out[i,"hp_Region_End_Relative"]<-re.site.end
        out[i,"hp_Region_Length"]<-re.site.end-hairpin.start+1
        temp.seq<-substr(seq,hairpin.start,re.site.end)  
        out[i,"hp_Sequence"]<-temp.seq
        out[i,"hp_Sequence_ReverseComplement"]<-reverse.complement(temp.seq)
        out[i,"sequence"]<-paste(temp.seq,linker.seq,reverse.complement(sequence = temp.seq),sep="")
        out[i,"inputID"]<-inputID
        out[i,"RE"]<-RE.id
        out[i,"REseq"]<-REseq
        out[i,"REpatterns"]<-REpatterns
        out[i,"linker_Sequence"]<-linker.seq
        out[i,"sequenceID"]<-paste(inputID,"_",out[i,"RE"],"_",out[i,"hp_Region_Start_Relative"],"_",
                                   out[i,"hp_Region_End_Relative"],sep="")
      }#for
    }#else1
  }#else2
  
  out$linker_Sequence<-linker.seq
  out$RE<-RE.id
  out$REseq<-REseq
  out$REpatterns<-REpatterns
  print(paste(date(),".....Done.",sep=""))
  return(out)
}#function

#####EXAMPLE. Not Run
#seq<-"ctgggacggcgaggacgcgaacctgcggagccgcgggtgtctgaggagccgcagggaacccccggcctagccgcgcccgcgtgtggccggagctgcgggccgggactgtgtccaggacagagccacaagcttgtccccagctcagggaggtccaggggcggcagagggagcgacaggctgcgaagcccaccggtgaccacgtgtgaacccgcgtgcgcccccagctcggccactccgtgcgggtctgccctcaccgcagctccggcctgccggccctgcctgctcccgtggtctgggatgtggccccggtgaggacccggccccatcaggcacagggtggatgtctgtggagtgaggtgtgtgtgacatattcatgtgaccacccgtgcagcgtcacgcgcctggccctgccgatgacaagggtgtgggcctgcgtgggcatgactgtgtgtgtgacacagagtgatgttgctgtgacccgtggctgcactccccacatcaccggctttcacagccttccggtaaagtgctgtgttctcccttctgtgtcttcgctgggacctggggcaagggtgggtgtggcccccacagctggagtcagcttctgtggggccttcccgagccctccccaccctggaccagaggcccagctggttggagcaggaagtacctgggctctggggtcagggatgggaaggctgaggaggcctgcgtgagctggacctggcctgggccctcctggccgtgcctgcctggtggtgcaggattcctggggctgatgacagacggggtagggctggggttggcgagcctcctgccgatacctcacgtagctgacctctgactcttccccagccaggctggccctgggagttgccggagagtcagtggatctgcaggctgcacgctggctgttacctttgcttctgggttcccacaggggtcatggttctgtggttctccagtcagggaccctagcagggccatggggcgtgacttcctggaggtgtggcctagtatggccacggcagaggatgggggaagagaaaggccccctttgtcagcccccgggctctgaaccaagctgaagccctccccctg"
#hp<-hairpinizer(seq,REseq = "AGCT",RE.id = "AluI",length.hairpin.region = 50) 
##########################################################################################################################################
#a table with some RE ids and sites...
#REset<-data.frame(RE.id=c("MspI","AluI","TaqI"),
#                  RE.seq=c("CCGG","AGCT","TCGA"))
##################################################################################################################################
#use hairpinizer to create hp molecules for several different REestriction enzymes.
#for (i in 1:nrow(REset)){
#  if (i==1){
#   hp<-hairpinizer(seq = seq,REseq = REset[i,"RE.seq"],RE.id = REset[i,"RE.id"],length.hairpin.region = 300) 
#  }
#  if(i>1){
#  temp<-hairpinizer(seq = seq,REseq = REset[i,"RE.seq"],RE.id = REset[i,"RE.id"],length.hairpin.region = 300) 
#  hp<-rbind(hp,temp)  
#  }
#}

#########################################################################################################
#########################################################################################################
#########################################################################################################
#BISULFITE PRIMER DESIGN

bisulfite.primer.design<-function(sequence,
                                  sequence.id,
                                  min.length.primer=26, 
                                  max.length.primer=32, 
                                  min.Tm.primer=52,
                                  max.Tm.primer=60, 
                                  max.Tm.difference.primer=2, 
                                  low.complexity.primer.removal=TRUE, 
                                  max.bins.low.complexity=7,
                                  remove.primers.with.n=TRUE,
                                  min.C2T.primer1=3,
                                  min.G2A.primer2=3,
                                  min.length.amplicon=200,
                                  max.length.amplicon=500, 
                                  min.number.gc.amplicon=0, 
                                  min.number.cg.amplicon=7,
                                  primer.align.binsize=8,
                                  strand="top",
                                  mode="exact"){
  
  version.id<-2.0
  
  ########################################################################
  
  #Load helper scripts
  #if("pathScripts" %in% ls(name=".GlobalEnv")){
  #  source(paste(pathScripts,"_Functions\\calculate.Tm_jil_v1.0.r",sep=""))
  #  source(paste(pathScripts,"_Functions\\read.my.fasta_jil_v1.0.r",sep=""))
  #  source(paste(pathScripts,"_Functions\\reverse.complement_jil_v1.0.r",sep=""))
  #  source(paste(pathScripts,"_Functions\\bisulfite.conversion_jil_V1.0.r",sep=""))
  #}
  
  #######################################################################
  
  #initiate sequence
  seq<-as.character(sequence)
  seq.id<-as.character(sequence.id)
  print(paste("Start Primer Design for Sequence: ",seq.id," on ",strand," strand...[",date(),"]",sep=""))
  print(paste("Length of Sequence: ",nchar(seq),sep=""))
  
  ########################################################################
  
  #perform bisulfite conversion of genomic sequence
  nome<-bisulfite.conversion(seq,strand = strand)
  nome<-tolower(nome)
  
  #########################################################################
  
  #split the nome sequence by all "c"s.
  c.split<-unlist(strsplit(x=nome,split="c"))
  names(c.split)<-paste("fragment",1:length(c.split),sep="")
  n.fragments<-length(c.split)
  
  ########################################################################
  
  #Remove splitted sequences that are shorter than the minimum length of the desired oligos.
  c.split.length<-lapply(c.split,nchar)
  long.enough<-c.split[c.split.length>=min.length.primer]
  
  ########################################################################
  
  #Prepare NEGATIVE RESULT OUTPUT
  
  out<-data.frame(sequence.id=NA,amplicon.id=NA,amplicon.length=NA,amplicon.start.relative=NA,amplicon.end.relative=NA,
                  dna.strand=NA,nGCs=NA,nCGs=NA,primer1.sequence=NA,primer2.sequence=NA, 
                  primer1.length=NA,primer2.length=NA,primer1.tm=NA,primer2.tm=NA,tm.difference=NA,
                  primer1.start.relative=NA,primer1.end.relative=NA,primer2.start.relative=NA,
                  primer2.end.relative=NA,primer.pair.id=NA,primer1.id=NA,primer2.id=NA,fragment1.id=NA,
                  fragment2.id=NA,fragment12.id=NA,fragment12.primer.ids=NA,primer1.gc.content=NA,
                  primer2.gc.content=NA,amplicon.gc.content=NA,amplicon.genomic.gc.content=NA,
                  primer1.c2t.conversion=NA,primer2.g2a.conversion=NA,primer1.self.alignment=NA,
                  primer2.self.alignment=NA,primer1.primer2.alignment=NA,amplicon.sequence=NA,
                  primer1.sequence.genomic=NA,primer2.sequence.genomic=NA,amplicon.sequence.genomic=NA,
                  index=NA)
  
  #######################################################################
  
  if (length(long.enough)<2){
    print("No Primer Design possible. Try to reduce min.length.primer.")
    print(paste(date()))
    return(out)
  }
  
  #######################################################################
  
  #CalcUlate the Tm for the remaining ones...
  tms<-lapply(long.enough,calculate.Tm)
  #Remove those oligos with too low Tm...
  tm.good<-as.list(long.enough[tms>min.Tm.primer])
  
  if (length(tm.good)<2){
    print("No Primer Design possible. Try to lower min.Tm.primer.")
    print(paste(date()))
    return(out)
  }
  
  ########################################################################
  
  #establish primers and organize them
  print("Generate & analyze individual primer sequences for each subfragment...")
  
  tmg<-data.frame(fragment.sequence=as.character(tm.good),fragment.id=names(tm.good))
  tmg$fragment.sequence<-as.character(tmg$fragment.sequence)
  tmg$fragment.length<-nchar(tmg$fragment.sequence)
  tmg$primer.start.starts<-1
  tmg$primer.start.ends<-tmg$fragment.length-min.length.primer
  tmg<-tmg[tmg$primer.start.ends>0,]
  tmg$primer.end.starts<-1+min.length.primer-1
  
  subfrags<-names(tm.good)[names(tm.good) %in% tmg$fragment.id]
  
  #very powerfull way! no loops here!
  #it is much quicker to use this way though it creates some wrong primers.
  #these can be removed in a later step easily.
  
  #primer length ranges to build start and end points for subfragments...
  plr<-max.length.primer-min.length.primer
  pr.starts<-1: (max(tmg$fragment.length)-min.length.primer+1)
  
  #
  if(mode=="exact"){
    st<-rep(pr.starts,each=plr+1)
    pr.dist<-min.length.primer:max.length.primer
    distrep<- rep(pr.dist,plr+1)
    en<- st+distrep-1
  }
  
  if(mode=="fast"){
    st<-rep(pr.starts,each=3)
    pr.dist<-c(min.length.primer,ceiling(mean(c(min.length.primer,max.length.primer),na.rm=TRUE)),max.length.primer)
    distrep<- rep(pr.dist,plr+1)
    en<- st+distrep-1
  }
  
  ############very powerfull!!!!
  pr.seqs<-sapply(X = tmg$fragment.sequence, function(x){
    prs<-as.data.frame(mapply(FUN=substr,start=st, stop=en, 
                              MoreArgs=list(x=x)))
  })
  ############
  pr.seqs<-sapply(pr.seqs,as.character)
  colnames(pr.seqs)<-subfrags
  pr.char<-sapply(pr.seqs,nchar)
  pr.seqs[pr.char<min.length.primer]<-"z"
  
  #convert dataframe to character for further handling. also assign primer.names.
  pr.sq<-as.character(pr.seqs)
  names(pr.sq)<-paste(rep(subfrags,each=nrow(pr.seqs)),"___",1:nrow(pr.seqs),sep="")
  
  #remove 'artificial' primers by deleting primers with length<min.length.primer
  prl<-sapply(pr.sq,nchar)
  pr.sq<-pr.sq[prl >= min.length.primer]
  print("Done.")
  
  ##################################################################
  
  #calculate primer melting temps.
  print ("Calculate primer melting temperatures...")
  pr.tm<-sapply(X = pr.sq,calculate.Tm)
  print("Done.")
  
  #################################################################
  
  #remove primers with inappropriate melting temps.
  print ("Remove primers with Tm too high or low...")
  f.good.tm<- pr.tm >= min.Tm.primer & pr.tm <= max.Tm.primer
  pr.sq<-pr.sq[f.good.tm]
  pr.tm<-pr.tm[f.good.tm]
  print("Done.")
  
  ########################################################################
  
  print ("Remove primer combinations that do not cover enough sites of interest...")
  #establish primer combinations that cover enough sites of interest
  pr.pairs<-expand.grid(primer2.sequence=pr.sq,primer1.sequence=pr.sq)
  pr.pairs<-cbind(pr.pairs,expand.grid(primer2.tm=pr.tm,primer1.tm=pr.tm),
                  expand.grid(primer2.id=names(pr.tm),primer1.id=names(pr.tm)))
  print (paste("Done. Found ",nrow(pr.pairs)," primer combinations. Analyze them...",sep=""))
  
  pr.pairs$primer1.sequence<-as.character(pr.pairs$primer1.sequence)
  pr.pairs$primer2.sequence<-as.character(pr.pairs$primer2.sequence)
  pr.pairs$primer1.tm<-as.numeric(pr.pairs$primer1.tm)
  pr.pairs$primer2.tm<-as.numeric(pr.pairs$primer2.tm)
  
  #Some Basic calculations for all potential primers
  pr.pairs$primer1.length<-nchar(pr.pairs$primer1.sequence)
  pr.pairs$primer2.length<-nchar(pr.pairs$primer2.sequence)
  pr.pairs$tm.difference<-abs(pr.pairs$primer1.tm - pr.pairs$primer2.tm)
  pr.pairs$fragment1.id<-gsub("___.*","",pr.pairs$primer1.id)
  pr.pairs$fragment2.id<-gsub("___.*","",pr.pairs$primer2.id)
  
  #######
  
  # keep only if enough features are present between the two subfragments....
  print ("Remove primer combinations with too few sites of interest...")
  pr.pairs<-pr.pairs[as.numeric(gsub("fragment","",pr.pairs$fragment2.id))-
                       as.numeric(gsub("fragment","",pr.pairs$fragment1.id))>= min.number.gc.amplicon+min.number.cg.amplicon,]
  
  print("Done.")
  
  if(nrow(pr.pairs)<1){
    print("Could not identify a primer pair that includes sites of interest.")
    print(paste(date()))
    return(out)
  }
  
  ######################################################################
  
  #remove primer combinations that do not fullfil tm.delta criteria
  print ("Remove primer combinations with too high Tm difference...")
  pr.pairs2<-pr.pairs[pr.pairs$tm.difference <= max.Tm.difference.primer,]
  print("Done.")
  
  if(nrow(pr.pairs2)<1){
    print("Could not identify a primer pair that includes sites of interest.")
    print(paste(date()))
    return(out)
  }
  
  ######################################################################
  
  #Remove primers that have "N"s
  if(remove.primers.with.n){
    print("Remove.primers with 'N's")
    p1.n<-grepl("N",pr.pairs2[,"primer1.sequence"],ignore.case=TRUE)
    p2.n<-grepl("N",pr.pairs2[,"primer2.sequence"],ignore.case=TRUE)
    pr.pairs2<-pr.pairs2[!p1.n & !p2.n,]
  }
  
  if(nrow(pr.pairs2)<1){
    print("After removal of 'N' containing primers nothing is left.")
    print(paste(date()))
    return(out)
  }
  
  #####################################################################
  
  #Match relative position of individual primers
  print("Calculate relative primer Positions...")
  len.1<-unlist(sapply(X = pr.pairs2$primer1.sequence,regexpr,nome,useBytes=TRUE,fixed=TRUE))
  len.2<-unlist(sapply(X = pr.pairs2$primer2.sequence,regexpr,nome,useBytes=TRUE,fixed=TRUE))
  
  len.2e<-len.2+nchar(names(len.2))-1
  len<-len.2e-len.1 #v1.7: ~7 sec
  
  pr.pairs2$amplicon.length<-as.numeric(len)
  pr.pairs2$amplicon.start.relative<-len.1
  pr.pairs2$amplicon.end.relative<-len.2e
  print("Done.")
  
  #########################################################
  
  #Select amplicons by AmpliconLength
  print("Filter amplicons by min/max length...")
  amp.sel<-pr.pairs2[pr.pairs2$amplicon.length>=min.length.amplicon & 
                       pr.pairs2$amplicon.length<=max.length.amplicon,]
  
  if (nrow(amp.sel)==0){
    print("No primers found. Try to change parameters [Amplicon Length].")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel), " primer pairs survived...",sep=""))
  
  ########################################################
  
  #fetch amplicon sequences for each primer combination.
  amp.sel$amplicon.sequence<-mapply(substr,nome,
                                    start = amp.sel$amplicon.start.relative,
                                    stop = amp.sel$amplicon.end.relative)
  amp.sel2 <- amp.sel #for historic reason, actually nothing happens here.
  amp.sel2$fragment12.id<-paste(amp.sel2$fragment1.id,"_",amp.sel2$fragment2.id,sep="")
  amp.sel2$fragment12.primer.ids<-paste(amp.sel2$primer1.id,"_",amp.sel2$primer2.id,sep="")
  
  #######################################################
  
  #select best tm primer combination for individual amplicons....
  print("Select optimal primer pair combination for Tm...")
  
  #VERY POWERFULL FUNCTION:
  lvl.frgs<-levels(factor(amp.sel2$fragment12.id))
  mins<-sapply(X = lvl.frgs,FUN = function(x) {y=which.min(amp.sel2[amp.sel2$fragment12.id==x,
                                                                    "tm.difference"])
  return(amp.sel2[amp.sel2$fragment12.id==x,"fragment12.primer.ids"][y])}
  )
  amp.sel2<-amp.sel2[match(mins,amp.sel2$fragment12.primer.ids),]
  
  ###########################################################   
  
  # check amplicons for gc/cg contents
  print("Scan amplicons for information content (GpCs/CpGs)...")
  
  amp.sel2$nCGs<-sapply(X=amp.sel2$amplicon.sequence,FUN = function(x){
    cgs<-as.numeric(unlist(gregexpr("cg",text = x,ignore.case=TRUE)))
    ncgs<-length(cgs)
    if(cgs[1]== -1){
      ncgs<-0
    }
    return(ncgs)
  })
  
  amp.sel2$nGCs<-sapply(X=amp.sel2$amplicon.sequence,FUN = function(x){
    gcs<-as.numeric(unlist(gregexpr("gc",text = x,ignore.case=TRUE)))
    ngcs<-length(gcs)
    if(gcs[1]== -1){
      ngcs<-0
    }
    return(ngcs)
  })
  
  #########################################################################
  
  #filter for gc and gc contents
  print("Filter amplicons for CpGs content...")
  amp.sel3<-amp.sel2[amp.sel2$nCGs>=min.number.cg.amplicon,]
  
  if(nrow(amp.sel3)==0){
    print("No Amplicon survived...Try to lower number of analyzable CpGs.")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel3)," amplicons survived...")) 
  
  #######################################################################################
  
  #check for numbers of CpG in amplicons
  print("Filter amplicons for GpC content...")
  
  amp.sel4<-amp.sel3[amp.sel3$nGCs>=min.number.gc.amplicon,]
  
  if(nrow(amp.sel4)==0){
    print("No Amplicon survived...Try to lower number of analyzable GpCs.")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel4)," amplicons survived..."))  
  
  #######################################################################################
  
  #LOW COMPLEXITY PRIMER REMOVAL
  if(low.complexity.primer.removal) {
    print("Analyse primers for low complexities...")
    
    p1<-c("a","c","g","t")
    p1m<-sapply(p1,rep,each=max.bins.low.complexity)
    p1m<-apply(p1m,2,paste,collapse="")
    p2<-as.character(sapply(p1,FUN = paste,c(p1,p1),sep=""))
    p2m<-as.character(sapply(p2,FUN = function(x) paste(rep(x,max.bins.low.complexity),collapse="")))
    p3<-as.character(sapply(p2,FUN = paste,c(p1,p1,p1),sep=""))
    p3m<-as.character(sapply(p3,FUN = function(x) paste(rep(x,max.bins.low.complexity),collapse="")))
    
    lowlist<-c(p1m,p2m,p3m)# all evil patterns that will be checked in all primers
    
    #check primer 1 & 2
    if(nrow(amp.sel4)>1){
      check1<-as.data.frame(mapply(FUN=regexpr,pattern=lowlist, 
                                   MoreArgs=list(text=amp.sel4$primer1.sequence)))
      check12<-cbind(check1,mapply(FUN=regexpr,pattern=lowlist, 
                                   MoreArgs=list(text=amp.sel4$primer2.sequence)))
      max.check<-apply(check12,1,max,na.rm=T)
      
      #remove amplicons with evil patterns in primers.
      amp.sel4<-amp.sel4[max.check<0,]
    }
    
    if(nrow(amp.sel4)==1){
      check1<-max(sapply(X = lowlist,FUN = regexpr, amp.sel4[1,"primer1.sequence"]))
      check2<-max(sapply(X = lowlist,FUN = regexpr, amp.sel4[1,"primer2.sequence"]))
      max.check<-max(c(check1,check2))
      amp.sel4<-amp.sel4[max.check<0,]    
    }
    
    if(nrow(amp.sel4)==0){
      print("No Amplicon survived...Try to lower number of analyzable CpGs.")
      print(paste(date()))
      return(out)
    }
    
    print(paste("Done. ",nrow(amp.sel4), " amplicons survived.", sep=""))
  }
  
  ########################################################################################
  
  #Scan Primer:Primer Alignments      
  print("Analyze primer1 self-alignments...")
  
  #check primer1 self-alignments    
  #make all bins for each individual primer
  starts<-1:(max(nchar(amp.sel4$primer1.sequence))-primer.align.binsize+1)
  ends<-starts+primer.align.binsize
  allbins<-as.data.frame(mapply(FUN=substr,start=starts, stop=ends, 
                                MoreArgs=list(x=amp.sel4$primer1.sequence)))
  allbins<-apply(allbins,2,as.character)
  allbins<-apply(allbins,2,tolower)
  allbins<-apply(allbins,2,reverse.complement)
  ab.char<-apply(allbins,2,nchar)
  ab.char[ab.char<primer.align.binsize]<-FALSE
  ab.char[!ab.char<primer.align.binsize]<-TRUE
  allbins[ab.char==0]<-"z"
  allbins<-t(allbins)
  
  #check the overlapp of all bins to the primer1.
  #To this end the reverse complement of bins will be used.
  for (ip1 in 1:nrow(amp.sel4)){
    amp.sel4[ip1,"primer1.self.alignment"]<-max(as.numeric(sapply(allbins[,ip1],gregexpr,amp.sel4[ip1,"primer1.sequence"],
                                                                  ignore.case=TRUE)),na.rm=TRUE)
  }
  print("Done.")
  
  ######################################################################################
  
  #Scan Primer:Primer Alignments      
  print("Analyze primer2 self-alignments...")
  
  #check primer2 self-alignments    
  #make all bins for each individual primer
  starts<-1:(max(nchar(amp.sel4$primer2.sequence))-primer.align.binsize+1)
  ends<-starts+primer.align.binsize
  allbins<-as.data.frame(mapply(FUN=substr,start=starts, stop=ends, 
                                MoreArgs=list(x=amp.sel4$primer2.sequence)))
  allbins<-apply(allbins,2,as.character)
  allbins<-apply(allbins,2,tolower)
  allbins<-apply(allbins,2,reverse.complement)
  ab.char<-apply(allbins,2,nchar)
  ab.char[ab.char<primer.align.binsize]<-FALSE
  ab.char[!ab.char<primer.align.binsize]<-TRUE
  allbins[ab.char==0]<-"z"
  allbins<-t(allbins)
  
  #check the overlapp of all bins to the primer1.
  #To this end the reverse complement of primer1 will be used.
  for (ip2 in 1:nrow(amp.sel4)){
    amp.sel4[ip2,"primer2.self.alignment"]<-max(as.numeric(sapply(allbins[,ip2],gregexpr,amp.sel4[ip2,"primer2.sequence"],
                                                                  ignore.case=TRUE)),na.rm=TRUE)
  }
  print("Done.")
  
  ######################################################################################
  
  #Scan Primer:Primer Alignments      
  print("Analyze primer1:primer2 alignments...")
  
  #check primer1:primer2 alignments    
  #make all bins for each individual primer2
  starts<-1:(max(nchar(amp.sel4$primer2.sequence))-primer.align.binsize+1)
  ends<-starts+primer.align.binsize
  allbins<-as.data.frame(mapply(FUN=substr,start=starts, stop=ends, 
                                MoreArgs=list(x=amp.sel4$primer2.sequence)))
  allbins<-apply(allbins,2,as.character)
  allbins<-apply(allbins,2,tolower)
  allbins<-apply(allbins,2,reverse.complement)
  ab.char<-apply(allbins,2,nchar)
  ab.char[ab.char<primer.align.binsize]<-FALSE
  ab.char[!ab.char<primer.align.binsize]<-TRUE
  allbins[ab.char==0]<-"z"
  allbins<-t(allbins)
  
  #check the overlapp of all bins to the primer1.
  #To this end the reverse complement of primer1 will be used.
  for (ip2 in 1:nrow(amp.sel4)){
    amp.sel4[ip2,"primer1.primer2.alignment"]<-max(as.numeric(sapply(allbins[,ip2],gregexpr,amp.sel4[ip2,"primer1.sequence"],
                                                                     ignore.case=TRUE)),na.rm=TRUE)
  }
  print("Done.")
  
  ##################################################################################
  
  print("Remove amplicons that have primer with high risk for primer-dimers...")
  selection<-amp.sel4[amp.sel4$primer1.self.alignment== -1 &
                        amp.sel4$primer2.self.alignment== -1 &
                        amp.sel4$primer1.primer2.alignment== -1 &
                        (amp.sel4$amplicon.start.relative != 0 | amp.sel4$amplicon.end.relative!=0),]
  print("Done.")
  
  #####################################################################################
  
  if(nrow(selection)==0){
    print("No Primers survived filtering.")
    return(out)
    
  }
  
  #calculate some additional properties.
   selection$amplicon.sequence.genomic<-mapply(FUN=substr,seq,
                                               start = selection$amplicon.start.relative,
                                               stop = selection$amplicon.end.relative)
  selection$primer1.start.relative<-selection$amplicon.start.relative
  selection$primer1.end.relative<-selection$amplicon.start.relative + selection$primer1.length -1
  selection$primer2.start.relative<-selection$amplicon.end.relative - selection$primer1.length +1
  selection$primer2.end.relative<-selection$amplicon.end.relative 
  selection$primer1.sequence.genomic<-mapply(FUN = substr,seq,
                                             selection$primer1.start.relative,
                                             selection$primer1.end.relative)
  selection$primer2.sequence.genomic<-mapply(FUN=substr,seq,
                                             selection$primer2.start.relative,
                                             selection$primer2.end.relative) 
  selection$primer2.sequence<-reverse.complement(selection$primer2.sequence)
  selection$primer2.sequence.genomic<-reverse.complement(selection$primer2.sequence.genomic)
  selection$sequence.id<-seq.id
  selection$amplicon.id<-paste(selection$sequence.id,"_",selection$amplicon.start.relative,"_",selection$amplicon.end.relative,sep="")
  selection$dna.strand<-strand
  selection$index<-1:nrow(selection)
  selection$primer.pair.id<-paste(selection$amplicon.id,".",selection$dna.strand,".",selection$index,sep="")
  selection$primer1.id<-paste(selection$amplicon.id,".",selection$dna.strand,".",selection$index,".1",sep="")
  selection$primer2.id<-paste(selection$amplicon.id,".",selection$dna.strand,".",selection$index,".2",sep="")
  selection<-selection[order(selection$nCGs,decreasing=TRUE),]
  selection$primer1.sequence<-tolower(selection$primer1.sequence)
  selection$primer2.sequence<-tolower(selection$primer2.sequence)
  selection$amplicon.sequence<-tolower(selection$amplicon.sequence)
  
  #gc contents for primers and amplicons
  print("Calculate gc contents of primer and amplicons...")
  selection$primer1.gc.content<-nchar(gsub("[a|t]","",selection$primer1.sequence))/nchar(selection$primer1.sequence)
  selection$primer2.gc.content<-nchar(gsub("[a|t]","",selection$primer2.sequence))/nchar(selection$primer2.sequence)
  selection$amplicon.gc.content<-nchar(gsub("[a|t]","",selection$amplicon.sequence))/nchar(selection$amplicon.sequence)
  selection$amplicon.genomic.gc.content<-nchar(gsub("[a|t]","",selection$amplicon.sequence.genomic))/nchar(selection$amplicon.sequence.genomic)
  print("Done.")
  
  #c2t conversion and g2a conversion of primer1&2
  print("Calculate c2t and g2a conversion numbers of primers...")
  selection$primer1.c2t.conversion<-nchar(gsub("[a|g|t]","",selection$primer1.sequence.genomic))
  selection$primer2.g2a.conversion<-nchar(gsub("[a|c|t]","",selection$primer2.sequence.genomic))
  print("Done.")
  
  #Remove primers with low numbers of c2t/g2a conversion sites.
  print("Remove primers with few conversion sites...")
  selection<-selection[selection$primer1.c2t.conversion>=min.C2T.primer1,]
  selection<-selection[selection$primer2.g2a.conversion>=min.G2A.primer2,]
  print("Done.")
  
  #####################################################################################
  
  #finalize out object
  out<-as.data.frame(selection[,c("sequence.id","amplicon.id",
                                  "amplicon.length","amplicon.start.relative","amplicon.end.relative",
                                  "dna.strand","nGCs","nCGs",
                                  "primer1.sequence","primer2.sequence", 
                                  "primer1.length","primer2.length",
                                  "primer1.tm", "primer2.tm","tm.difference",
                                  "primer1.start.relative","primer1.end.relative",
                                  "primer2.start.relative","primer2.end.relative",
                                  "primer.pair.id",
                                  "primer1.id","primer2.id",
                                  "fragment1.id","fragment2.id","fragment12.id",
                                  "fragment12.primer.ids",
                                  "primer1.gc.content","primer2.gc.content",
                                  "amplicon.gc.content",
                                  #"amplicon.genomic.gc.content",
                                  "primer1.c2t.conversion","primer2.g2a.conversion",
                                  "primer1.self.alignment",
                                  "primer2.self.alignment",
                                  "primer1.primer2.alignment",
                                  "amplicon.sequence",
                                  "primer1.sequence.genomic","primer2.sequence.genomic",
                                  "amplicon.sequence.genomic",
                                  "index")])
  
  ##############################################################################################      
  print(paste("Completed bisulfite primer analysis for ",seq.id," on ",strand," strand: found ",nrow(out)," primer pairs.",sep=""))     
  print(paste("Finished bisulfite primer analysis for ",seq.id," on ",strand," strand.   [",date(),"].",sep="")) 
  return(out)
}#function

##############################################################################################
##############################################################################################
##############################################################################################

#NOME PRIMER DESIGN
#
nome.primer.design<-function(sequence,
                             sequence.id,
                             min.length.primer=26,
                             max.length.primer=32, 
                             min.Tm.primer=52,
                             max.Tm.primer=60, 
                             max.Tm.difference.primer=3, 
                             low.complexity.primer.removal=TRUE, 
                             max.bins.low.complexity=7,
                             remove.primers.with.n=TRUE,
                             min.C2T.primer1=3,
                             min.G2A.primer2=3,
                             min.length.amplicon=200,
                             max.length.amplicon=500, 
                             min.number.gc.amplicon=5, 
                             min.number.cg.amplicon=0,
                             primer.align.binsize=8,
                             strand="top",
                             mode="exact"){
  
  version.id<-2.0
  
  ########################################################################
  
  #Load helper scripts
  #if("pathScripts" %in% ls(name=".GlobalEnv")){
  #  source(paste(pathScripts,"_Functions\\calculate.Tm_jil_v1.0.r",sep=""))
  #  source(paste(pathScripts,"_Functions\\read.my.fasta_jil_v1.0.r",sep=""))
  #  source(paste(pathScripts,"_Functions\\reverse.complement_jil_v1.0.r",sep=""))
  #  source(paste(pathScripts,"_Functions\\gc_Conversion_jil_V1.2.r",sep=""))
  #}
  
  #######################################################################
  
  #initiate sequence
  seq<-as.character(sequence)
  print(seq)
  seq.id<-as.character(sequence.id)
  print(seq.id)
  print(paste("Start primer design for sequence: ",seq.id," on ",strand," strand...[",date(),"]",sep=""))
  print(paste("Length of Sequence: ",nchar(seq),sep=""))
  
  ########################################################################
  
  #perform bisulfite conversion of genomic sequence
  nome<-gcCon(seq,strand = strand)
  nome<-tolower(nome)
  
  #########################################################################
  
  #split the nome sequence by all "C"s.
  if(strand=="top"){
    c.split<-unlist(strsplit(x=nome,split="c"))
  } #if strand plus
  
  if (strand=="bottom"){
    split1<-strsplit(x=nome,split="gc")
    c.split<-unlist(lapply(split1,strsplit,"cg"))
  }#if    strand minus
  
  names(c.split)<-paste("fragment",1:length(c.split),sep="")
  n.fragments<-length(c.split)
  
  ########################################################################
  
  #Remove splitted sequences that are shorter than the minimum length of the desired oligos.
  c.split.length<-lapply(c.split,nchar)
  long.enough<-c.split[c.split.length>=min.length.primer]
  
  ########################################################################
  
  #Prepare NEGATIVE RESULT OUTPUT
  
  out<-data.frame(sequence.id=NA,amplicon.id=NA,amplicon.length=NA,amplicon.start.relative=NA,amplicon.end.relative=NA,
                  dna.strand=NA,nGCs=NA,nCGs=NA,primer1.sequence=NA,primer2.sequence=NA, 
                  primer1.length=NA,primer2.length=NA,primer1.tm=NA,primer2.tm=NA,tm.difference=NA,
                  primer1.start.relative=NA,primer1.end.relative=NA,primer2.start.relative=NA,
                  primer2.end.relative=NA,primer.pair.id=NA,primer1.id=NA,primer2.id=NA,fragment1.id=NA,
                  fragment2.id=NA,fragment12.id=NA,fragment12.primer.ids=NA,primer1.gc.content=NA,
                  primer2.gc.content=NA,amplicon.gc.content=NA,amplicon.genomic.gc.content=NA,
                  primer1.c2t.conversion=NA,primer2.g2a.conversion=NA,primer1.self.alignment=NA,
                  primer2.self.alignment=NA,primer1.primer2.alignment=NA,amplicon.sequence=NA,
                  primer1.sequence.genomic=NA,primer2.sequence.genomic=NA,amplicon.sequence.genomic=NA,
                  index=NA)
  
  #######################################################################
  
  if (length(long.enough)<2){
    print("No primer design possible. Try to reduce min.length.primer.")
    print(paste(date()))
    return(out)
  }
  
  #######################################################################
  
  #Caluclate the Tm for the remaining ones...
  tms<-lapply(long.enough,calculate.Tm)
  #Remove those oligos with too low Tm...
  tm.good<-as.list(long.enough[tms>min.Tm.primer])
  
  if (length(tm.good)<2){
    print("No primer design possible. Try to lower min.Tm.primer.")
    print(paste(date()))
    return(out)
  }
  
  ########################################################################
  
  #establish primers and organize them
  print("Generate & analyze individual primer sequences for each subfragment...")
  
  tmg<-data.frame(fragment.sequence=as.character(tm.good),fragment.id=names(tm.good))
  tmg$fragment.sequence<-as.character(tmg$fragment.sequence)
  tmg$fragment.length<-nchar(tmg$fragment.sequence)
  
  subfrags<-names(tm.good)[names(tm.good) %in% tmg$fragment.id]
  
  #very powerfull way! no loops here!
  #it is much quicker to use this way though it creates some wrong primers.
  #these can be removed in a later step easily.
  
  #primer length ranges to build start and end points for subfragments...
  plr<-max.length.primer-min.length.primer
  pr.starts<-1: (max(tmg$fragment.length)-min.length.primer+1)
  
  #
  if(mode=="exact"){
    st<-rep(pr.starts,each=plr+1)
    pr.dist<-min.length.primer:max.length.primer
    distrep<- rep(pr.dist,plr+1)
    en<- st+distrep-1
  }
  
  if(mode=="fast"){
    st<-rep(pr.starts,each=3)
    pr.dist<-c(min.length.primer,ceiling(mean(c(min.length.primer,max.length.primer),na.rm=TRUE)),max.length.primer)
    distrep<- rep(pr.dist,plr+1)
    en<- st+distrep-1
  }
  
  ############very powerfull!!!!
  pr.seqs<-sapply(X = tmg$fragment.sequence, function(x){
    prs<-as.data.frame(mapply(FUN=substr,start=st, stop=en, 
                              MoreArgs=list(x=x)))
  })
  ############
  pr.seqs<-sapply(pr.seqs,as.character)
  colnames(pr.seqs)<-subfrags
  pr.char<-sapply(pr.seqs,nchar)
  pr.seqs[pr.char<min.length.primer]<-"z"
  
  #convert dataframe to character for further handling. also assign primer.names.
  pr.sq<-as.character(pr.seqs)
  names(pr.sq)<-paste(rep(subfrags,each=nrow(pr.seqs)),"___",1:nrow(pr.seqs),sep="")
  
  #remove 'artificial' primers by deleting primers with length<min.length.primer
  prl<-sapply(pr.sq,nchar)
  pr.sq<-pr.sq[prl >= min.length.primer]
  print("Done.")
  
  ##################################################################
  
  #calculate primer melting temps.
  print ("Calculate primer melting temperatures...")
  pr.tm<-sapply(X = pr.sq,calculate.Tm)
  print("Done.")
  
  #################################################################
  
  #remove primers with inappropriate melting temps.
  print ("Remove primers with Tm too high or low...")
  f.good.tm<- pr.tm >= min.Tm.primer & pr.tm <= max.Tm.primer
  pr.sq<-pr.sq[f.good.tm]
  pr.tm<-pr.tm[f.good.tm]
  print("Done.")
  
  ########################################################################
  
  print ("Establish all possible primer combinations...")
  #establish primer combinations that cover enough sites of interest
  pr.pairs<-expand.grid(primer2.sequence=pr.sq,primer1.sequence=pr.sq)
  pr.pairs<-cbind(pr.pairs,expand.grid(primer2.tm=pr.tm,primer1.tm=pr.tm),
                  expand.grid(primer2.id=names(pr.tm),primer1.id=names(pr.tm)))
  print (paste("Done. Found ",nrow(pr.pairs)," primer combinations. Analyze them...",sep=""))
  
  pr.pairs$primer1.sequence<-as.character(pr.pairs$primer1.sequence)
  pr.pairs$primer2.sequence<-as.character(pr.pairs$primer2.sequence)
  pr.pairs$primer1.tm<-as.numeric(pr.pairs$primer1.tm)
  pr.pairs$primer2.tm<-as.numeric(pr.pairs$primer2.tm)
  
  #Some Basic calculations for all potential primers
  pr.pairs$primer1.length<-nchar(pr.pairs$primer1.sequence)
  pr.pairs$primer2.length<-nchar(pr.pairs$primer2.sequence)
  pr.pairs$tm.difference<-abs(pr.pairs$primer1.tm - pr.pairs$primer2.tm)
  pr.pairs$fragment1.id<-gsub("___.*","",pr.pairs$primer1.id)
  pr.pairs$fragment2.id<-gsub("___.*","",pr.pairs$primer2.id)
  print("Done.")
  
  #######
  
  # keep only if enough features are present between the two subfragments....
  print ("Remove primer combinations with too few sites of interest...")
  pr.pairs<-pr.pairs[as.numeric(gsub("fragment","",pr.pairs$fragment2.id))-
                       as.numeric(gsub("fragment","",pr.pairs$fragment1.id))>= min.number.gc.amplicon+min.number.cg.amplicon,]
  print("Done.")
  
  if(nrow(pr.pairs)<1){
    print("Could not identify a primer pair that includes sites of interest.")
    print(paste(date()))
    return(out)
  }
  
  ######################################################################
  
  #remove primer combinations that do not fullfil tm.delta criteria
  print ("Remove primer combinations with too high Tm difference...")
  pr.pairs2<-pr.pairs[pr.pairs$tm.difference <= max.Tm.difference.primer,]
  print("Done.")
  
  if(nrow(pr.pairs2)<1){
    print("Could not identify a primer pair that includes sites of interest.")
    print(paste(date()))
    return(out)
  }
  
  ######################################################################
  
  #Remove primers that have "N"s
  if(remove.primers.with.n){
    print("Remove.primers with 'N's")
    p1.n<-grepl("N",pr.pairs2[,"primer1.sequence"],ignore.case=TRUE)
    p2.n<-grepl("N",pr.pairs2[,"primer2.sequence"],ignore.case=TRUE)
    pr.pairs2<-pr.pairs2[!p1.n & !p2.n,]
  }
  
  if(nrow(pr.pairs2)<1){
    print("After removal of 'N' containing primers nothing is left.")
    print(paste(date()))
    return(out)
  }
  
  #####################################################################
  
  #Match relative position of individual primers
  print("Calculate relative primer Positions...")
  len.1<-unlist(sapply(X = pr.pairs2$primer1.sequence,regexpr,nome,useBytes=TRUE,fixed=TRUE))
  len.2<-unlist(sapply(X = pr.pairs2$primer2.sequence,regexpr,nome,useBytes=TRUE,fixed=TRUE))
  
  len.2e<-len.2+nchar(names(len.2))-1
  len<-len.2e-len.1 #v1.7: ~7 sec
  
  pr.pairs2$amplicon.length<-as.numeric(len)
  pr.pairs2$amplicon.start.relative<-len.1
  pr.pairs2$amplicon.end.relative<-len.2e
  print("Done.")
  
  #########################################################
  
  #Select amplicons by AmpliconLength
  print("Filter amplicons by min/max length...")
  amp.sel<-pr.pairs2[pr.pairs2$amplicon.length>=min.length.amplicon & 
                       pr.pairs2$amplicon.length<=max.length.amplicon,]
  
  if (nrow(amp.sel)==0){
    print("No primers found. Try to change parameters [Amplicon Length].")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel), " primer pairs survived...",sep=""))
  
  ########################################################
  
  #fetch amplicon sequences for each primer combination.
  amp.sel$amplicon.sequence<-mapply(substr,nome,
                                    start = amp.sel$amplicon.start.relative,
                                    stop = amp.sel$amplicon.end.relative)
  amp.sel2 <- amp.sel #for historic reason, actually nothing happens here.
  amp.sel2$fragment12.id<-paste(amp.sel2$fragment1.id,"_",amp.sel2$fragment2.id,sep="")
  amp.sel2$fragment12.primer.ids<-paste(amp.sel2$primer1.id,"_",amp.sel2$primer2.id,sep="")
  
  #######################################################
  
  #select best tm primer combination for individual amplicons....
  print("Select optimal primer pair combination for Tm...")
  
  #VERY POWERFULL FUNCTION:
  lvl.frgs<-levels(factor(amp.sel2$fragment12.id))
  mins<-sapply(X = lvl.frgs,FUN = function(x) {y=which.min(amp.sel2[amp.sel2$fragment12.id==x,
                                                                    "tm.difference"])
  return(amp.sel2[amp.sel2$fragment12.id==x,"fragment12.primer.ids"][y])}
  )
  amp.sel2<-amp.sel2[match(mins,amp.sel2$fragment12.primer.ids),]
  
  ###########################################################   
  
  # check amplicons for gc/cg contents
  print("Scan amplicons for information content (GpCs/CpGs)...")
  
  amp.sel2$nCGs<-sapply(X=amp.sel2$amplicon.sequence,FUN = function(x){
    cgs<-as.numeric(unlist(gregexpr("cg",text = x,ignore.case=TRUE)))
    ncgs<-length(cgs)
    if(cgs[1]== -1){
      ncgs<-0
    }
    return(ncgs)
  })
  
  amp.sel2$nGCs<-sapply(X=amp.sel2$amplicon.sequence,FUN = function(x){
    gcs<-as.numeric(unlist(gregexpr("gc",text = x,ignore.case=TRUE)))
    ngcs<-length(gcs)
    if(gcs[1]== -1){
      ngcs<-0
    }
    return(ngcs)
  })
  
  #########################################################################
  
  #filter for gc and gc contents
  print("Filter amplicons for CpGs content...")
  amp.sel3<-amp.sel2[amp.sel2$nCGs>=min.number.cg.amplicon,]
  
  if(nrow(amp.sel3)==0){
    print("No Amplicon survived...Try to lower number of analyzable CpGs.")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel3)," amplicons survived...")) 
  
  #######################################################################################
  
  #check for numbers of CpG in amplicons
  print("Filter amplicons for GpC content...")
  
  amp.sel4<-amp.sel3[amp.sel3$nGCs>=min.number.gc.amplicon,]
  
  if(nrow(amp.sel4)==0){
    print("No Amplicon survived...Try to lower number of analyzable GpCs.")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel4)," amplicons survived..."))  
  
  #######################################################################################
  
  #LOW COMPLEXITY PRIMER REMOVAL
  if(low.complexity.primer.removal) {
    print("Analyse primers for low complexities...")
    
    p1<-c("a","c","g","t")
    p1m<-sapply(p1,rep,each=max.bins.low.complexity)
    p1m<-apply(p1m,2,paste,collapse="")
    p2<-as.character(sapply(p1,FUN = paste,c(p1,p1),sep=""))
    p2m<-as.character(sapply(p2,FUN = function(x) paste(rep(x,max.bins.low.complexity),collapse="")))
    p3<-as.character(sapply(p2,FUN = paste,c(p1,p1,p1),sep=""))
    p3m<-as.character(sapply(p3,FUN = function(x) paste(rep(x,max.bins.low.complexity),collapse="")))
    
    lowlist<-c(p1m,p2m,p3m)# all evil patterns that will be checked in all primers
    #check primer 1 & 2
    if(nrow(amp.sel4)>1){
      check1<-as.data.frame(mapply(FUN=regexpr,pattern=lowlist, 
                                   MoreArgs=list(text=amp.sel4$primer1.sequence)))
      check12<-cbind(check1,mapply(FUN=regexpr,pattern=lowlist, 
                                   MoreArgs=list(text=amp.sel4$primer2.sequence)))
      max.check<-apply(check12,1,max,na.rm=T)
      
      #remove amplicons with evil patterns in primers.
      amp.sel4<-amp.sel4[max.check<0,]
    }
    
    if(nrow(amp.sel4)==1){
      check1<-max(sapply(X = lowlist,FUN = regexpr, amp.sel4[1,"primer1.sequence"]))
      check2<-max(sapply(X = lowlist,FUN = regexpr, amp.sel4[1,"primer2.sequence"]))
      max.check<-max(c(check1,check2))
      amp.sel4<-amp.sel4[max.check<0,]    
    }
    
    if(nrow(amp.sel4)==0){
      print("No Amplicon survived...Try to lower number of analyzable CpGs.")
      print(paste(date()))
      return(out)
    }
    
    print(paste("Done. ",nrow(amp.sel4), " amplicons survived.", sep=""))
  }
  
  ########################################################################################
  
  #Scan Primer:Primer Alignments      
  print("Analyze primer1 self-alignments...")
  
  #check primer1 self-alignments    
  #make all bins for each individual primer
  starts<-1:(max(nchar(amp.sel4$primer1.sequence))-primer.align.binsize+1)
  ends<-starts+primer.align.binsize
  allbins<-as.data.frame(mapply(FUN=substr,start=starts, stop=ends, 
                                MoreArgs=list(x=amp.sel4$primer1.sequence)))
  allbins<-apply(allbins,2,as.character)
  allbins<-apply(allbins,2,tolower)
  allbins<-apply(allbins,2,reverse.complement)
  ab.char<-apply(allbins,2,nchar)
  ab.char[ab.char<primer.align.binsize]<-FALSE
  ab.char[!ab.char<primer.align.binsize]<-TRUE
  allbins[ab.char==0]<-"z"
  allbins<-t(allbins)
  
  #check the overlapp of all bins to the primer1.
  #To this end the reverse complement of bins will be used.
  for (ip1 in 1:nrow(amp.sel4)){
    amp.sel4[ip1,"primer1.self.alignment"]<-max(as.numeric(sapply(allbins[,ip1],gregexpr,amp.sel4[ip1,"primer1.sequence"],
                                                                  ignore.case=TRUE)),na.rm=TRUE)
  }
  print("Done.")
  
  ######################################################################################
  
  #Scan Primer:Primer Alignments      
  print("Analyze primer2 self-alignments...")
  
  #check primer2 self-alignments    
  #make all bins for each individual primer
  starts<-1:(max(nchar(amp.sel4$primer2.sequence))-primer.align.binsize+1)
  ends<-starts+primer.align.binsize
  allbins<-as.data.frame(mapply(FUN=substr,start=starts, stop=ends, 
                                MoreArgs=list(x=amp.sel4$primer2.sequence)))
  allbins<-apply(allbins,2,as.character)
  allbins<-apply(allbins,2,tolower)
  allbins<-apply(allbins,2,reverse.complement)
  ab.char<-apply(allbins,2,nchar)
  ab.char[ab.char<primer.align.binsize]<-FALSE
  ab.char[!ab.char<primer.align.binsize]<-TRUE
  allbins[ab.char==0]<-"z"
  allbins<-t(allbins)
  
  #check the overlapp of all bins to the primer1.
  #To this end the reverse complement of primer1 will be used.
  for (ip2 in 1:nrow(amp.sel4)){
    amp.sel4[ip2,"primer2.self.alignment"]<-max(as.numeric(sapply(allbins[,ip2],gregexpr,amp.sel4[ip2,"primer2.sequence"],
                                                                  ignore.case=TRUE)),na.rm=TRUE)
  }
  print("Done.")
  
  ######################################################################################
  
  #Scan Primer:Primer Alignments      
  print("Analyze primer1:primer2 alignments...")
  
  #check primer1:primer2 alignments    
  #make all bins for each individual primer2
  starts<-1:(max(nchar(amp.sel4$primer2.sequence))-primer.align.binsize+1)
  ends<-starts+primer.align.binsize
  allbins<-as.data.frame(mapply(FUN=substr,start=starts, stop=ends, 
                                MoreArgs=list(x=amp.sel4$primer2.sequence)))
  allbins<-apply(allbins,2,as.character)
  allbins<-apply(allbins,2,tolower)
  allbins<-apply(allbins,2,reverse.complement)
  ab.char<-apply(allbins,2,nchar)
  ab.char[ab.char<primer.align.binsize]<-FALSE
  ab.char[!ab.char<primer.align.binsize]<-TRUE
  allbins[ab.char==0]<-"z"
  allbins<-t(allbins)
  
  #check the overlapp of all bins to the primer1.
  #To this end the reverse complement of primer1 will be used.
  for (ip2 in 1:nrow(amp.sel4)){
    amp.sel4[ip2,"primer1.primer2.alignment"]<-max(as.numeric(sapply(allbins[,ip2],gregexpr,amp.sel4[ip2,"primer1.sequence"],
                                                                     ignore.case=TRUE)),na.rm=TRUE)
  }
  print("Done.")
  
  ##################################################################################
  
  print("Remove amplicons that have primer with high risk for primer-dimers...")
  selection<-amp.sel4[amp.sel4$primer1.self.alignment== -1 &
                        amp.sel4$primer2.self.alignment== -1 &
                        amp.sel4$primer1.primer2.alignment== -1,]
  print("Done.")
  
  if(nrow(selection)==0){
    print("No Primers survived filtering.")
    return(out)

  }
  
  #####################################################################################
  #calculate some additional properties.
  selection$amplicon.sequence.genomic<-mapply(FUN=substr,seq,
                                               start = selection$amplicon.start.relative,
                                               stop = selection$amplicon.end.relative)
  selection$primer1.start.relative<-selection$amplicon.start.relative
  selection$primer1.end.relative<-selection$amplicon.start.relative + selection$primer1.length -1
  selection$primer2.start.relative<-selection$amplicon.end.relative - selection$primer1.length +1
  selection$primer2.end.relative<-selection$amplicon.end.relative 
  selection$primer1.sequence.genomic<-mapply(FUN = substr,seq,
                                             selection$primer1.start.relative,
                                             selection$primer1.end.relative)
  selection$primer2.sequence.genomic<-mapply(FUN=substr,seq,
                                             selection$primer2.start.relative,
                                             selection$primer2.end.relative) 
  selection$primer2.sequence<-reverse.complement(selection$primer2.sequence)
  selection$primer2.sequence.genomic<-reverse.complement(selection$primer2.sequence.genomic)
  selection$sequence.id<-seq.id
  selection$amplicon.id<-paste(selection$sequence.id,"_",selection$amplicon.start.relative,"_",selection$amplicon.end.relative,sep="")
  selection$dna.strand<-strand
  selection$index<-1:nrow(selection)
  selection$primer.pair.id<-paste(selection$amplicon.id,".",selection$dna.strand,".",selection$index,sep="")
  selection$primer1.id<-paste(selection$amplicon.id,".",selection$dna.strand,".",selection$index,".1",sep="")
  selection$primer2.id<-paste(selection$amplicon.id,".",selection$dna.strand,".",selection$index,".2",sep="")
  selection<-selection[order(selection$nCGs,decreasing=TRUE),]
  selection$primer1.sequence<-tolower(selection$primer1.sequence)
  selection$primer2.sequence<-tolower(selection$primer2.sequence)
  selection$amplicon.sequence<-tolower(selection$amplicon.sequence)
  
  #gc contents for primers and amplicons
  print("Calculate gc contents of primer and amplicons...")
  selection$primer1.gc.content<-nchar(gsub("[a|t]","",selection$primer1.sequence))/nchar(selection$primer1.sequence)
  selection$primer2.gc.content<-nchar(gsub("[a|t]","",selection$primer2.sequence))/nchar(selection$primer2.sequence)
  selection$amplicon.gc.content<-nchar(gsub("[a|t]","",selection$amplicon.sequence))/nchar(selection$amplicon.sequence)
  selection$amplicon.genomic.gc.content<-nchar(gsub("[a|t]","",selection$amplicon.sequence.genomic))/nchar(selection$amplicon.sequence.genomic)
  print("Done.")
  
  #c2t conversion and g2a conversion of primer1&2
  print("Calculate c2t and g2a conversion numbers of primers...")
  selection$primer1.c2t.conversion<-nchar(gsub("[a|g|t]","",selection$primer1.sequence.genomic))
  selection$primer2.g2a.conversion<-nchar(gsub("[a|c|t]","",selection$primer2.sequence.genomic))
  print("Done.")
  
  #Remove primers with low numbers of c2t/g2a conversion sites.
  print("Remove primers with few conversion sites...")
  selection<-selection[selection$primer1.c2t.conversion>=min.C2T.primer1,]
  selection<-selection[selection$primer2.g2a.conversion>=min.G2A.primer2,]
  print("Done.")
  
  #####################################################################################
  
  #finalize out object
  out<-as.data.frame(selection[,c("sequence.id","amplicon.id",
                                  "amplicon.length","amplicon.start.relative","amplicon.end.relative",
                                  "dna.strand","nGCs","nCGs",
                                  "primer1.sequence","primer2.sequence", 
                                  "primer1.length","primer2.length",
                                  "primer1.tm", "primer2.tm","tm.difference",
                                  "primer1.start.relative","primer1.end.relative",
                                  "primer2.start.relative","primer2.end.relative",
                                  "primer.pair.id",
                                  "primer1.id","primer2.id",
                                  "fragment1.id","fragment2.id","fragment12.id",
                                  "fragment12.primer.ids",
                                  "primer1.gc.content","primer2.gc.content",
                                  "amplicon.gc.content",
                                  #"amplicon.genomic.gc.content",
                                  "primer1.c2t.conversion","primer2.g2a.conversion",
                                  "primer1.self.alignment",
                                  "primer2.self.alignment",
                                  "primer1.primer2.alignment",
                                  "amplicon.sequence",
                                  "primer1.sequence.genomic","primer2.sequence.genomic",
                                  "amplicon.sequence.genomic",
                                  "index")])
  
  ##############################################################################################      
  print(paste(date(),".....Completed NOMe primer analysis for ",seq.id," on ",strand," strand: found ",nrow(out)," primer pairs.",sep=""))     
  
  return(out)
  
  print(paste(date(),".....Job completed.",sep="")) 
}#function

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
genomic.primer.design<-function(sequence,
                                sequence.id,
                                min.length.primer=26, 
                                max.length.primer=32, 
                                min.Tm.primer=52,
                                max.Tm.primer=60, 
                                max.Tm.difference.primer=2, 
                                low.complexity.primer.removal=TRUE, 
                                max.bins.low.complexity=7,
                                remove.primers.with.n=TRUE,
                                min.C2T.primer1=3,
                                min.G2A.primer2=3,
                                min.length.amplicon=200,
                                max.length.amplicon=500, 
                                min.number.gc.amplicon=0, 
                                min.number.cg.amplicon=7,
                                primer.align.binsize=8,
                                strand="top",
                                mode="exact",
                                chop.size=30){#chop.size splits the input by defined sizes
  
  version.id<-1.0
  
  ########################################################################
  
  #Load helper scripts
  #if("pathScripts" %in% ls(name=".GlobalEnv")){
  #  source(paste(pathScripts,"_Functions\\calculate.Tm_jil_v1.0.r",sep=""))
  #  source(paste(pathScripts,"_Functions\\read.my.fasta_jil_v1.0.r",sep=""))
  #  source(paste(pathScripts,"_Functions\\reverse.complement_jil_v1.0.r",sep=""))
  #  source(paste(pathScripts,"_Functions\\bisulfite.conversion_jil_V1.0.r",sep=""))
  #}
  
  #######################################################################
  
  #initiate sequence
  seq<-as.character(sequence)
  print(seq)
  seq.id<-as.character(sequence.id)
  print(seq.id)
  print(paste("Start Primer Design for Sequence: ",seq.id," on ",strand," strand...[",date(),"]",sep=""))
  print(paste("Length of Sequence: ",nchar(seq),sep=""))
  
  ########################################################################
  
  #perform bisulfite conversion of genomic sequence
  if(strand=="top"){nome<-tolower(seq)}
  if(strand=="bottom"){nome<-reverse.complement(tolower(seq))}
  
  #########################################################################
  
  #split the nome sequence by defined sizes.
  
  c.split <- strsplit(nome, paste("(?<=.{",chop.size,"})",sep=""), perl = TRUE)[[1]]
  names(c.split)<-paste("fragment",1:length(c.split),sep="")
  n.fragments<-length(c.split)
  
  ########################################################################
  
  #Remove splitted sequences that are shorter than the minimum length of the desired oligos.
  c.split.length<-lapply(c.split,nchar)
  long.enough<-c.split[c.split.length>=min.length.primer]
  
  ########################################################################
  
  #Prepare NEGATIVE RESULT OUTPUT
  
  out<-data.frame(sequence.id=NA,amplicon.id=NA,amplicon.length=NA,amplicon.start.relative=NA,amplicon.end.relative=NA,
                  dna.strand=NA,nGCs=NA,nCGs=NA,primer1.sequence=NA,primer2.sequence=NA, 
                  primer1.length=NA,primer2.length=NA,primer1.tm=NA,primer2.tm=NA,tm.difference=NA,
                  primer1.start.relative=NA,primer1.end.relative=NA,primer2.start.relative=NA,
                  primer2.end.relative=NA,primer.pair.id=NA,primer1.id=NA,primer2.id=NA,fragment1.id=NA,
                  fragment2.id=NA,fragment12.id=NA,fragment12.primer.ids=NA,primer1.gc.content=NA,
                  primer2.gc.content=NA,amplicon.gc.content=NA,amplicon.genomic.gc.content=NA,
                  primer1.c2t.conversion=NA,primer2.g2a.conversion=NA,primer1.self.alignment=NA,
                  primer2.self.alignment=NA,primer1.primer2.alignment=NA,amplicon.sequence=NA,
                  primer1.sequence.genomic=NA,primer2.sequence.genomic=NA,amplicon.sequence.genomic=NA,
                  index=NA)
  
  #######################################################################
  
  if (length(long.enough)<2){
    print("No Primer Design possible. Try to reduce min.length.primer.")
    print(paste(date()))
    return(out)
  }
  
  #######################################################################
  
  #CalcUlate the Tm for the remaining ones...
  tms<-lapply(long.enough,calculate.Tm)
  #Remove those oligos with too low Tm...
  tm.good<-as.list(long.enough[tms>min.Tm.primer])
  
  if (length(tm.good)<2){
    print("No Primer Design possible. Try to lower min.Tm.primer.")
    print(paste(date()))
    return(out)
  }
  
  ########################################################################
  
  #establish primers and organize them
  print("Generate & analyze individual primer sequences for each subfragment...")
  
  tmg<-data.frame(fragment.sequence=as.character(tm.good),fragment.id=names(tm.good))
  tmg$fragment.sequence<-as.character(tmg$fragment.sequence)
  tmg$fragment.length<-nchar(tmg$fragment.sequence)
  tmg$primer.start.starts<-1
  tmg$primer.start.ends<-tmg$fragment.length-min.length.primer
  tmg<-tmg[tmg$primer.start.ends>0,]
  tmg$primer.end.starts<-1+min.length.primer-1
  
  subfrags<-names(tm.good)[names(tm.good) %in% tmg$fragment.id]
  
  #very powerfull way! no loops here!
  #it is much quicker to use this way though it creates some wrong primers.
  #these can be removed in a later step easily.
  
  #primer length ranges to build start and end points for subfragments...
  plr<-max.length.primer-min.length.primer
  pr.starts<-1: (max(tmg$fragment.length)-min.length.primer+1)
  
  #
  if(mode=="exact"){
    st<-rep(pr.starts,each=plr+1)
    pr.dist<-min.length.primer:max.length.primer
    distrep<- rep(pr.dist,plr+1)
    en<- st+distrep-1
  }
  
  if(mode=="fast"){
    st<-rep(pr.starts,each=3)
    pr.dist<-c(min.length.primer,ceiling(mean(c(min.length.primer,max.length.primer),na.rm=TRUE)),max.length.primer)
    distrep<- rep(pr.dist,plr+1)
    en<- st+distrep-1
  }
  
  ############very powerfull!!!!
  pr.seqs<-sapply(X = tmg$fragment.sequence, function(x){
    prs<-as.data.frame(mapply(FUN=substr,start=st, stop=en, 
                              MoreArgs=list(x=x)))
  })
  ############
  pr.seqs<-sapply(pr.seqs,as.character)
  colnames(pr.seqs)<-subfrags
  pr.char<-sapply(pr.seqs,nchar)
  pr.seqs[pr.char<min.length.primer]<-"z"
  
  #convert dataframe to character for further handling. also assign primer.names.
  pr.sq<-as.character(pr.seqs)
  names(pr.sq)<-paste(rep(subfrags,each=nrow(pr.seqs)),"___",1:nrow(pr.seqs),sep="")
  
  #remove 'artificial' primers by deleting primers with length<min.length.primer
  prl<-sapply(pr.sq,nchar)
  pr.sq<-pr.sq[prl >= min.length.primer]
  print("Done.")
  
  ##################################################################
  
  #calculate primer melting temps.
  print ("Calculate primer melting temperatures...")
  pr.tm<-sapply(X = pr.sq,calculate.Tm)
  print("Done.")
  
  #################################################################
  
  #remove primers with inappropriate melting temps.
  print ("Remove primers with Tm too high or low...")
  f.good.tm<- pr.tm >= min.Tm.primer & pr.tm <= max.Tm.primer
  pr.sq<-pr.sq[f.good.tm]
  pr.tm<-pr.tm[f.good.tm]
  print("Done.")
  
  ########################################################################
  
  print ("Remove primer combinations that do not cover enough sites of interest...")
  #establish primer combinations that cover enough sites of interest
  pr.pairs<-expand.grid(primer2.sequence=pr.sq,primer1.sequence=pr.sq)
  pr.pairs<-cbind(pr.pairs,expand.grid(primer2.tm=pr.tm,primer1.tm=pr.tm),
                  expand.grid(primer2.id=names(pr.tm),primer1.id=names(pr.tm)))
  print (paste("Done. Found ",nrow(pr.pairs)," primer combinations. Analyze them...",sep=""))
  
  pr.pairs$primer1.sequence<-as.character(pr.pairs$primer1.sequence)
  pr.pairs$primer2.sequence<-as.character(pr.pairs$primer2.sequence)
  pr.pairs$primer1.tm<-as.numeric(pr.pairs$primer1.tm)
  pr.pairs$primer2.tm<-as.numeric(pr.pairs$primer2.tm)
  
  #Some Basic calculations for all potential primers
  pr.pairs$primer1.length<-nchar(pr.pairs$primer1.sequence)
  pr.pairs$primer2.length<-nchar(pr.pairs$primer2.sequence)
  pr.pairs$tm.difference<-abs(pr.pairs$primer1.tm - pr.pairs$primer2.tm)
  pr.pairs$fragment1.id<-gsub("___.*","",pr.pairs$primer1.id)
  pr.pairs$fragment2.id<-gsub("___.*","",pr.pairs$primer2.id)
  
  #######
  
  # keep only if enough features are present between the two subfragments....
  print ("Remove primer combinations with too few sites of interest...")
  pr.pairs<-pr.pairs[as.numeric(gsub("fragment","",pr.pairs$fragment2.id))-
                       as.numeric(gsub("fragment","",pr.pairs$fragment1.id))>= min.number.gc.amplicon+min.number.cg.amplicon,]
  
  print("Done.")
  
  if(nrow(pr.pairs)<1){
    print("Could not identify a primer pair that includes sites of interest.")
    print(paste(date()))
    return(out)
  }
  
  ######################################################################
  
  #remove primer combinations that do not fullfil tm.delta criteria
  print ("Remove primer combinations with too high Tm difference...")
  pr.pairs2<-pr.pairs[pr.pairs$tm.difference <= max.Tm.difference.primer,]
  print("Done.")
  
  if(nrow(pr.pairs2)<1){
    print("Could not identify a primer pair that includes sites of interest.")
    print(paste(date()))
    return(out)
  }
  
  ######################################################################
  
  #Remove primers that have "N"s
  if(remove.primers.with.n){
    print("Remove.primers with 'N's")
    p1.n<-grepl("N",pr.pairs2[,"primer1.sequence"],ignore.case=TRUE)
    p2.n<-grepl("N",pr.pairs2[,"primer2.sequence"],ignore.case=TRUE)
    pr.pairs2<-pr.pairs2[!p1.n & !p2.n,]
  }
  
  if(nrow(pr.pairs2)<1){
    print("After removal of 'N' containing primers nothing is left.")
    print(paste(date()))
    return(out)
  }
  
  #####################################################################
  
  #Match relative position of individual primers
  print("Calculate relative primer Positions...")
  len.1<-unlist(sapply(X = pr.pairs2$primer1.sequence,regexpr,nome,useBytes=TRUE,fixed=TRUE))
  len.2<-unlist(sapply(X = pr.pairs2$primer2.sequence,regexpr,nome,useBytes=TRUE,fixed=TRUE))
  
  len.2e<-len.2+nchar(names(len.2))-1
  len<-len.2e-len.1 #v1.7: ~7 sec
  
  pr.pairs2$amplicon.length<-as.numeric(len)
  pr.pairs2$amplicon.start.relative<-len.1
  pr.pairs2$amplicon.end.relative<-len.2e
  print("Done.")
  
  #########################################################
  
  #Select amplicons by AmpliconLength
  print("Filter amplicons by min/max length...")
  amp.sel<-pr.pairs2[pr.pairs2$amplicon.length>=min.length.amplicon & 
                       pr.pairs2$amplicon.length<=max.length.amplicon,]
  
  if (nrow(amp.sel)==0){
    print("No primers found. Try to change parameters [Amplicon Length].")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel), " primer pairs survived...",sep=""))
  
  ########################################################
  
  #fetch amplicon sequences for each primer combination.
  amp.sel$amplicon.sequence<-mapply(substr,nome,
                                    start = amp.sel$amplicon.start.relative,
                                    stop = amp.sel$amplicon.end.relative)
  amp.sel2 <- amp.sel #for historic reason, actually nothing happens here.
  amp.sel2$fragment12.id<-paste(amp.sel2$fragment1.id,"_",amp.sel2$fragment2.id,sep="")
  amp.sel2$fragment12.primer.ids<-paste(amp.sel2$primer1.id,"_",amp.sel2$primer2.id,sep="")
  
  #######################################################
  
  #select best tm primer combination for individual amplicons....
  print("Select optimal primer pair combination for Tm...")
  
  #VERY POWERFUL FUNCTION:
  lvl.frgs<-levels(factor(amp.sel2$fragment12.id))
  mins<-sapply(X = lvl.frgs,FUN = function(x) {y=which.min(amp.sel2[amp.sel2$fragment12.id==x,
                                                                    "tm.difference"])
  return(amp.sel2[amp.sel2$fragment12.id==x,"fragment12.primer.ids"][y])}
  )
  amp.sel2<-amp.sel2[match(mins,amp.sel2$fragment12.primer.ids),]
  
  ###########################################################   
  
  # check amplicons for gc/cg contents
  print("Scan amplicons for information content (GpCs/CpGs)...")
  
  amp.sel2$nCGs<-sapply(X=amp.sel2$amplicon.sequence,FUN = function(x){
    cgs<-as.numeric(unlist(gregexpr("cg",text = x,ignore.case=TRUE)))
    ncgs<-length(cgs)
    if(cgs[1]== -1){
      ncgs<-0
    }
    return(ncgs)
  })
  
  amp.sel2$nGCs<-sapply(X=amp.sel2$amplicon.sequence,FUN = function(x){
    gcs<-as.numeric(unlist(gregexpr("gc",text = x,ignore.case=TRUE)))
    ngcs<-length(gcs)
    if(gcs[1]== -1){
      ngcs<-0
    }
    return(ngcs)
  })
  
  #########################################################################
  
  #filter for gc and gc contents
  print("Filter amplicons for CpGs content...")
  amp.sel3<-amp.sel2[amp.sel2$nCGs>=min.number.cg.amplicon,]
  
  if(nrow(amp.sel3)==0){
    print("No Amplicon survived...Try to lower number of analyzable CpGs.")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel3)," amplicons survived...")) 
  
  #######################################################################################
  
  #check for numbers of CpG in amplicons
  print("Filter amplicons for GpC content...")
  
  amp.sel4<-amp.sel3[amp.sel3$nGCs>=min.number.gc.amplicon,]
  
  if(nrow(amp.sel4)==0){
    print("No Amplicon survived...Try to lower number of analyzable GpCs.")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel4)," amplicons survived..."))  
  
  #######################################################################################
  
  #LOW COMPLEXITY PRIMER REMOVAL
  if(low.complexity.primer.removal) {
    print("Analyse primers for low complexities...")
    
    p1<-c("a","c","g","t")
    p1m<-sapply(p1,rep,each=max.bins.low.complexity)
    p1m<-apply(p1m,2,paste,collapse="")
    p2<-as.character(sapply(p1,FUN = paste,c(p1,p1),sep=""))
    p2m<-as.character(sapply(p2,FUN = function(x) paste(rep(x,max.bins.low.complexity),collapse="")))
    p3<-as.character(sapply(p2,FUN = paste,c(p1,p1,p1),sep=""))
    p3m<-as.character(sapply(p3,FUN = function(x) paste(rep(x,max.bins.low.complexity),collapse="")))
    
    lowlist<-c(p1m,p2m,p3m)# all evil patterns that will be checked in all primers
    
    #check primer 1 & 2
    if(nrow(amp.sel4)>1){
      check1<-as.data.frame(mapply(FUN=regexpr,pattern=lowlist, 
                                   MoreArgs=list(text=amp.sel4$primer1.sequence)))
      check12<-cbind(check1,mapply(FUN=regexpr,pattern=lowlist, 
                                   MoreArgs=list(text=amp.sel4$primer2.sequence)))
      max.check<-apply(check12,1,max,na.rm=T)
      
      #remove amplicons with evil patterns in primers.
      amp.sel4<-amp.sel4[max.check<0,]
    }
    
    if(nrow(amp.sel4)==1){
      check1<-max(sapply(X = lowlist,FUN = regexpr, amp.sel4[1,"primer1.sequence"]))
      check2<-max(sapply(X = lowlist,FUN = regexpr, amp.sel4[1,"primer2.sequence"]))
      max.check<-max(c(check1,check2))
      amp.sel4<-amp.sel4[max.check<0,]    
    }
    
    if(nrow(amp.sel4)==0){
      print("No Amplicon survived...Try to lower number of analyzable CpGs.")
      print(paste(date()))
      return(out)
    }
    
    print(paste("Done. ",nrow(amp.sel4), " amplicons survived.", sep=""))
  }
  
  ########################################################################################
  
  #Scan Primer:Primer Alignments      
  print("Analyze primer1 self-alignments...")
  
  #check primer1 self-alignments    
  #make all bins for each individual primer
  starts<-1:(max(nchar(amp.sel4$primer1.sequence))-primer.align.binsize+1)
  ends<-starts+primer.align.binsize
  allbins<-as.data.frame(mapply(FUN=substr,start=starts, stop=ends, 
                                MoreArgs=list(x=amp.sel4$primer1.sequence)))
  allbins<-apply(allbins,2,as.character)
  allbins<-apply(allbins,2,tolower)
  allbins<-apply(allbins,2,reverse.complement)
  ab.char<-apply(allbins,2,nchar)
  ab.char[ab.char<primer.align.binsize]<-FALSE
  ab.char[!ab.char<primer.align.binsize]<-TRUE
  allbins[ab.char==0]<-"z"
  allbins<-t(allbins)
  
  #check the overlapp of all bins to the primer1.
  #To this end the reverse complement of bins will be used.
  for (ip1 in 1:nrow(amp.sel4)){
    amp.sel4[ip1,"primer1.self.alignment"]<-max(as.numeric(sapply(allbins[,ip1],gregexpr,amp.sel4[ip1,"primer1.sequence"],
                                                                  ignore.case=TRUE)),na.rm=TRUE)
  }
  print("Done.")
  
  ######################################################################################
  
  #Scan Primer:Primer Alignments      
  print("Analyze primer2 self-alignments...")
  
  #check primer2 self-alignments    
  #make all bins for each individual primer
  starts<-1:(max(nchar(amp.sel4$primer2.sequence))-primer.align.binsize+1)
  ends<-starts+primer.align.binsize
  allbins<-as.data.frame(mapply(FUN=substr,start=starts, stop=ends, 
                                MoreArgs=list(x=amp.sel4$primer2.sequence)))
  allbins<-apply(allbins,2,as.character)
  allbins<-apply(allbins,2,tolower)
  allbins<-apply(allbins,2,reverse.complement)
  ab.char<-apply(allbins,2,nchar)
  ab.char[ab.char<primer.align.binsize]<-FALSE
  ab.char[!ab.char<primer.align.binsize]<-TRUE
  allbins[ab.char==0]<-"z"
  allbins<-t(allbins)
  
  #check the overlapp of all bins to the primer1.
  #To this end the reverse complement of primer1 will be used.
  for (ip2 in 1:nrow(amp.sel4)){
    amp.sel4[ip2,"primer2.self.alignment"]<-max(as.numeric(sapply(allbins[,ip2],gregexpr,amp.sel4[ip2,"primer2.sequence"],
                                                                  ignore.case=TRUE)),na.rm=TRUE)
  }
  print("Done.")
  
  ######################################################################################
  
  #Scan Primer:Primer Alignments      
  print("Analyze primer1:primer2 alignments...")
  
  #check primer1:primer2 alignments    
  #make all bins for each individual primer2
  starts<-1:(max(nchar(amp.sel4$primer2.sequence))-primer.align.binsize+1)
  ends<-starts+primer.align.binsize
  allbins<-as.data.frame(mapply(FUN=substr,start=starts, stop=ends, 
                                MoreArgs=list(x=amp.sel4$primer2.sequence)))
  allbins<-apply(allbins,2,as.character)
  allbins<-apply(allbins,2,tolower)
  allbins<-apply(allbins,2,reverse.complement)
  ab.char<-apply(allbins,2,nchar)
  ab.char[ab.char<primer.align.binsize]<-FALSE
  ab.char[!ab.char<primer.align.binsize]<-TRUE
  allbins[ab.char==0]<-"z"
  allbins<-t(allbins)
  
  #check the overlapp of all bins to the primer1.
  #To this end the reverse complement of primer1 will be used.
  for (ip2 in 1:nrow(amp.sel4)){
    amp.sel4[ip2,"primer1.primer2.alignment"]<-max(as.numeric(sapply(allbins[,ip2],gregexpr,amp.sel4[ip2,"primer1.sequence"],
                                                                     ignore.case=TRUE)),na.rm=TRUE)
  }
  print("Done.")
  
  ##################################################################################
  
  print("Remove amplicons that have primer with high risk for primer-dimers...")
  selection<-amp.sel4[amp.sel4$primer1.self.alignment== -1 &
                        amp.sel4$primer2.self.alignment== -1 &
                        amp.sel4$primer1.primer2.alignment== -1,]
  print("Done.")
  
  if(nrow(selection)==0){
    print("No Primers survived filtering.")
    return(out)
    
  }
  
  #####################################################################################
  
  #calculate some additional properties.
  selection$amplicon.sequence.genomic<-mapply(FUN=substr,seq,
                                               start = selection$amplicon.start.relative,
                                               stop = selection$amplicon.end.relative)
  selection$primer1.start.relative<-selection$amplicon.start.relative
  selection$primer1.end.relative<-selection$amplicon.start.relative + selection$primer1.length -1
  selection$primer2.start.relative<-selection$amplicon.end.relative - selection$primer1.length +1
  selection$primer2.end.relative<-selection$amplicon.end.relative 
  selection$primer1.sequence.genomic<-mapply(FUN = substr,seq,
                                             selection$primer1.start.relative,
                                             selection$primer1.end.relative)
  selection$primer2.sequence.genomic<-mapply(FUN=substr,seq,
                                             selection$primer2.start.relative,
                                             selection$primer2.end.relative) 
  selection$primer2.sequence<-reverse.complement(selection$primer2.sequence)
  selection$primer2.sequence.genomic<-reverse.complement(selection$primer2.sequence.genomic)
  selection$sequence.id<-seq.id
  selection$amplicon.id<-paste(selection$sequence.id,"_",selection$amplicon.start.relative,"_",selection$amplicon.end.relative,sep="")
  selection$dna.strand<-strand
  selection$index<-1:nrow(selection)
  selection$primer.pair.id<-paste(selection$amplicon.id,".",selection$dna.strand,".",selection$index,sep="")
  selection$primer1.id<-paste(selection$amplicon.id,".",selection$dna.strand,".",selection$index,".1",sep="")
  selection$primer2.id<-paste(selection$amplicon.id,".",selection$dna.strand,".",selection$index,".2",sep="")
  selection<-selection[order(selection$nCGs,decreasing=TRUE),]
  selection$primer1.sequence<-tolower(selection$primer1.sequence)
  selection$primer2.sequence<-tolower(selection$primer2.sequence)
  selection$amplicon.sequence<-tolower(selection$amplicon.sequence)
  
  #gc contents for primers and amplicons
  print("Calculate gc contents of primer and amplicons...")
  selection$primer1.gc.content<-nchar(gsub("[a|t]","",selection$primer1.sequence))/nchar(selection$primer1.sequence)
  selection$primer2.gc.content<-nchar(gsub("[a|t]","",selection$primer2.sequence))/nchar(selection$primer2.sequence)
  selection$amplicon.gc.content<-nchar(gsub("[a|t]","",selection$amplicon.sequence))/nchar(selection$amplicon.sequence)
  selection$amplicon.genomic.gc.content<-nchar(gsub("[a|t]","",selection$amplicon.sequence.genomic))/nchar(selection$amplicon.sequence.genomic)
  print("Done.")
  
  #c2t conversion and g2a conversion of primer1&2
  print("Calculate c2t and g2a conversion numbers of primers...")
  selection$primer1.c2t.conversion<-nchar(gsub("[a|g|t]","",selection$primer1.sequence.genomic))
  selection$primer2.g2a.conversion<-nchar(gsub("[a|c|t]","",selection$primer2.sequence.genomic))
  print("Done.")
  
  #Remove primers with low numbers of c2t/g2a conversion sites.
  print("Remove primers with few conversion sites...")
  selection<-selection[selection$primer1.c2t.conversion>=min.C2T.primer1,]
  selection<-selection[selection$primer2.g2a.conversion>=min.G2A.primer2,]
  print("Done.")
  
  #####################################################################################
  
  #finalize out object

  out<-as.data.frame(selection[,c("sequence.id","amplicon.id",#vorhanden
                                  "amplicon.length","amplicon.start.relative","amplicon.end.relative",#vorhanden
                                  "dna.strand","nGCs","nCGs",#vorhanden
                                  "primer1.sequence","primer2.sequence",#vorhanden 
                                  "primer1.length","primer2.length",#vorhanden
                                  "primer1.tm", "primer2.tm","tm.difference",#vorhanden
                                  "primer1.start.relative","primer1.end.relative",#vorhanden
                                  "primer2.start.relative","primer2.end.relative",#vorhanden
                                  "primer.pair.id",#vorhanden
                                  "primer1.id","primer2.id",#vorhanden
                                  "fragment1.id","fragment2.id","fragment12.id",#vorhanden
                                  "fragment12.primer.ids",#vorhanden
                                  "primer1.gc.content","primer2.gc.content",#vorhanden
                                  "amplicon.gc.content",#vorhanden
                                  #"amplicon.genomic.gc.content",
                                  "primer1.c2t.conversion","primer2.g2a.conversion",#vorhanden
                                  "primer1.self.alignment",#vorhanden
                                  "primer2.self.alignment",#vorhanden
                                  "primer1.primer2.alignment",#vorhanden
                                  "amplicon.sequence",#vorhanden
                                  "primer1.sequence.genomic","primer2.sequence.genomic",#vorhanden
                                  "amplicon.sequence.genomic",
                                  "index")])#vorhanden
  
  ##############################################################################################      
  print(paste("Completed genomic primer analysis for ",seq.id," on ",strand," strand: found ",nrow(out)," primer pairs.",sep=""))     
  print(paste("Finished genomic primer analysis for ",seq.id," on ",strand," strand.   [",date(),"].",sep="")) 
  return(out)
}#function

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################


CLEVER.primer.design<-function(sequence,
                               sequence.id,
                               min.length.primer=26, 
                               max.length.primer=32, 
                               min.Tm.primer=52,
                               max.Tm.primer=60, 
                               max.Tm.difference.primer=2, 
                               low.complexity.primer.removal=TRUE, 
                               max.bins.low.complexity=7,
                               remove.primers.with.n=TRUE,
                               min.C2T.primer1=3,
                               min.G2A.primer2=3,
                               min.length.amplicon=200,
                               max.length.amplicon=500, 
                               min.number.gc.amplicon=0, 
                               min.number.cg.amplicon=7,
                               primer.align.binsize=8,
                               strand="top",
                               mode="exact"){
  
  version.id<-2.0
  
  ########################################################################
  
  #Load helper scripts
  if("pathScripts" %in% ls(name=".GlobalEnv")){
    source(paste(pathScripts,"_Functions\\calculate.Tm_jil_v1.0.r",sep=""))
    source(paste(pathScripts,"_Functions\\read.my.fasta_jil_v1.0.r",sep=""))
    source(paste(pathScripts,"_Functions\\reverse.complement_jil_v1.0.r",sep=""))
    #  source(paste(pathScripts,"_Functions\\bisulfite.conversion_jil_V1.0.r",sep=""))
  }
  
  #######################################################################
  
  #initiate sequence
  seq<-as.character(sequence)
  print(seq)
  seq.id<-as.character(sequence.id)
  print(seq.id)
  print(paste("Start Primer Design for Sequence: ",seq.id," on ",strand," strand...[",date(),"]",sep=""))
  print(paste("Length of Sequence: ",nchar(seq),sep=""))
  
  ########################################################################
  
  #do NOT! perform bisulfite conversion of genomic sequence
  nome<-seq
  nome<-tolower(nome)
  
  #########################################################################
  
  #split the nome sequence by all "cg"s.
  c.split<-unlist(strsplit(x=nome,split="cg"))
  names(c.split)<-paste("fragment",1:length(c.split),sep="")
  n.fragments<-length(c.split)
  
  ########################################################################
  
  #Remove splitted sequences that are shorter than the minimum length of the desired oligos.
  c.split.length<-lapply(c.split,nchar)
  long.enough<-c.split[c.split.length>=min.length.primer]
  
  ########################################################################
  
  #Prepare NEGATIVE RESULT OUTPUT
  
  out<-data.frame(sequence.id=NA,amplicon.id=NA,amplicon.length=NA,amplicon.start.relative=NA,amplicon.end.relative=NA,
                  dna.strand=NA,nGCs=NA,nCGs=NA,primer1.sequence=NA,primer2.sequence=NA, 
                  primer1.length=NA,primer2.length=NA,primer1.tm=NA,primer2.tm=NA,tm.difference=NA,
                  primer1.start.relative=NA,primer1.end.relative=NA,primer2.start.relative=NA,
                  primer2.end.relative=NA,primer.pair.id=NA,primer1.id=NA,primer2.id=NA,fragment1.id=NA,
                  fragment2.id=NA,fragment12.id=NA,fragment12.primer.ids=NA,primer1.gc.content=NA,
                  primer2.gc.content=NA,amplicon.gc.content=NA,amplicon.genomic.gc.content=NA,
                  primer1.c2t.conversion=NA,primer2.g2a.conversion=NA,primer1.self.alignment=NA,
                  primer2.self.alignment=NA,primer1.primer2.alignment=NA,amplicon.sequence=NA,
                  primer1.sequence.genomic=NA,primer2.sequence.genomic=NA,amplicon.sequence.genomic=NA,
                  index=NA)
  
  #######################################################################
  
  if (length(long.enough)<2){
    print("No Primer Design possible. Try to reduce min.length.primer.")
    print(paste(date()))
    return(out)
  }
  
  #######################################################################
  
  #CalcUlate the Tm for the remaining ones...
  tms<-lapply(long.enough,calculate.Tm)
  #Remove those oligos with too low Tm...
  tm.good<-as.list(long.enough[tms>min.Tm.primer])
  
  if (length(tm.good)<2){
    print("No Primer Design possible. Try to lower min.Tm.primer.")
    print(paste(date()))
    return(out)
  }
  
  ########################################################################
  
  #establish primers and organize them
  print("Generate & analyze individual primer sequences for each subfragment...")
  
  tmg<-data.frame(fragment.sequence=as.character(tm.good),fragment.id=names(tm.good))
  tmg$fragment.sequence<-as.character(tmg$fragment.sequence)
  tmg$fragment.length<-nchar(tmg$fragment.sequence)
  tmg$primer.start.starts<-1
  tmg$primer.start.ends<-tmg$fragment.length-min.length.primer
  tmg<-tmg[tmg$primer.start.ends>0,]
  tmg$primer.end.starts<-1+min.length.primer-1
  
  subfrags<-names(tm.good)[names(tm.good) %in% tmg$fragment.id]
  
  #very powerfull way! no loops here!
  #it is much quicker to use this way though it creates some wrong primers.
  #these can be removed in a later step easily.
  
  #primer length ranges to build start and end points for subfragments...
  plr<-max.length.primer-min.length.primer
  pr.starts<-1: (max(tmg$fragment.length)-min.length.primer+1)
  
  #
  if(mode=="exact"){
    st<-rep(pr.starts,each=plr+1)
    pr.dist<-min.length.primer:max.length.primer
    distrep<- rep(pr.dist,plr+1)
    en<- st+distrep-1
  }
  
  if(mode=="fast"){
    st<-rep(pr.starts,each=3)
    pr.dist<-c(min.length.primer,ceiling(mean(c(min.length.primer,max.length.primer),na.rm=TRUE)),max.length.primer)
    distrep<- rep(pr.dist,plr+1)
    en<- st+distrep-1
  }
  
  ############very powerfull!!!!
  pr.seqs<-sapply(X = tmg$fragment.sequence, function(x){
    prs<-as.data.frame(mapply(FUN=substr,start=st, stop=en, 
                              MoreArgs=list(x=x)))
  })
  ############
  pr.seqs<-sapply(pr.seqs,as.character)
  colnames(pr.seqs)<-subfrags
  pr.char<-sapply(pr.seqs,nchar)
  pr.seqs[pr.char<min.length.primer]<-"z"
  
  #convert dataframe to character for further handling. also assign primer.names.
  pr.sq<-as.character(pr.seqs)
  names(pr.sq)<-paste(rep(subfrags,each=nrow(pr.seqs)),"___",1:nrow(pr.seqs),sep="")
  
  #remove 'artificial' primers by deleting primers with length<min.length.primer
  prl<-sapply(pr.sq,nchar)
  pr.sq<-pr.sq[prl >= min.length.primer]
  print("Done.")
  
  ##################################################################
  
  #calculate primer melting temps.
  print ("Calculate primer melting temperatures...")
  pr.tm<-sapply(X = pr.sq,calculate.Tm)
  print("Done.")
  
  #################################################################
  
  #remove primers with inappropriate melting temps.
  print ("Remove primers with Tm too high or low...")
  f.good.tm<- pr.tm >= min.Tm.primer & pr.tm <= max.Tm.primer
  pr.sq<-pr.sq[f.good.tm]
  pr.tm<-pr.tm[f.good.tm]
  print("Done.")
  
  ########################################################################
  
  print ("Remove primer combinations that do not cover enough sites of interest...")
  #establish primer combinations that cover enough sites of interest
  pr.pairs<-expand.grid(primer2.sequence=pr.sq,primer1.sequence=pr.sq)
  pr.pairs<-cbind(pr.pairs,expand.grid(primer2.tm=pr.tm,primer1.tm=pr.tm),
                  expand.grid(primer2.id=names(pr.tm),primer1.id=names(pr.tm)))
  print (paste("Done. Found ",nrow(pr.pairs)," primer combinations. Analyze them...",sep=""))
  
  pr.pairs$primer1.sequence<-as.character(pr.pairs$primer1.sequence)
  pr.pairs$primer2.sequence<-as.character(pr.pairs$primer2.sequence)
  pr.pairs$primer1.tm<-as.numeric(pr.pairs$primer1.tm)
  pr.pairs$primer2.tm<-as.numeric(pr.pairs$primer2.tm)
  
  #Some Basic calculations for all potential primers
  pr.pairs$primer1.length<-nchar(pr.pairs$primer1.sequence)
  pr.pairs$primer2.length<-nchar(pr.pairs$primer2.sequence)
  pr.pairs$tm.difference<-abs(pr.pairs$primer1.tm - pr.pairs$primer2.tm)
  pr.pairs$fragment1.id<-gsub("___.*","",pr.pairs$primer1.id)
  pr.pairs$fragment2.id<-gsub("___.*","",pr.pairs$primer2.id)
  
  #######
  
  # keep only if enough features are present between the two subfragments....
  print ("Remove primer combinations with too few sites of interest...")
  pr.pairs<-pr.pairs[as.numeric(gsub("fragment","",pr.pairs$fragment2.id))-
                       as.numeric(gsub("fragment","",pr.pairs$fragment1.id))>= min.number.gc.amplicon+min.number.cg.amplicon,]
  
  print("Done.")
  
  if(nrow(pr.pairs)<1){
    print("Could not identify a primer pair that includes sites of interest.")
    print(paste(date()))
    return(out)
  }
  
  ######################################################################
  
  #remove primer combinations that do not fullfil tm.delta criteria
  print ("Remove primer combinations with too high Tm difference...")
  pr.pairs2<-pr.pairs[pr.pairs$tm.difference <= max.Tm.difference.primer,]
  print("Done.")
  
  if(nrow(pr.pairs2)<1){
    print("Could not identify a primer pair that includes sites of interest.")
    print(paste(date()))
    return(out)
  }
  
  ######################################################################
  
  #Remove primers that have "N"s
  if(remove.primers.with.n){
    print("Remove.primers with 'N's")
    p1.n<-grepl("N",pr.pairs2[,"primer1.sequence"],ignore.case=TRUE)
    p2.n<-grepl("N",pr.pairs2[,"primer2.sequence"],ignore.case=TRUE)
    pr.pairs2<-pr.pairs2[!p1.n & !p2.n,]
  }
  
  if(nrow(pr.pairs2)<1){
    print("After removal of 'N' containing primers nothing is left.")
    print(paste(date()))
    return(out)
  }
  
  #####################################################################
  
  #Match relative position of individual primers
  print("Calculate relative primer Positions...")
  len.1<-unlist(sapply(X = pr.pairs2$primer1.sequence,regexpr,nome,useBytes=TRUE,fixed=TRUE))
  len.2<-unlist(sapply(X = pr.pairs2$primer2.sequence,regexpr,nome,useBytes=TRUE,fixed=TRUE))
  
  len.2e<-len.2+nchar(names(len.2))-1
  len<-len.2e-len.1 #v1.7: ~7 sec
  
  pr.pairs2$amplicon.length<-as.numeric(len)
  pr.pairs2$amplicon.start.relative<-len.1
  pr.pairs2$amplicon.end.relative<-len.2e
  print("Done.")
  
  #########################################################
  
  #Select amplicons by AmpliconLength
  print("Filter amplicons by min/max length...")
  amp.sel<-pr.pairs2[pr.pairs2$amplicon.length>=min.length.amplicon & 
                       pr.pairs2$amplicon.length<=max.length.amplicon,]
  
  if (nrow(amp.sel)==0){
    print("No primers found. Try to change parameters [Amplicon Length].")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel), " primer pairs survived...",sep=""))
  
  ########################################################
  
  #fetch amplicon sequences for each primer combination.
  amp.sel$amplicon.sequence<-mapply(substr,nome,
                                    start = amp.sel$amplicon.start.relative,
                                    stop = amp.sel$amplicon.end.relative)
  amp.sel2 <- amp.sel #for historic reason, actually nothing happens here.
  amp.sel2$fragment12.id<-paste(amp.sel2$fragment1.id,"_",amp.sel2$fragment2.id,sep="")
  amp.sel2$fragment12.primer.ids<-paste(amp.sel2$primer1.id,"_",amp.sel2$primer2.id,sep="")
  
  #######################################################
  
  #select best tm primer combination for individual amplicons....
  print("Select optimal primer pair combination for Tm...")
  
  #VERY POWERFULL FUNCTION:
  lvl.frgs<-levels(factor(amp.sel2$fragment12.id))
  mins<-sapply(X = lvl.frgs,FUN = function(x) {y=which.min(amp.sel2[amp.sel2$fragment12.id==x,
                                                                    "tm.difference"])
  return(amp.sel2[amp.sel2$fragment12.id==x,"fragment12.primer.ids"][y])}
  )
  amp.sel2<-amp.sel2[match(mins,amp.sel2$fragment12.primer.ids),]
  
  ###########################################################   
  
  # check amplicons for gc/cg contents
  print("Scan amplicons for information content (GpCs/CpGs)...")
  
  amp.sel2$nCGs<-sapply(X=amp.sel2$amplicon.sequence,FUN = function(x){
    cgs<-as.numeric(unlist(gregexpr("cg",text = x,ignore.case=TRUE)))
    ncgs<-length(cgs)
    if(cgs[1]== -1){
      ncgs<-0
    }
    return(ncgs)
  })
  
  amp.sel2$nGCs<-sapply(X=amp.sel2$amplicon.sequence,FUN = function(x){
    gcs<-as.numeric(unlist(gregexpr("gc",text = x,ignore.case=TRUE)))
    ngcs<-length(gcs)
    if(gcs[1]== -1){
      ngcs<-0
    }
    return(ngcs)
  })
  
  #########################################################################
  
  #filter for gc and gc contents
  print("Filter amplicons for CpGs content...")
  amp.sel3<-amp.sel2[amp.sel2$nCGs>=min.number.cg.amplicon,]
  
  if(nrow(amp.sel3)==0){
    print("No Amplicon survived...Try to lower number of analyzable CpGs.")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel3)," amplicons survived...")) 
  
  #######################################################################################
  
  #check for numbers of CpG in amplicons
  print("Filter amplicons for GpC content...")
  
  amp.sel4<-amp.sel3[amp.sel3$nGCs>=min.number.gc.amplicon,]
  
  if(nrow(amp.sel4)==0){
    print("No Amplicon survived...Try to lower number of analyzable GpCs.")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel4)," amplicons survived..."))  
  
  #######################################################################################
  
  #LOW COMPLEXITY PRIMER REMOVAL
  if(low.complexity.primer.removal) {
    print("Analyse primers for low complexities...")
    
    p1<-c("a","c","g","t")
    p1m<-sapply(p1,rep,each=max.bins.low.complexity)
    p1m<-apply(p1m,2,paste,collapse="")
    p2<-as.character(sapply(p1,FUN = paste,c(p1,p1),sep=""))
    p2m<-as.character(sapply(p2,FUN = function(x) paste(rep(x,max.bins.low.complexity),collapse="")))
    p3<-as.character(sapply(p2,FUN = paste,c(p1,p1,p1),sep=""))
    p3m<-as.character(sapply(p3,FUN = function(x) paste(rep(x,max.bins.low.complexity),collapse="")))
    
    lowlist<-c(p1m,p2m,p3m)# all evil patterns that will be checked in all primers
    
    #check primer 1 & 2
    if(nrow(amp.sel4)>1){
      check1<-as.data.frame(mapply(FUN=regexpr,pattern=lowlist, 
                                   MoreArgs=list(text=amp.sel4$primer1.sequence)))
      check12<-cbind(check1,mapply(FUN=regexpr,pattern=lowlist, 
                                   MoreArgs=list(text=amp.sel4$primer2.sequence)))
      max.check<-apply(check12,1,max,na.rm=T)
      
      #remove amplicons with evil patterns in primers.
      amp.sel4<-amp.sel4[max.check<0,]
    }
    
    if(nrow(amp.sel4)==1){
      check1<-max(sapply(X = lowlist,FUN = regexpr, amp.sel4[1,"primer1.sequence"]))
      check2<-max(sapply(X = lowlist,FUN = regexpr, amp.sel4[1,"primer2.sequence"]))
      max.check<-max(c(check1,check2))
      amp.sel4<-amp.sel4[max.check<0,]    
    }
    
    if(nrow(amp.sel4)==0){
      print("No Amplicon survived...Try to lower number of analyzable CpGs.")
      print(paste(date()))
      return(out)
    }
    
    print(paste("Done. ",nrow(amp.sel4), " amplicons survived.", sep=""))
  }
  
  ########################################################################################
  
  #Scan Primer:Primer Alignments      
  print("Analyze primer1 self-alignments...")
  
  #check primer1 self-alignments    
  #make all bins for each individual primer
  starts<-1:(max(nchar(amp.sel4$primer1.sequence))-primer.align.binsize+1)
  ends<-starts+primer.align.binsize
  allbins<-as.data.frame(mapply(FUN=substr,start=starts, stop=ends, 
                                MoreArgs=list(x=amp.sel4$primer1.sequence)))
  allbins<-apply(allbins,2,as.character)
  allbins<-apply(allbins,2,tolower)
  allbins<-apply(allbins,2,reverse.complement)
  ab.char<-apply(allbins,2,nchar)
  ab.char[ab.char<primer.align.binsize]<-FALSE
  ab.char[!ab.char<primer.align.binsize]<-TRUE
  allbins[ab.char==0]<-"z"
  allbins<-t(allbins)
  
  #check the overlapp of all bins to the primer1.
  #To this end the reverse complement of bins will be used.
  for (ip1 in 1:nrow(amp.sel4)){
    amp.sel4[ip1,"primer1.self.alignment"]<-max(as.numeric(sapply(allbins[,ip1],gregexpr,amp.sel4[ip1,"primer1.sequence"],
                                                                  ignore.case=TRUE)),na.rm=TRUE)
  }
  print("Done.")
  
  ######################################################################################
  
  #Scan Primer:Primer Alignments      
  print("Analyze primer2 self-alignments...")
  
  #check primer2 self-alignments    
  #make all bins for each individual primer
  starts<-1:(max(nchar(amp.sel4$primer2.sequence))-primer.align.binsize+1)
  ends<-starts+primer.align.binsize
  allbins<-as.data.frame(mapply(FUN=substr,start=starts, stop=ends, 
                                MoreArgs=list(x=amp.sel4$primer2.sequence)))
  allbins<-apply(allbins,2,as.character)
  allbins<-apply(allbins,2,tolower)
  allbins<-apply(allbins,2,reverse.complement)
  ab.char<-apply(allbins,2,nchar)
  ab.char[ab.char<primer.align.binsize]<-FALSE
  ab.char[!ab.char<primer.align.binsize]<-TRUE
  allbins[ab.char==0]<-"z"
  allbins<-t(allbins)
  
  #check the overlapp of all bins to the primer1.
  #To this end the reverse complement of primer1 will be used.
  for (ip2 in 1:nrow(amp.sel4)){
    amp.sel4[ip2,"primer2.self.alignment"]<-max(as.numeric(sapply(allbins[,ip2],gregexpr,amp.sel4[ip2,"primer2.sequence"],
                                                                  ignore.case=TRUE)),na.rm=TRUE)
  }
  print("Done.")
  
  ######################################################################################
  
  #Scan Primer:Primer Alignments      
  print("Analyze primer1:primer2 alignments...")
  
  #check primer1:primer2 alignments    
  #make all bins for each individual primer2
  starts<-1:(max(nchar(amp.sel4$primer2.sequence))-primer.align.binsize+1)
  ends<-starts+primer.align.binsize
  allbins<-as.data.frame(mapply(FUN=substr,start=starts, stop=ends, 
                                MoreArgs=list(x=amp.sel4$primer2.sequence)))
  allbins<-apply(allbins,2,as.character)
  allbins<-apply(allbins,2,tolower)
  allbins<-apply(allbins,2,reverse.complement)
  ab.char<-apply(allbins,2,nchar)
  ab.char[ab.char<primer.align.binsize]<-FALSE
  ab.char[!ab.char<primer.align.binsize]<-TRUE
  allbins[ab.char==0]<-"z"
  allbins<-t(allbins)
  
  #check the overlapp of all bins to the primer1.
  #To this end the reverse complement of primer1 will be used.
  for (ip2 in 1:nrow(amp.sel4)){
    amp.sel4[ip2,"primer1.primer2.alignment"]<-max(as.numeric(sapply(allbins[,ip2],gregexpr,amp.sel4[ip2,"primer1.sequence"],
                                                                     ignore.case=TRUE)),na.rm=TRUE)
  }
  print("Done.")
  
  ##################################################################################
  
  print("Remove amplicons that have primer with high risk for primer-dimers...")
  selection<-amp.sel4[amp.sel4$primer1.self.alignment== -1 &
                        amp.sel4$primer2.self.alignment== -1 &
                        amp.sel4$primer1.primer2.alignment== -1,]
  print("Done.")
  
  #####################################################################################
  
  #calculate some additional properties.
  #calculations are only possible, if the selection list contains more than just the names of the variables
   
  if(length(selection) > 22){
      selection$amplicon.sequence.genomic<-mapply(FUN=substr,seq,
                                                  start = selection$amplicon.start.relative,
                                                  stop = selection$amplicon.end.relative)
  }
  
      selection$primer1.start.relative<-selection$amplicon.start.relative
      selection$primer1.end.relative<-selection$amplicon.start.relative + selection$primer1.length -1
      selection$primer2.start.relative<-selection$amplicon.end.relative - selection$primer1.length +1
      selection$primer2.end.relative<-selection$amplicon.end.relative 
    if (length(selection$primer1.start.relative) != 0 && length(selection$primer1.start.relative) != 0){
      
      selection$primer1.sequence.genomic<-mapply(FUN = substr,seq,
                                                 selection$primer1.start.relative,
                                                 selection$primer1.end.relative)
      selection$primer2.sequence.genomic<-mapply(FUN=substr,seq,
                                                 selection$primer2.start.relative,
                                                 selection$primer2.end.relative) 
    }
      selection$primer2.sequence<-reverse.complement(selection$primer2.sequence)
      selection$primer2.sequence.genomic<-reverse.complement(selection$primer2.sequence.genomic)
      #TODO: seq.id non-existent --> where to get the id?
      selection["sequence.id"] <- seq.id
      #selection$sequence.id<-seq.id 
      selection$amplicon.id<-paste(selection$sequence.id,"_",selection$amplicon.start.relative,"_",selection$amplicon.end.relative,sep="")
      selection$dna.strand<-strand
      selection$index<-1:nrow(selection)
      selection$primer.pair.id<-paste(selection$amplicon.id,".",selection$dna.strand,".",selection$index,sep="")
      selection$primer1.id<-paste(selection$amplicon.id,".",selection$dna.strand,".",selection$index,".1",sep="")
      selection$primer2.id<-paste(selection$amplicon.id,".",selection$dna.strand,".",selection$index,".2",sep="")
      selection<-selection[order(selection$nCGs,decreasing=TRUE),]
      selection$primer1.sequence<-tolower(selection$primer1.sequence)
      selection$primer2.sequence<-tolower(selection$primer2.sequence)
      selection$amplicon.sequence<-tolower(selection$amplicon.sequence)
      
      #gc contents for primers and amplicons
      print("Calculate gc contents of primer and amplicons...")
      selection$primer1.gc.content<-nchar(gsub("[a|t]","",selection$primer1.sequence))/nchar(selection$primer1.sequence)
      selection$primer2.gc.content<-nchar(gsub("[a|t]","",selection$primer2.sequence))/nchar(selection$primer2.sequence)
      selection$amplicon.gc.content<-nchar(gsub("[a|t]","",selection$amplicon.sequence))/nchar(selection$amplicon.sequence)
      #selection$amplicon.genomic.gc.content<-nchar(gsub("[a|t]","",selection$amplicon.sequence.genomic))/nchar(selection$amplicon.sequence.genomic)
      print("Done.")
      
      #c2t conversion and g2a conversion of primer1&2
      print("Calculate c2t and g2a conversion numbers of primers...")
      selection$primer1.c2t.conversion<-nchar(gsub("[a|g|t]","",selection$primer1.sequence.genomic))
      selection$primer2.g2a.conversion<-nchar(gsub("[a|c|t]","",selection$primer2.sequence.genomic))
      print("Done.")
      
      #Remove primers with low numbers of c2t/g2a conversion sites.
      print("Remove primers with few conversion sites...")
      selection<-selection[selection$primer1.c2t.conversion>=min.C2T.primer1,]
      selection<-selection[selection$primer2.g2a.conversion>=min.G2A.primer2,]
      print("Done.")

    
  
  #####################################################################################
  
  #finalize out object
  print(selection)
  out<-as.data.frame(selection[,c("sequence.id","amplicon.id",#vorhanden
                                  "amplicon.length","amplicon.start.relative","amplicon.end.relative",#vorhanden
                                  "dna.strand",#vorhanden
                                  "nGCs","nCGs",#vorhanden
                                  "primer1.sequence","primer2.sequence", #vorhanden
                                  "primer1.length","primer2.length",#vorhanden
                                  "primer1.tm", "primer2.tm","tm.difference",#vorhanden
                                  "primer1.start.relative","primer1.end.relative",#vorhanden
                                  "primer2.start.relative","primer2.end.relative",#vorhanden
                                  "primer.pair.id",#vorhanden
                                  "primer1.id","primer2.id",#vorhanden
                                  "fragment1.id","fragment2.id","fragment12.id",#vorhanden
                                  "fragment12.primer.ids",#vorhanden
                                  "primer1.gc.content","primer2.gc.content",#vorhanden
                                  "amplicon.gc.content",#vorhanden
                                  #"amplicon.genomic.gc.content",
                                  "primer1.c2t.conversion","primer2.g2a.conversion",#vorhanden
                                  "primer1.self.alignment",#vorhanden
                                  "primer2.self.alignment",#vorhanden
                                  "primer1.primer2.alignment",#vorhanden
                                  "amplicon.sequence",#vorhanden
                                  "primer1.sequence.genomic","primer2.sequence.genomic",#vorhanden
                                  #"amplicon.sequence.genomic",
                                  "index")])#vorhanden
  
  ##############################################################################################      
  print(paste("Completed CLEVER primer analysis for ",seq.id," on ",strand," strand: found ",nrow(out)," primer pairs.",sep=""))     
  print(paste("Finished CLEVER primer analysis for ",seq.id," on ",strand," strand.   [",date(),"].",sep="")) 
  return(out)
}#function

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################






#get SNP info for a genomic intervall 
# will call ensembl rest API (http://mar2017.rest.ensembl.org/)
#

fetch.snp.info.rest = function(assembly = NULL, #'hg19' or 'hg38'
                               chr = NULL, #'chr1' or '1'
                               start = NULL, #'12345678'
                               end = NULL,   #'12345679'
                               ...){
  require(httr)
  require(jsonlite)
  require(xml2)
  
  iassembly = switch(EXPR = assembly, 
                     "hg19" = "GRCh37",
                     "hg38" = "GRCh38",
                     "mm10" = "GRCm38")
  
  ispecies = switch(EXPR = assembly, 
                    "hg19" = "human",
                    "hg38" = "human",
                    "mm10" = "mus_musculus")
  
  server <- switch(EXPR = assembly, 
                   "hg19" = "http://grch37.rest.ensembl.org",
                   "hg38" = "http://rest.ensembl.org",
                   "mm10" = "http://rest.ensembl.org")
  
  ext <- paste0("/overlap/region/",ispecies,"/",gsub("chr","",chr),":",start,"-",end,":1?content-type=text/plain;feature=variation")
  
  if(length(ext) != 0){
    r <- GET(paste(server, ext[1], sep=""))
    s <- content(r)
  }
  
  if(length(ext) >= 2){
    for (i in 2:length(ext)){
      r <- GET(paste(server, ext[i], sep=""))
      stop_for_status(r)
      s <- rbind(s, content(r))
    }
  }
  
  s2 = lapply(s,function(x) {tdf = data.frame(c1 = names(unlist(x)), c2 = unlist(x))})
  s3 = sapply(1:length(s2), function(x) { s2[x][[1]] = as.data.frame(s2[x][[1]][c("seq_region_name","start","end","strand","assembly_name","id","feature_type","consequence_type","source"),"c2"])})
  s4 = as.data.frame(matrix(unlist(s3),nrow = length(s3), ncol = 9, byrow = TRUE))
  colnames(s4) = c("chr","start","end","strand","assembly","rs_id","feature_type","consequence_type","source")
  return(s4)
  
}

##########

##example
#snp = fetch.snp.info.rest(assembly = "hg38",
#                          chr = "chr7",
#                          start = "140424943",
#                          end = "140426943")

##########

#get repeat info for a genomic intervall 
# will call ensembl rest API (http://mar2017.rest.ensembl.org/)
#

fetch.repeat.info.rest = function(assembly = NULL, #'hg19' or 'hg38'
                                  chr = NULL, #'chr1' or '1'
                                  start = NULL, #'12345678'
                                  end = NULL,   #'12345679'
                                  ...){
  require(httr)
  require(jsonlite)
  require(xml2)
  
  iassembly = switch(EXPR = assembly, 
                     "hg19" = "GRCh37",
                     "hg38" = "GRCh38",
                     "mm10" = "GRCm38")
  
  ispecies = switch(EXPR = assembly, 
                    "hg19" = "human",
                    "hg38" = "human",
                    "mm10" = "mus_musculus")
  
  server <- switch(EXPR = assembly, 
                   "hg19" = "http://grch37.rest.ensembl.org",
                   "hg38" = "http://rest.ensembl.org",
                   "mm10" = "http://rest.ensembl.org")
  
  ext <- paste0("/overlap/region/",ispecies,"/",gsub("chr","",chr),":",start,"-",end,":1?content-type=text/plain;feature=repeat")
  
  if(length(ext) != 0){
    r <- GET(paster(server, ext[1], sep=""))
    s <- content(r)
  }
  
  if(length(ext) >=2){
    for(i in 2:length(ext)){
      r <- GET(paste(server, ext[i], sep=""))
      stop_for_status(r)
      s <- rbind(s, content(r))
    }
  }
  
  s2 = lapply(s,function(x) {tdf = data.frame(c1 = names(unlist(x)), c2 = unlist(x))})
  s3 = sapply(1:length(s2), function(x) { s2[x][[1]] = as.data.frame(s2[x][[1]][c("seq_region_name","start","end","strand","assembly_name","feature_type","description"),"c2"])})
  s4 = as.data.frame(matrix(unlist(s3),nrow = length(s3), ncol = 7, byrow = TRUE))
  colnames(s4) = c("chr","start","end","strand","assembly","feature_type","description")
  
  return(s4)
  
}

##########

##example
#reps = fetch.repeat.info.rest(assembly = "hg38",
#                             chr = "chr7",
#                             start = "140424943",
#                             end = "140426943")#

#reps
##########

#more info here: http://rest.ensembl.org/documentation/info/sequence_region

#returns a character with the requested sequence

fetch.dna.sequence= function(assembly = NULL, #'hg19' or 'hg38'
                             chr = NULL, #'chr1' or '1'
                             start = NULL, #'12345678'
                             end = NULL,   #'12345679'
                             ...){
  require(httr)
  require(jsonlite)
  require(xml2)
  
  iassembly = switch(EXPR = assembly, 
                     "hg19" = "GRCh37",
                     "hg38" = "GRCh38",
                     "mm10" = "GRCm38")
  ispecies = switch(EXPR = assembly, 
                    "hg19" = "human",
                    "hg38" = "human",
                    "mm10" = "mus_musculus")
  
  server <- "http://rest.ensembl.org"
  ext <- paste0("/sequence/region/",ispecies,"/",gsub("chr","",chr),":",start,"..",end,":1?;coord_system_version=",iassembly)
  
  if(length(ext) != 0){
    r <- GET(paste(server, ext[1], sep=""), content_type("text/plain"))
    s <- content(r)
    
  }
  
  for (i in 2:length(ext)){
    r <- GET(paste(server, ext[i], sep = ""), content_type("text/plain"))
    s <- rbind(s, content(r))
  }
  
  #r <- GET(paste(server, ext, sep = ""), content_type("text/plain"))
  
  #stop_for_status(r)
  return(s)
  
}


#example
# fetch.dna.sequence(assembly = "hg19",
#                    chr="chr1",
#                    start="100000000",
#                    end="100001000")





#get gene info for a genomic intervall 
# will call ensembl rest API (http://mar2017.rest.ensembl.org/)
#

fetch.gene.info.rest = function(assembly = NULL, #'hg19' or 'hg38'
                                chr = NULL, #'chr1' or '1'
                                start = NULL, #'12345678'
                                end = NULL,   #'12345679'
                                ...){
  require(httr)
  require(jsonlite)
  require(xml2)
  
  iassembly = switch(EXPR = assembly, 
                     "hg19" = "GRCh37",
                     "hg38" = "GRCh38",
                     "mm10" = "GRCm38")
  
  ispecies = switch(EXPR = assembly, 
                    "hg19" = "human",
                    "hg38" = "human",
                    "mm10" = "mus_musculus")
  
  server <- switch(EXPR = assembly, 
                   "hg19" = "http://grch37.rest.ensembl.org",
                   "hg38" = "http://rest.ensembl.org",
                   "mm10" = "http://rest.ensembl.org")
  
  ext <- paste0("/overlap/region/",ispecies,"/",gsub("chr","",chr),":",start,"-",end,":1?content-type=text/plain;feature=gene")
  if(length(ext) != 0){
    r <- GET(paste(server, ext[1], sep=""))
    s <- content(r)
  }
  
  if(length(ext) >= 2){
    for (i in 2:length(ext)){
      r <- GET(paste(server, ext[i], sep=""))
      stop_for_status(r)
      s <- rbind(s, content(r))
    }
  }
  
  if(length(s)>0){
    s2 = lapply(s,function(x) {tdf = data.frame(c1 = names(unlist(x)), c2 = unlist(x))})
    s3 = sapply(1:length(s2), function(x) {s2[x][[1]] = as.data.frame(s2[x][[1]][c("seq_region_name","start","end",
                                                                                   "strand","assembly_name","id",
                                                                                   "external_name","feature_type",
                                                                                   "biotype","version","description",
                                                                                   "source"),"c2"])})
    s4 = as.data.frame(matrix(unlist(s3),nrow = length(s3), ncol = 12, byrow = TRUE))
    colnames(s4) = c("chr","start","end","strand","assembly","id","gene_name","feature_type","biotype",
                     "version","description","source")
  }
  
  if(length(s) == 0){
    s4 = data.frame(chr=NULL, start= NULL, end=NULL, assembly=NULL, id=NULL, gene_name= NULL,
                    feature_type=NULL, biotype=NULL, description=NULL, source=NULL)
  }
  
  return(s4)
  
}

##########

##example
#gene = fetch.gene.info.rest(assembly = "hg38",
#                          chr = "chr7",
#                          start = "140424943",
#                          end = "140426943")

##########

