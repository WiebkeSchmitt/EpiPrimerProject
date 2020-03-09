source("primer.design.r")
source("HelperFunctions.R")

primer.design.pipeline<-function(table.in,#filename.in = NULL, # direct path to the file with i) regions (.bed format) or ii) sequences (filepath) (character)
                                 input.type="regions",#"regions" (for a file with regions in .bed format) or "sequences" (for a ".txt" file with sequences) (fixed vocabulary) (character)
                                 path.out = NULL, # directory where results will be stored (directory) (character)
                                 analysis.id="PrimerAnalysis",# Do not use'SonderZeichen'
                                 primer.type="bisulfite",#"NOME" or "bisulfite" or "hp_bisulfite" or "hp_NOME" or "genomic" or "hp_genomic" or "CLEVER" or "hp_CLEVER" (fixed vocabulary) (character)
                                 mode="fast", #"exact" or "fast" fixed vocabulary) (character)
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
                                 min.MAF.snp = 0.01,   #max allowed minor allele frequency for SNPs
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
  
  #######################################Create amplicon barplots##########################################
  create.amplicon.barplots=FALSE
  #######################################Parallelization##########################################

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
  #######################################External packages##########################################
  
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
  
  #######################################R2HTML##########################################
  require(R2HTML)
  
  #######################################Setting up working directory##########################################
  
  dir.create(path.out)
  path.wd<-paste(path.out,analysis.id,"/",sep="/")
  dir.create(path.wd)
  setwd(path.wd)
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
  
  #########################################Setting up logfile########################################
  fn.log<-paste(path.wd,analysis.id,"#logfile.txt",sep="")
  log<-function(text,add=TRUE){
    write.log(paste(text,sep=""),file=fn.log,add=add)
  }
  log("Logfile initiated.",add=FALSE)
  ##########################################Checking primer parameters#######################################

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
  
  
  ######################################Adjustments of chop size for genomic primers###########################################
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
  ##################################Create overview .html file###############################################

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
  #########################################Data import########################################
  
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
    
    ###########################################Removing entries w/o sequence retrieved############################################
    
    bad.seq.f<-bed$sequence.length<1  
    bad.out<-bed[bad.seq.f,]
    
    if(nrow(bad.out)>0){
      write.table(bad.out,file=paste(path.sequences,analysis.id,"EntriesWithNoSequences.txt",sep=""),sep="\t",dec=".",
                  col.names=T,row.names=F)  
    }
    
    bed<-bed[!bad.seq.f,]
    
    ###########################################perform input sequence characterization########################################
    
    log("Start characterization of input sequences...")
    
    seq.ana.all<-list()
    
    for (iiii in 1:nrow(bed)){
      
      seqanai<-sequence.characterization(paste(bed[iiii,"sequence"]))
      write.table(seqanai,file=paste(path.sequences,"sequence.characterization_",bed[iiii,"sequenceID"],".txt",sep=""),
                  sep="\t",dec=".",col.names=T,row.names = FALSE,quote=F)
      
      seq.ana.all[[paste(bed[iiii,"sequenceID"])]]<-seqanai
    }
    
    log("Done.")
    
    #############################################retrieve snp information for input sequences##########################################
    
    if (check4snps){
      
      supported.assemblies.snp<-c("hg18","hg19","hg38","mm9","mm10")
      
      log("Start SNP analysis...")
      
      if(input.type=="sequences"){
        #should check biomart....
        log(paste("WARNING: To check for SNPs input type must be regions."))
        fileConn<-file(paste(path.warnings,"Warning.SNPcheck.data.error.txt",sep=""))
        writeLines(c("WARNING: To check for SNPs input type must be regions."), fileConn)
        close(fileConn)
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
    
    
    ##############################################retrieve repeat information for input sequences#########################################
    
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
    
    #############################################retrieve gene information for input sequences  ##########################################
    
    if (annotate.genes){
      
      supported.assemblies.snp<-c("hg18","hg19","hg38","mm9","mm10")
      
      log("Start Gene annotation analysis...")
      
      if(input.type=="sequences"){
        
        #should check biomart....
        log(paste("WARNING: To check for gene annotations input type must be regions."))
        fileConn<-file(paste(path.warnings,"Warning.GeneAnnotationcheck.data.error.txt",sep=""))
        writeLines(c("WARNING: To check for Gene Annotation input type must be regions."), fileConn)
        close(fileConn)
        hp.initial.input.type<-input.type
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
    
    ###########################################retrieve cpg island information for input sequences  ############################################
    
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
            
            chrom<-paste(bedia[,"chr"])
            allstarts<-as.numeric(as.character(bedia[,"start"]))
            allends<-as.numeric(as.character(bedia[,"end"]))
            
            cpgitable<-switch(assem, 
                              "hg18" = "cpgIslandExt",
                              "hg19" = "cpgIslandExt",
                              "hg38" = "cpgIslandExt",
                              "mm9"  = "cpgIslandExt",
                              "mm10" = "cpgIslandExt")
            
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
    
    #################################changes according to primer type######################################################
    
    if(primer.type == "bisulfite" | primer.type == "hp_bisulfite" | primer.type == "CrispRCas9PCR"){
      bed$sequence.converted<-sapply(bed$sequence,bisulfite.conversion,strand=strand)
    }
    
    if(primer.type == "NOME" | primer.type == "hp_NOME"){
      bed$sequence.converted<-sapply(bed$sequence,gcCon,strand=strand)   
    }
    
    if(primer.type == "genomic" | primer.type == "hp_genomic" | primer.type == "CLEVER" | primer.type == "hp_CLEVER"){
      bed$sequence.converted<-bed$sequence
    }
    
    hp.initial.input.type<-input.type
    
    bed$Index<-1:nrow(bed)
    write.table(bed,file=paste(path.sequences,analysis.id,"_PrimerDesign_InputSequences.txt",sep=""),sep="\t",dec=".",
                col.names=TRUE,row.names=FALSE,quote=F)
    log("Done.")
  }
  
  if (input.type=="sequences"){
    
    log("Import Sequences from .txt File...")
    bed<-table.in #read.table(file=filename.in,sep="\t",dec=".",header=TRUE)
    #html.report(filenames =  filename.in,filename.out = paste(path.html,"bedfile_",analysis.id,".html",sep=""),txt.header = TRUE,txt.sep = "\t")
    log("Done.")
    
    
    ###########################################perform input sequence characterization########################################
    
    log("Start characterization of input sequences...")
    
    seq.ana.all<-list()
    
    for (iiii in 1:nrow(bed)){
      
      seqanai<-sequence.characterization(paste(bed[iiii,"sequence"]))
      write.table(seqanai,file=paste(path.sequences,"sequence.characterization_",bed[iiii,"sequenceID"],".txt",sep=""),
                  sep="\t",dec=".",col.names=T,row.names = FALSE,quote=F)
      
      seq.ana.all[[paste(bed[iiii,"sequenceID"])]]<-seqanai
    }
    
    log("Done.")
    
    #############################################retrieve snp information for input sequences##########################################
    
    
    
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
  ######################################hairpin specific section: calculate restriction digests##########################################
  
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
                                  "Bsu361"),
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
    
    hp.res.ids<-paste(hp.info[,"RE.id"])
    hp.res.seqs<-paste(hp.info[,"RE.seq"])
    hp.linker.seqs<-paste(hp.info[,"linker.sequence"])
    
    ##############################################Module: Hairpinizer######################################## 
    
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
    
    ################################################complete results################################
    
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
    
    #########################################pre filter for primer design.#######################################
    
    hpf<-hp[!is.na(hp$hp_Region_Length),]
    hpf<-hpf[hpf$hp_Region_Length>=hp.length.min & hpf$hp_Region_Length<=hp.length.max,]
    hpf$hp_Sequence<-as.character(hpf$hp_Sequence)
    hpf$linker_Sequence<-as.character(hpf$linker_Sequence)
    
    hpf$LinkerStartInMolecule<-nchar(hpf$hp_Sequence)+1
    hpf$LinkerEndInMolecule<-hpf$LinkerStartInMolecule+nchar(hpf$linker_Sequence)-1
    
    ##########################################export temporal results######################################
    
    log("Export temporal results.")
    write.table(hpf,file=paste(path.wd,analysis.id,"_Hairpinizer_Results_PreFiltered.txt",sep=""),
                sep="\t",dec=".",col.names=TRUE,row.names=FALSE,quote=F)
    
    ############################################Filter for target regions...####################################
    
    hpf2<-hpf
    
    ############################################if necessary re-assign target start/end column ids.####################################
    
    target.start.column.id<-"LinkerStartInMolecule"
    target.end.column.id<-"LinkerEndInMolecule"
    use.target.columns<-TRUE
    
    ############################################export of results####################################
    
    log("Export temporal results.")
    write.table(hpf2,file=paste(path.wd,analysis.id,"_Hairpinizer_Results_Final.txt",sep=""),
                sep="\t",dec=".",col.names=TRUE,row.names=FALSE,quote=F)
    
    log("Done.")
    
    ########################################Bedfiile swapping########################################
    
    log("Bedfile swapping...")
    write.table(bed,file=paste(path.sequences,analysis.id,"_PrimerDesign_Initial.InputSequences.txt",sep=""),sep="\t",dec=".",
                col.names=TRUE,row.names=FALSE,quote=F)
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
    
    ############################################perform input sequence characterization on hairpin regions####################################    
    
    log("Start characterization of input sequences...")
    
    for (iiii in 1:nrow(bed)){
      
      seqanai<-sequence.characterization(paste(bed[iiii,"sequence"]))
      write.table(seqanai,file=paste(path.sequences,"sequence.characterization_",bed[iiii,"sequenceID"],".txt",sep=""),
                  sep="\t",dec=".",col.names=T,row.names = FALSE,quote=F)
      
      seq.ana.all[[paste(bed[iiii,"sequenceID"])]]<-seqanai
    }
    
    log("Done.")
    
  }#if hp_bisulfite or hp_NOMe
  
  #############################################strand selection###################################
  
  if (strand=="both"){
    strand<-c("top","bottom")
  }
  
  # disable strand selection for hairpin analyses
  if(primer.type=="hp_bisulfite" | primer.type=="hp_NOME" | primer.type=="hp_genomic"){
    strand<-"top"
    log("Enforced strand-setting to top-strand since it is a hairpin analysis.")
  }
  
  ###############################################primer.design section####################################
  
  index<-1

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
      
      ######################################primer.design results post-processing######################################
      
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
  
  #############################################Processing results######################################
  
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
    
    ############################################Export results#####################################
    
    log(paste("Export final temporal results..."))
    write.table(results2,file=paste(path.temp,analysis.id,"_TempResults.final.txt",sep=""),
                col.names=T,row.names=T,sep="\t",dec=".",quote=F)
    log("Done.")
    
    #############################################Filter na long amplicons.#######################################
    
    results2<-results2[!is.na(results2$amplicon.length),]
    if(nrow(results2)>0){
      
      ###############################################Checking for SNPs####################################

      if (check4snps){
        
        log("Assign SNP information to designed amplicons...")
        
        if(exists(x = "all_my_snps")){
          
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
            
          }# iamplicons
          
          log("Done.")
          
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
      
      ##############################################Checking for Repeats#####################################
      
      if (check4repeats){
        
        log("Assign repeat information to designed amplicons...")
        
        if(exists(x = "all_my_repeats")){
          
          #results2$repeat.db<-"Repeatmasker"
          # results2$amplicon.n.repeats<-NA
          # results2$primer1.n.repeats<-NA
          # results2$primer2.n.repeats<-NA
          # results2$amplicon.repeat.ids<-NA
          # results2$primer1.repeat.ids<-NA
          # results2$primer2.repeat.ids<-NA
          
          if(nrow(results2) >= 1){
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
                results2[iamplicons,"amplicon.repeat.ids"] <- paste(amp.reps[,"description"],collapse=",")
                results2[iamplicons,"primer1.n.repeats"] <- nrow(p1.reps)
                results2[iamplicons,"primer1.repeat.ids"] <- paste(p1.reps[,"description"],collapse=",")
                results2[iamplicons,"primer2.n.repeats"] <- nrow(p2.reps)
                results2[iamplicons,"primer2.repeat.ids"] <- paste(p2.reps[,"description"],collapse=",")
                
              }# iamplicons
          }
          log("Done.")
          
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
      
      ##########################################Annotate genes#########################################
      
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
            print("amp.genes")
            print(amp.genes)
            results2[iamplicons,"amplicon.n.genes"] <- nrow(amp.genes)
            results2[iamplicons,"amplicon.gene.ids"] <- paste(amp.genes[,"gene_name"],collapse=",")
            
          }# iamplicons
          
          log("Done.")
          
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
      
      ###############################################Annotate CPG islands####################################
      
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
      
      ################################################Remove results that do not cover target region.###################################
      
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
      
      #####################################################remove double results#####################################
      
      #fix v9.8: remove double results
      temp.p1p2.ids<-paste0(results2$amplicon.id,"_",results2$primer1.length,"_",results2$primer2.length)
      results2<-results2[!duplicated(temp.p1p2.ids),]
      write.table(results2,file=paste(path.temp,analysis.id,"_TempResults.final_after.dup.removal.txt",sep=""),
                  col.names=T,row.names=T,sep="\t",dec=".",quote=F)
      
      #####################################################create columns with added NGS adaptors to primers.#####################################
      
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
      
      ################################################export results###################################
      
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
      
      ##############################################for all amplicons generate .fasta file#####################################
      
      afn<-paste(path.sequences,results2$amplicon.id,".FASTA",sep="")
      mapply(FUN = write.fasta,results2$amplicon.sequence.genomic,results2$amplicon.id,afn)
      results3<-results2
      
      ##############################################create toplist tables#####################################
      
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
        
        ##############################################Export toplist tables#####################################
        
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
      
      ################################################generate ucsc browser custom tracks for each sequence###################################
      
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
      ################################################generate custom tracks for BiSearch tool.###################################
      
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
      
      ###################################################generate .bed files for amplicon and/or primers################################
      
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
      
      #############################################Generate summary######################################
      
      log("Generate summary...")
      end.run<-paste(Sys.time())
      smry<-data.frame(PARAMETERS=c("Analysis Start","Analysis End","Number of Input Sequences",
                                    "Number of PrimerPairs [Total]","Number of Sequences with Primers","Number of Sequences without Primers"),
                       SETTINGS=c(start.run,end.run, paste(nrow(bed)),paste(n.total.primer),
                                  paste(nrow(toplist)),paste(nrow(bed)-nrow(toplist))))
      write.table(smry,file=paste(path.wd,"summary_",analysis.id,".txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",dec=".")
      #html.report(filenames =  paste(path.html,"summary_",analysis.id,".txt",sep=""),filename.out = paste(path.html,"summary_",analysis.id,".html",sep=""),txt.header = TRUE,txt.sep = "\t")
      
      # primer designs by input sequence
      if(!is.na(results2$sequence.id)){
        tbl<-data.frame(table(results2$sequence.id))
        colnames(tbl)<-c("sequence.id","amplicons[n]")
        write.table(tbl,file=paste(path.wd,"primerdesigns.by.sequence_",analysis.id,".txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",dec=".")
        #html.report(filenames =  paste(path.html,"primerdesigns.by.sequence_",analysis.id,".txt",sep=""),filename.out = paste(path.html,"summary_",analysis.id,".html",sep=""),txt.header = TRUE,txt.sep = "\t")
      }
      log("Done.")
      
      ###################################################Write Whitelist################################
      
      log("Write whitelist...")
      white<-bed[paste(bed$sequenceID) %in% paste(toplist[,1]),]
      
      if(!file.exists(paste(path.wd,"primer_",analysis.id,"_whitelist.txt",sep=""))){
        write.table(white,file=paste(path.wd,"primer_",analysis.id,"_whitelist.txt",sep=""),
                    col.names=TRUE,row.names=FALSE,sep="\t",dec=".") 
      }
      
      if(file.exists(paste(path.wd,"primer_",analysis.id,"_whitelist.txt",sep=""))){
        white.old<-read.table(paste(path.wd,"primer_",analysis.id,"_whitelist.txt",sep=""),header=TRUE,sep="\t",dec=".")
        white<-rbind(white.old,white)
        write.table(white,file=paste(path.wd,"primer_",analysis.id,"_whitelist.txt",sep=""),
                    col.names=TRUE,row.names=FALSE,sep="\t",dec=".") 
      }
      
      log("Done.")
      
      log("Write blacklist...")
      black<-bed[!paste(bed$sequenceID) %in% paste(white$sequenceID),]
      write.table(black,file=paste(path.wd,"primer_",analysis.id,"_blacklist.txt",sep=""),
                  col.names=TRUE,row.names=FALSE,sep="\t",dec=".") 
      log("Done.")
      
      ###################################################Create graphics################################

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
              
              ####################################################per sequence id primer design overview plot###########################################################
              
              if(primer.type == "bisulfite" | primer.type == "NOME" | primer.type == "genomic" | primer.type == "CLEVER" | primer.type=="CrispRCas9PCR"){
                
                #per sequence id primer design overview plot
                
                #bed[ir,c("start","end","assembly","sequenceID","sequence.length","sequence.adress","sequence")]
                ipt<-NA
                bed.length<-nchar(as.character(bed[as.character(bed$sequenceID)==ibps,"sequence"]))
                sa<-seq.ana.all[[ibps]]#read.table(file="E:\\Science_Jil\\Projects\\SB_PrimerDesign\\#LiveDemo\\html_bis_withSNPanalysis\\sequences\\sequence.characterization_sequence1.txt",header=T,sep="\t",dec=".")
                sels<-results2[as.character(results2$sequence.id)==ibps,]
                
              }# if(not hairpin)
              
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
                
                if((primer.type == "hp_bisulfite" | primer.type == "hp_NOME" | primer.type == "hp_genomic" | primer.type == "hp_CLEVER") &&
                   hp.initial.input.type == "regions"){
                  
                  #trick to also show SNPs on the second part of the hairpin
                  tolo.temp<-tolo
                  tolo.temp$end.relative <- lnk.eim + as.numeric(as.character(tolo$start.relative)) - 1 
                  tolo.temp$start.relative <- lnk.eim + as.numeric(as.character(tolo$end.relative)) - 1
                  tolo<-rbind(tolo,tolo.temp)
                  
                }#if hp
                
                #debug output for tolo table
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
                
                if((primer.type == "bisulfite" | primer.type == "NOME" | primer.type == "genomic" | primer.type == "CLEVER") &&
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
                
                if(nrow(tolo[as.numeric(as.character(tolo$start.relative)) < 0,])>0 && !is.na(tolo[as.numeric(as.character(tolo$start.relative)) < 0,"start.relative"])){
                  tolo[as.numeric(as.character(tolo$start.relative)) < 0,"start.relative"]<-0
                }

                if(nrow(tolo[as.numeric(as.character(tolo$end.relative)) < 0,])>0 && !is.na(tolo[as.numeric(as.character(tolo$end.relative)) > bed.length,"end.relative"])){
                  tolo[as.numeric(as.character(tolo$end.relative)) > bed.length,"end.relative"]<-bed.length
                }
                
                #debug tolo
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
                
                if((primer.type == "bisulfite" | primer.type == "NOME" | primer.type == "genomic" | primer.type == "CLEVER") &&
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
      }# if create.graphics
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
    
    log("reset working directory...")
    setwd(start.wd)
    log("done.")
    
    log("end.")
    
    return(black)#output for primerwheel
    
  }#if(exists("results))
  
}#end of primer design pipeline