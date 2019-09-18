processInput <-function(table.in,#filename.in = NULL, # direct path to the file with i) regions (.bed format) or ii) sequences (filepath) (character)
                                 input.type="regions",#"regions" (for a file with regions in .bed format) or "sequences" (for a ".txt" file with sequences) (fixed vocabulary) (character)
                                 path.out = NULL, # directory where results will be stored (directory) (character)
                                 analysis.id="PrimerAnalysis",# Do not use'SonderZeichen'
                                 primer.type="bisulfite",#"NOME" or "bisulfite" or "hp_bisulfite" or "hp_NOME" or "genomic" or "hp_genomic" or "CLEVER" or "hp_CLEVER" (fixed vocabulary) (character)
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
                        ){  
  
  #################### Check kind of input: BED file or Sequence file ################################
  if (input.type == "sequence") {
    
  }
  
  if (input.type == "regions"){
    
  }
  
  
}