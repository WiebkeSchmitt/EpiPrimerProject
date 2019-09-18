source("C:\Users\Wiebk\Documents\epiprimer\R\Primerpair.r")

generalDesign<-function(seq,
                        seq.id,
                        nome,
                        tmg,
                        subfrags,
                        min.length.primer, 
                        max.length.primer, 
                        min.Tm.primer,
                        max.Tm.primer, 
                        max.Tm.difference.primer, 
                        low.complexity.primer.removal, 
                        max.bins.low.complexity,
                        remove.primers.with.n,
                        min.C2T.primer1,
                        min.G2A.primer2,
                        min.length.amplicon,
                        max.length.amplicon, 
                        min.number.gc.amplicon, 
                        min.number.cg.amplicon,
                        primer.align.binsize,
                        strand,
                        mode){
  
  
  out<-as.data.frame(new("PrimerPair"))
  
  #very powerfull way! no loops here!
  #it is much quicker to use this way though it creates some wrong primers.
  #these can be removed in a later step easily.
  
  ##### primer length ranges to build start and end points for subfragments... #####
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
  

  
  ####### calculate primer melting temps. #####
  print ("Calculate primer melting temperatures...")
  pr.tm<-sapply(X = pr.sq,calculate.Tm)
  print("Done.")
  
  
  ######## remove primers with inappropriate melting temps. #####
  print ("Remove primers with Tm too high or low...")
  f.good.tm<- pr.tm >= min.Tm.primer & pr.tm <= max.Tm.primer
  pr.sq<-pr.sq[f.good.tm]
  pr.tm<-pr.tm[f.good.tm]
  print("Done.")
  
  ###### establish primer combinations that cover enough sites of interest #####
  print ("Remove primer combinations that do not cover enough sites of interest...")
  
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
  
  ######## keep only if enough features are present between the two subfragments.... #####
  print ("Remove primer combinations with too few sites of interest...")
  pr.pairs<-pr.pairs[as.numeric(gsub("fragment","",pr.pairs$fragment2.id))-
                       as.numeric(gsub("fragment","",pr.pairs$fragment1.id))>= min.number.gc.amplicon+min.number.cg.amplicon,]
  print("Done.")
  
  if(nrow(pr.pairs)<1){
    print("Could not identify a primer pair that includes sites of interest.")
    print(paste(date()))
    return(out)
  }
  
  #### remove primer combinations that do not fullfil tm.delta criteria ####
  print ("Remove primer combinations with too high Tm difference...")
  pr.pairs2<-pr.pairs[pr.pairs$tm.difference <= max.Tm.difference.primer,]
  print("Done.")
  
  if(nrow(pr.pairs2)<1){
    print("Could not identify a primer pair that includes sites of interest.")
    print(paste(date()))
    return(out)
  }
  
  
  ##### Remove primers that have "N"s #####
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
  
  
  ###### Match relative position of individual primers #####
  print("Calculate relative primer Positions...")
  len.1<-unlist(sapply(X = pr.pairs2$primer1.sequence,regexpr,nome,useBytes=TRUE,fixed=TRUE))
  len.2<-unlist(sapply(X = pr.pairs2$primer2.sequence,regexpr,nome,useBytes=TRUE,fixed=TRUE))
  
  len.2e<-len.2+nchar(names(len.2))-1
  len<-len.2e-len.1 #v1.7: ~7 sec
  
  pr.pairs2$amplicon.length<-as.numeric(len)
  pr.pairs2$amplicon.start.relative<-len.1
  pr.pairs2$amplicon.end.relative<-len.2e
  print("Done.")
  
  ##### Select amplicons by AmpliconLength #####
  print("Filter amplicons by min/max length...")
  amp.sel<-pr.pairs2[pr.pairs2$amplicon.length>=min.length.amplicon & 
                       pr.pairs2$amplicon.length<=max.length.amplicon,]
  
  if (nrow(amp.sel)==0){
    print("No primers found. Try to change parameters [Amplicon Length].")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel), " primer pairs survived...",sep=""))
  
  ###### fetch amplicon sequences for each primer combination. ####
  amp.sel$amplicon.sequence<-mapply(substr,nome,
                                    start = amp.sel$amplicon.start.relative,
                                    stop = amp.sel$amplicon.end.relative)
  amp.sel2 <- amp.sel #for historic reason, actually nothing happens here.
  amp.sel2$fragment12.id<-paste(amp.sel2$fragment1.id,"_",amp.sel2$fragment2.id,sep="")
  amp.sel2$fragment12.primer.ids<-paste(amp.sel2$primer1.id,"_",amp.sel2$primer2.id,sep="")
  
  ###### select best tm primer combination for individual amplicons.... #####
  print("Select optimal primer pair combination for Tm...")
  
  #VERY POWERFULL FUNCTION:
  lvl.frgs<-levels(factor(amp.sel2$fragment12.id))
  mins<-sapply(X = lvl.frgs, function(x) {y=which.min(amp.sel2[amp.sel2$fragment12.id==x,
                                                                    "tm.difference"])
  return(amp.sel2[amp.sel2$fragment12.id==x,"fragment12.primer.ids"][y])}
  )
  amp.sel2<-amp.sel2[match(mins,amp.sel2$fragment12.primer.ids),]
  
  ###### check amplicons for gc/cg contents #####
  print("Scan amplicons for information content (GpCs/CpGs)...")
  
  amp.sel2$nCGs<-sapply(X=amp.sel2$amplicon.sequence, function(x){
    cgs<-as.numeric(unlist(gregexpr("cg",text = x,ignore.case=TRUE)))
    ncgs<-length(cgs)
    if(cgs[1]== -1){
      ncgs<-0
    }
    return(ncgs)
  })
  
  amp.sel2$nGCs<-sapply(X=amp.sel2$amplicon.sequence, function(x){
    gcs<-as.numeric(unlist(gregexpr("gc",text = x,ignore.case=TRUE)))
    ngcs<-length(gcs)
    if(gcs[1]== -1){
      ngcs<-0
    }
    return(ngcs)
  })
  
  ###### filter for gc and gc contents #####
  print("Filter amplicons for CpGs content...")
  amp.sel3<-amp.sel2[amp.sel2$nCGs>=min.number.cg.amplicon,]
  
  if(nrow(amp.sel3)==0){
    print("No Amplicon survived...Try to lower number of analyzable CpGs.")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel3)," amplicons survived...")) 
  
  ###### check for numbers of CpG in amplicons #####
  print("Filter amplicons for GpC content...")
  
  amp.sel4<-amp.sel3[amp.sel3$nGCs>=min.number.gc.amplicon,]
  
  if(nrow(amp.sel4)==0){
    print("No Amplicon survived...Try to lower number of analyzable GpCs.")
    print(paste(date()))
    return(out)
  }
  
  print(paste(nrow(amp.sel4)," amplicons survived..."))  
  
  ###### LOW COMPLEXITY PRIMER REMOVAL #####
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
  
  ###### Scan Primer:Primer Alignments ####
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
  
  ###### Scan Primer:Primer Alignments ####
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
  
  ##### Scan Primer:Primer Alignments ####
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
  
  ######### Remove amplicons that have primer with high risk for primer-dimers... ##############
  
  print("Remove amplicons that have primer with high risk for primer-dimers...")
  selection<-amp.sel4[amp.sel4$primer1.self.alignment== -1 &
                        amp.sel4$primer2.self.alignment== -1 &
                        amp.sel4$primer1.primer2.alignment== -1,]
  print("Done.")
  
  
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
  
  ##### finalize out object #####
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
                                  "amplicon.gc.content","amplicon.genomic.gc.content",
                                  "primer1.c2t.conversion","primer2.g2a.conversion",
                                  "primer1.self.alignment",
                                  "primer2.self.alignment",
                                  "primer1.primer2.alignment",
                                  "amplicon.sequence",
                                  "primer1.sequence.genomic","primer2.sequence.genomic",
                                  "amplicon.sequence.genomic",
                                  "index")])
  
  ################################# Done! ############################################      
  
  return(out)
}#generalDesign


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
  
  #initiate sequence
  seq<-as.character(sequence)
  seq.id<-as.character(sequence.id)
  
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
  
  return(generalDesign(seq,
                       seq.id,
                       nome,
                       tmg,
                       subfrags,
                       min.length.primer, 
                       max.length.primer, 
                       min.Tm.primer,
                       max.Tm.primer, 
                       max.Tm.difference.primer, 
                       low.complexity.primer.removal, 
                       max.bins.low.complexity,
                       remove.primers.with.n,
                       min.C2T.primer1,
                       min.G2A.primer2,
                       min.length.amplicon,
                       max.length.amplicon, 
                       min.number.gc.amplicon, 
                       min.number.cg.amplicon,
                       primer.align.binsize,
                       strand,
                       mode))
}#nome.primer.design


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
  source("Primerpair.r")
  
  version.id<-3.0
  
  #initiate sequence
  seq<-as.character(sequence)
  seq.id<-as.character(sequence.id)
  
  ########################################################################
  
  #perform bisulfite conversion of genomic sequence
  nome<-bisulfite.conversion(seq,strand = strand)
  nome<-tolower(nome)
  
  #########################################################################
  
  #split the nome sequence by all "C"s.
  c.split<-unlist(strsplit(x=nome,split="c"))
  names(c.split)<-paste("fragment",seq_along(c.split),sep="")
  n.fragments<-length(c.split)
  
  ########################################################################
  
  #Remove splitted sequences that are shorter than the minimum length of the desired oligos.
  #c.split.length<-lapply(c.split,nchar)
  long.enough<-c.split[nchar(c.split)>=min.length.primer]
  
  ########################################################################
  
  #Prepare NEGATIVE RESULT OUTPUT
  
  out<-as.data.frame(new("PrimerPair"))
    # data.frame(sequence.id=NA,amplicon.id=NA,amplicon.length=NA,amplicon.start.relative=NA,amplicon.end.relative=NA,
    #               dna.strand=NA,nGCs=NA,nCGs=NA,primer1.sequence=NA,primer2.sequence=NA, 
    #               primer1.length=NA,primer2.length=NA,primer1.tm=NA,primer2.tm=NA,tm.difference=NA,
    #               primer1.start.relative=NA,primer1.end.relative=NA,primer2.start.relative=NA,
    #               primer2.end.relative=NA,primer.pair.id=NA,primer1.id=NA,primer2.id=NA,fragment1.id=NA,
    #               fragment2.id=NA,fragment12.id=NA,fragment12.primer.ids=NA,primer1.gc.content=NA,
    #               primer2.gc.content=NA,amplicon.gc.content=NA,amplicon.genomic.gc.content=NA,
    #               primer1.c2t.conversion=NA,primer2.g2a.conversion=NA,primer1.self.alignment=NA,
    #               primer2.self.alignment=NA,primer1.primer2.alignment=NA,amplicon.sequence=NA,
    #               primer1.sequence.genomic=NA,primer2.sequence.genomic=NA,amplicon.sequence.genomic=NA,
    #               index=NA)
  
  
  
  if (length(long.enough)<2){
    print("No primer design possible. Try to reduce min.length.primer.")
    print(paste(date()))
    return(out)
  }
  
 
  
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
  tmg$primer.start.starts<-1
  tmg$primer.start.ends<-tmg$fragment.length-min.length.primer
  tmg<-tmg[tmg$primer.start.ends>0,]
  tmg$primer.end.starts<-1+min.length.primer-1
  
  subfrags<-names(tm.good)[names(tm.good) %in% tmg$fragment.id]
  
  return(generalDesign(seq,
                       seq.id,
                       nome,
                       tmg,
                       subfrags,
                       min.length.primer, 
                       max.length.primer, 
                       min.Tm.primer,
                       max.Tm.primer, 
                       max.Tm.difference.primer, 
                       low.complexity.primer.removal, 
                       max.bins.low.complexity,
                       remove.primers.with.n,
                       min.C2T.primer1,
                       min.G2A.primer2,
                       min.length.amplicon,
                       max.length.amplicon, 
                       min.number.gc.amplicon, 
                       min.number.cg.amplicon,
                       primer.align.binsize,
                       strand,
                       mode))
}#bisulfite.primer.design


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
  
  #initiate sequence
  seq<-as.character(sequence)
  seq.id<-as.character(sequence.id)
  
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
  tmg$primer.start.starts<-1
  tmg$primer.start.ends<-tmg$fragment.length-min.length.primer
  tmg<-tmg[tmg$primer.start.ends>0,]
  tmg$primer.end.starts<-1+min.length.primer-1
  
  subfrags<-names(tm.good)[names(tm.good) %in% tmg$fragment.id]
  
  return(generalDesign(seq,
                       seq.id,
                       nome,
                       tmg,
                       subfrags,
                       min.length.primer, 
                       max.length.primer, 
                       min.Tm.primer,
                       max.Tm.primer, 
                       max.Tm.difference.primer, 
                       low.complexity.primer.removal, 
                       max.bins.low.complexity,
                       remove.primers.with.n,
                       min.C2T.primer1,
                       min.G2A.primer2,
                       min.length.amplicon,
                       max.length.amplicon, 
                       min.number.gc.amplicon, 
                       min.number.cg.amplicon,
                       primer.align.binsize,
                       strand,
                       mode))
}#genomic.primer.design


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
  
  #initiate sequence
  seq<-as.character(sequence)
  seq.id<-as.character(sequence.id)
  
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
  tmg$primer.start.starts<-1
  tmg$primer.start.ends<-tmg$fragment.length-min.length.primer
  tmg<-tmg[tmg$primer.start.ends>0,]
  tmg$primer.end.starts<-1+min.length.primer-1
  
  subfrags<-names(tm.good)[names(tm.good) %in% tmg$fragment.id]
  
  return(generalDesign(seq,
                       seq.id,
                       nome,
                       tmg,
                       subfrags,
                       min.length.primer, 
                       max.length.primer, 
                       min.Tm.primer,
                       max.Tm.primer, 
                       max.Tm.difference.primer, 
                       low.complexity.primer.removal, 
                       max.bins.low.complexity,
                       remove.primers.with.n,
                       min.C2T.primer1,
                       min.G2A.primer2,
                       min.length.amplicon,
                       max.length.amplicon, 
                       min.number.gc.amplicon, 
                       min.number.cg.amplicon,
                       primer.align.binsize,
                       strand,
                       mode))
}#CLEVER.primer.design