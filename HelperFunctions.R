source("generalDesign.R")

sequence.characterization<-function(sequence,
                                    feature1="cg",
                                    feature1.subpos=1,#1st position = 1, 2nd position =2,...
                                    feature2="gc",
                                    feature2.subpos=2){#1st position = 1, 2nd position =2,...
  myseq<-tolower(sequence)
  feature1<-tolower(feature1)
  feature2<-tolower(feature2)
  
  ###############################all pos splits###########################################################################
  s0<-strsplit(x = myseq,split = "")
  
  
  ############all feature1
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
  
  ############all feature2
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
  
  ###########feature.overlapps
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
  
  ############results as data.frame
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
  
  if(!is.na(ol1[1])){
    df[ol1,paste(feature1,"_overlapps_",feature2,sep="")]<-TRUE
  }
  if(!is.na(ol2[1])){
    df[ol2,paste(feature2,"_overlapps_",feature1,sep="")]<-TRUE
  }
  
  df[,paste(feature1,"_count",sep="")]<-NA
  df[subpos1,paste(feature1,"_count",sep="")]<-names(subpos1)
  df[,paste(feature2,"_count",sep="")]<-NA
  df[subpos2,paste(feature2,"_count",sep="")]<-names(subpos2)
  
  df[is.na(df)]<-""
  
  return(df)
}

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

#write a sequence in .fasta format
write.fasta<-function(sequence="Example",sequence.id="my.fasta",filename=paste(getwd(),"/my.fasta.FASTA",sep="")){
  cat(paste('>',sequence.id, '\n', sequence,sep=""),file=filename)
}

#read.fasta file
read.my.fasta<-function(filename){
  require("seqinr")
  fasta<-read.fasta(file = filename, 
                    seqtype = "DNA", as.string = TRUE, forceDNAtolower = TRUE,
                    set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = FALSE)
  return(fasta)  
}

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

#Export some text to a file...
#Warning: If file already exists it will be overwritten!
#Optional: add=TRUE, opens an existing file and adds the text without overwriting original contents.

write.text<-function(text="BLAH!",file="text.txt",add=FALSE){
  cat(text,file=file,append=add)
}

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
}#function hairpinizer

##############repeat info via Enseml API#########################
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
  
  for (i in 1:length(ext)){
    r <- GET(paste(server, ext[i], sep = ""))
  }
  
  stop_for_status(r)
  
  s = content(r)
  s2 = lapply(s,function(x) {tdf = data.frame(c1 = names(unlist(x)), c2 = unlist(x))})
  s3 = sapply(1:length(s2), function(x) { s2[x][[1]] = as.data.frame(s2[x][[1]][c("seq_region_name","start","end","strand","assembly_name","feature_type","description"),"c2"])})
  s4 = as.data.frame(matrix(unlist(s3),nrow = length(s3), ncol = 7, byrow = TRUE))
  colnames(s4) = c("chr","start","end","strand","assembly","feature_type","description")
  
  return(s4)
  
}

##############DNA annotation info via Enseml API#########################
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
  
  for (i in 1:length(ext)){
    r <- GET(paste(server, ext[i], sep = ""), content_type("text/plain"))
  }
  
  #r <- GET(paste(server, ext, sep = ""), content_type("text/plain"))
  
  stop_for_status(r)
  return(content(r))
  
}

##############SNP info info via Enseml API#########################
#get SNP info for a genomic intervall 
# will call ensembl rest API (http://mar2017.rest.ensembl.org/)

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
  
  r <- GET(paste(server, ext[1], sep=""))
  
  s = content(r)
  s2 = lapply(s,function(x) {tdf = data.frame(c1 = names(unlist(x)), c2 = unlist(x))})
  s3 = sapply(1:length(s2), function(x) { s2[x][[1]] = as.data.frame(s2[x][[1]][c("seq_region_name","start","end","strand","assembly_name","id","feature_type","consequence_type","source"),"c2"])})
  s4 = as.data.frame(matrix(unlist(s3),nrow = length(s3), ncol = 9, byrow = TRUE))
  colnames(s4) = c("chr","start","end","strand","assembly","rs_id","feature_type","consequence_type","source")
  return(s4)
  
}

#this function is deprecated and should not be used anymore
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




#get SNP info for a genomic intervall 
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
  r <- GET(paste(server, ext, sep = ""))
  
  stop_for_status(r)
  
  s = content(r)
  s2 = lapply(s,function(x) {tdf = data.frame(c1 = names(unlist(x)), c2 = unlist(x))})
  s3 = sapply(1:length(s2), function(x) {s2[x][[1]] = as.data.frame(s2[x][[1]][c("seq_region_name","start","end",
                                                                                 "strand","assembly_name","id",
                                                                                 "external_name","feature_type",
                                                                                 "biotype","version","description",
                                                                                 "source"),"c2"])})
  s4 = as.data.frame(matrix(unlist(s3),nrow = length(s3), ncol = 12, byrow = TRUE))
  colnames(s4) = c("chr","start","end","strand","assembly","id","gene_name","feature_type","biotype",
                   "version","description","source")
  
  return(s4)
  
}