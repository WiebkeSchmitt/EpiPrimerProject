


##mimiq-biq tool
#Rbiq<-function(filename.biq){
#}




####################################################################################################################

# sequences (can be multiple) need to be character and have names.
# reference (only one) needs to be character and have name.
# alignment.type: one of 'global', 'local', 'global-local' or 'local-global'.
# mode 'bisulfite' does not give penalties on c:t mismatches in the alignments.
# subsitutionMatrix: optionally you can deliver a matrix for alignment penalties.

#idea: prefilter reads by length? would lead to speed-up
#idea: paralellization of alignments?

pairwiseSequenceAlignment.biqFormat<-function(sequences.path, #path to reads file
                                    reference.path, #path to reference file
                                    mode = "bisulfite", 
                                    alignment.type="global-local",
                                    subsitutionMatrix=NULL, #should it be added to GUI ?
									minimum.alignment.score=NULL, #should be positive
                                    minimum.percentage.identity=NULL,# value 0-100, should be >50
                                      # reference.sequence=NA, #string
                                      reference.sequence.id, #optional name 
                                      sample.id, #optional name
                                      feature="c",
                                      a.value="x",
                                      c.value="1",
                                      g.value="x",
                                      t.value="0",
                                      n.value="x"
									
									
									
									
									
									
									
									){

require(Biostrings)
require(seqinr)
sequences <- sequences.path
reference <- reference.path
reference.sequence <- paste0(reference)



  seqs<-sapply(sequences,toupper)
  ref<-toupper(reference)
  
#set up nucleotide substitution matrix
  if(is.null(subsitutionMatrix)){
    #mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
    
    mat <- data.frame(A=c(1,-3,-3,-3,-2), 
                      C=c(-3,1,-3,-3,-2),
                      G=c(-3,-3,1,-3,-2),
                      T=c(-3,-3,-3,1,-2),
                      N=c(-2,-2,-2,-2,-2))
    rownames(mat) <-c("A","C","G","T","N")
    mat<-as.matrix(mat)
  
    if(mode=="bisulfite"){
      mat[2,4] <- 1
    }
    
  }
  
  if(!is.null(subsitutionMatrix)){
    mat<-subsitutionMatrix
  }


#perform pairwise alignments
  psa<-sapply(seqs, function(x) {
       pairwiseAlignment(pattern = ref, 
                         subject=x, 
                         substitutionMatrix=mat, 
                         scoreOnly=FALSE, 
                         type=alignment.type)
   }
  )
  
  
#alignment score
  score.all<-sapply(psa,score)
  
  if(!is.null(minimum.alignment.score)){
    
    psa<-psa[score.all>= minimum.alignment.score]
    
  }
  
  #percentage of identity  
  perc.identity<-sapply(psa,pid)
  
  if(!is.null(minimum.percentage.identity)){
    
    psa<-psa[perc.identity >= minimum.percentage.identity]
    
  }
  
  refs<-sapply(psa,pattern)
  rds<-sapply(psa,subject)
  
  get.cs<-sapply(refs,function(x)gregexpr(feature,x,ignore.case=TRUE))
  
  massive<-sapply(1:length(rds), function(x) c2s(substring(text=rds[[x]],
                                                           first=as.numeric(get.cs[[x]]), 
                                                           last=as.numeric(get.cs[[x]]))))
  
  massive<-gsub("A",a.value,gsub("T",t.value,gsub("C",c.value,gsub("G",g.value,
           gsub("N",n.value,gsub("-","x",massive))))))
  
  res<-data.frame(ID=names(psa),
                  Alignment.score=as.numeric(sapply(psa,score)),
                  Sequence.identity=as.numeric(sapply(psa,pid)),
                  Methylation.pattern=massive)
  
  nc.master.ref<-as.numeric(nchar(gsub("[A|G|T|a|g|t]","",reference.sequence)))
  
  #remove reads with more/less features than master ref.
  res<-res[nchar(as.character(res$Methylation.pattern))==nc.master.ref,]
  
  #produce more info
  res$Mean.methylation<-as.numeric(sapply(gsub("[0|x]","",res$Methylation.pattern),nchar)/nc.master.ref)
  res$Missing.sites<-as.numeric(sapply(gsub("[0|1]","",res$Methylation.pattern),nchar))
  res$Missing.sites.percentage<-res$Missing.sites/nc.master.ref
  res$Conversion.rate<-NA
  res$Reference.sequence<-reference.sequence.id
  res$Sample.name<-sample.id
  
  save_table <- write.table(res,"C:\\Users\\Wiebk\\OneDrive\\Masterthesis\\R\\epiprimer\\analysis\\result_biq.txt",row.names=FALSE,sep="\t")
  
  return(res)
  
  }
  