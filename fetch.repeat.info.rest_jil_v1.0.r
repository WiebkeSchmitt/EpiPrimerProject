
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
  r <- GET(paste(server, ext, sep = ""))
  stop_for_status(r)
  
  s = content(r)
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
