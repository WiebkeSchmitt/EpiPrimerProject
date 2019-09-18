

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
  
  r <- GET(paste(server, ext, sep = ""), content_type("text/plain"))
  
  stop_for_status(r)
  return(content(r))
  
}


#example
# fetch.dna.sequence(assembly = "hg19",
#                    chr="chr1",
#                    start="100000000",
#                    end="100001000")


