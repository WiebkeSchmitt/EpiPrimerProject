

#retrieves SNP meta informations from REST database, requires SNP ID (rs#)

fetch.snp.stats.rest = function(id = NULL, # e.g. "rs539335951"
                          assembly = "hg38", #"hg19","hg38" or "mm10"
                          ...){

  require(httr)
  require(jsonlite)
  require(xml2)
  
  server <- switch(EXPR = assembly, 
                   "hg19" = "http://grch37.rest.ensembl.org",
                   "hg38" = "http://rest.ensembl.org",
                   "mm10" = "http://rest.ensembl.org")
  
  species = switch(EXPR = assembly, 
                  "hg19" = "homo_sapiens",
                  "hg38" = "homo_sapiens",
                  "mm10" = "mus_musculus")
  
  ext <- paste0('/variation/',species)
  
  r <- POST(paste(server, ext, sep = ''), content_type('application/json'), accept('application/json'), body = paste0('{ "ids" : ["',id,'" ] }'))
  
  stop_for_status(r)
  
 results = data.frame(t(sapply(content(r),c)))
  
  results$MAF = as.numeric(as.character(results$MAF))
  results$most_severe_consequence = paste0(unlist(results$most_severe_consequence),collapse=",")
  results$source = paste0(unlist(results$source),collapse=",")
  results$evidence = paste0(unlist(results$evidence),collapse=",")
  results$ancestral_allele = paste0(unlist(results$ancestral_allele),collapse=",")
  results$name = paste0(unlist(results$name),collapse=",")
  results$minor_allele = paste0(unlist(results$minor_allele),collapse=",")
  results$var_class = paste0(unlist(results$var_class),collapse=",")
  results$mappings = paste0(unlist(results$mappings),collapse=",")
  results$ambiguity = paste0(unlist(results$ambiguity),collapse=",")
  results$synonyms = paste0(unlist(results$synonyms),collapse=",")
  
  results = results[,order(colnames(results))]
  
  return(results)
   
}

#
#example
#ids = c("rs1305355659","rs539335951")
#sa = as.data.frame(t(sapply(ids,fetch.snp.stats.rest,assembly = "hg19)))
#
#
#
#










