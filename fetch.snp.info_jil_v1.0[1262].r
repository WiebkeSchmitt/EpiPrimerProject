

#more info here: https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
#returns a data.frame with requested snp info

fetch.snp.info = function(assembly = NULL,
                          chr = NULL,
                          start = NULL,
                          end = NULL,
                          attributes = c('refsnp_id','chrom_start','chrom_strand','allele','minor_allele','minor_allele_freq'),
                          ...){
  
  supported.assemblies = c("hg19","hg38","mm10")
  
  if(! assembly %in% supported.assemblies){
    
    stop(paste0("Assembly not supported. Currently supported assemblies are: ",supported.assemblies))
    
  }
  
  
  require(biomaRt)
  
  mart.dataset = switch(assembly, "hg19" = "hsapiens_snp", 
                        "hg38" = "hsapiens_snp",
                        "mm10" = "mmusculus_snp")
  
  mart.host = switch(assembly,
                     "hg19" = "grch37.ensembl.org", 
                     "hg38" = "www.ensembl.org",
                     "mm10" = "www.ensembl.org")
  
  #configure mart
  #snpmart = useEnsembl(biomart = "snp", dataset = database)
  
  
  snp.mart = useMart(biomart = "ENSEMBL_MART_SNP", 
                     host = mart.host, 
                     path = "/biomart/martservice", 
                     dataset = mart.dataset) #mouse:mmusculus_snp
  
  
  #call snp info from mart
  snp.results = getBM(attributes = attributes, 
                      filters = c('chr_name','start','end'), 
                      values = list(as.numeric(gsub("chr","",chr)),
                                    as.numeric(gsub("chr","",start)),
                                    as.numeric(gsub("chr","",end))), 
                      mart = snp.mart)
  
  return(snp.results)
  
}

