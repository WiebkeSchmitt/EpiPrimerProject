

rtl.get.snp.info.by.region <- function(assembly = NULL,
                                   chr = NULL,
                                   start = NULL,
                                   end = NULL,
                                   ...){
  require(rtracklayer)
  
  session <- browserSession()
  genome(session) <- assembly
  #trackNames(session) ## list the track names
  trackID = switch(EXPR = assembly,
                   'mm9' = 'snp128',
                   'mm10' = 'snp142Common',
                   'hg18' = 'snp130',
                   'hg19' = 'snp151Common',
                   'hg38' = 'snp151Common')
  
  myranges = GRangesForUCSCGenome(assembly, chr,IRanges(start, end))
  
  query <- ucscTableQuery(session, trackID, range = myranges)
  
  #tableNames(query)
  tableID = switch(EXPR = assembly,
                   'mm9' = 'snp128',
                   'mm10' = 'snp142Common',
                   'hg18' = 'snp130',
                   'hg19' = 'snp151Common',
                   'hg38' = 'snp151Common')
  
  tableName (query) <- tableID
  
  getab = getTable(query)
  return(getab)
  
}

#example
#res = rtl.get.snp.info.by.region("hg19",chr = "chr10",start = 8095800, end = 8099900)
#res







