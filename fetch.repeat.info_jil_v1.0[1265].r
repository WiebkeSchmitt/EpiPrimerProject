

#returns a annotation for repeats

#http://bioconductor.org/packages/release/bioc/vignettes/rtracklayer/inst/doc/rtracklayer.R

#returns a data.frame with repeat info with 17 columns: 
#"bin"       "swScore"   "milliDiv"  "milliDel"  "milliIns"  "genoName"  "genoStart" "genoEnd"   "genoLeft"  "strand"    "repName"   "repClass" 
#"repFamily" "repStart"  "repEnd"    "repLeft"   "id"

fetch.repeat.info= function(assembly = NULL, #'hg19','hg38', 'mm9','mm10'
                             chr = NULL, #'chr1' or '1'
                             start = NULL, #'12345678'
                             end = NULL,   #'12345679'
                             ...){
  require (rtracklayer)
  mySession = browserSession("UCSC")
  genome(mySession) <- assembly
  myrange <- GRanges(paste0("chr",gsub("chr","",chr)), IRanges(as.numeric(start),as.numeric(end)))
  tbl.rmsk <- getTable(
              ucscTableQuery(mySession, track="rmsk", 
              range=myrange, table="rmsk"))
 
  return(tbl.rmsk)
}


#example
#fetch.repeat.info(assembly = "mm10", chr = "chr1", start = 100000000, end = 100010000)

