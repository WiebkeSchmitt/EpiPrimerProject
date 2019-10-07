library(methods)
library(Biostrings)
library(GenomicRanges)
library(BSgenome)

#S4 version of referenceSequence
setClass("referenceSequence",
  contains = "BSgenome"
)

#TODO: clean out redundant variables
#TODO: replace with bioConductor equivalents like Bstrings and GRanges.
 setClass("PrimerPair",
   representation(
     
     ## Name of the sequence
     sequence.id = "character",
     
     ## Name of this particular primer pair
     amplicon.id = "character",
     
     ## Amplicon length
     amplicon.length = "integer",
     
     ## Relative amplicon start
     amplicon.start.relative = "integer",
     
     ## Relative amplicon end
     amplicon.end.relative = "integer",
     
     ## Strand top/bottom
     dna.strand = "character",
     
     ## Number of GCs
     nGCs = "integer",
     
     ## Number of CGs
     nCGs = "integer",
     
     ## Sequence of the first primer
     primer1.sequence = "character",
     
     ## Sequence of the second primer
     primer2.sequence = "character",
     
     ## Length of first primer
     primer1.length = "integer",
     
     ## Length of second primer
     primer2.length = "integer",
     
     ## Melting temperature of first primer
     primer1.tm = "numeric",
     
     ## Melting temperature of second primer
     primer2.tm = "numeric",
     
     ## Difference in melting temperature
     tm.difference = "numeric",
     
     ## Relative start of first primer
    primer1.start.relative = "integer",
     
     ## Relative end of first primer
     primer1.end.relative = "integer",
     
     ## Relative start of second primer
     primer2.start.relative = "integer",
     
     ## Relative end of second primer
     primer2.end.relative = "integer",
     
     ## primer pair id
     primer.pair.id = "character",
     
     ## ID of first primer
     primer1.id = "character",
     
     ## ID of second primer
     primer2.id = "character",
     
     ## ID of first fragment
     fragment1.id = "character",
     
     ## ID of second fragment
     fragment2.id = "character",
   
     ## ID for both fragments
     fragment12.id = "character",
     
     ## ID for both fragments and primers
     fragment12.primer.ids = "character",
     
     ## GC content of first primer
     primer1.gc.content = "numeric",
     
     ## GC content of second primer
     primer2.gc.content = "numeric",
     
     ## GC content in amplicon
     amplicon.gc.content = "numeric",
     
     ## genomic GC content in amplicon
     amplicon.genomic.gc.content = "numeric",
     
     ## Number of C-to-T conversions in first primer??
     primer1.c2t.conversion = "integer",
     
     ## Number of G-to-A conversions in second primer??
     primer2.g2a.conversion = "integer",
     
     ## Self-alignment of first primer
     primer1.self.alignment = "integer",
     
    ## Self-alignment of second primer
     primer2.self.alignment = "integer",
     
     ## alignments between primers
     primer1.primer2.alignment = "integer",
   
    ## Amplicon sequence
    amplicon.sequence = "character",

    ## Genomic sequence of first primer
    primer1.sequence.genomic = "character",

    ## Genomic sequence fo second primer
    primer2.sequence.genomic = "character",

    ## Genomic amplicon sequence
    amplicon.sequence.genomic = "character",

    ## index
    index = "integer"
  ),
  prototype(
    sequence.id=NA_character_,
    amplicon.id=NA_character_,
    amplicon.length=NA_integer_,
    amplicon.start.relative=NA_integer_,
    amplicon.end.relative=NA_integer_,
    dna.strand=NA_character_,
    nGCs=NA_integer_,
    nCGs=NA_integer_,
    primer1.sequence=NA_character_,
    primer2.sequence=NA_character_,
    primer1.length=NA_integer_,
    primer2.length=NA_integer_,
    primer1.tm=NA_real_,
    primer2.tm=NA_real_,
    tm.difference=NA_real_,
    primer1.start.relative=NA_integer_,
    primer1.end.relative=NA_integer_,
    primer2.start.relative=NA_integer_,
    primer2.end.relative=NA_integer_,
    primer.pair.id=NA_character_,
    primer1.id=NA_character_,
    primer2.id=NA_character_,
    fragment1.id=NA_character_,
    fragment2.id=NA_character_,
    fragment12.id=NA_character_,
    fragment12.primer.ids=NA_character_,
    primer1.gc.content=NA_real_,
    primer2.gc.content=NA_real_,
    amplicon.gc.content=NA_real_,
    amplicon.genomic.gc.content=NA_real_,
    primer1.c2t.conversion=NA_integer_,
    primer2.g2a.conversion=NA_integer_,
    primer1.self.alignment=NA_integer_,
    primer2.self.alignment=NA_integer_,
    primer1.primer2.alignment=NA_integer_,
    amplicon.sequence=NA_character_,
    primer1.sequence.genomic=NA_character_,
    primer2.sequence.genomic=NA_character_,
    amplicon.sequence.genomic=NA_character_,
    index=NA_integer_
  )
)

setClass("bisulifitePrimer",
 contains = "PrimerPair",
 representation(
 )
)

setClass("NOMePrimer",
  contains = "PrimerPair",
  representation(
  )
)

setClass("genomicPrimer",
  contains = "PrimerPair",
  representation(
  )
)

setClass("CLEVERPrimer",
  contains = "PrimerPair",
  representation(
  )
)

as.data.frame.PrimerPair <- function(object) {
  return(
    data.frame(
      sequence.id=object@sequence.id,
      amplicon.id=object@amplicon.id,
      amplicon.length=object@amplicon.length,
      amplicon.start.relative=object@amplicon.start.relative,
      amplicon.end.relative=object@amplicon.end.relative,
      dna.strand=object@dna.strand,
      nGCs=object@nGCs,
      nCGs=object@nCGs,
      primer1.sequence=object@primer1.sequence,
      primer2.sequence=object@primer2.sequence,
      primer1.length=object@primer1.length,
      primer2.length=object@primer2.length,
      primer1.tm=object@primer1.tm,
      primer2.tm=object@primer2.tm,
      tm.difference=object@tm.difference,
      primer1.start.relative=object@primer1.start.relative,
      primer1.end.relative=object@primer1.end.relative,
      primer2.start.relative=object@primer2.start.relative,
      primer2.end.relative=object@primer2.end.relative,
      primer.pair.id=object@primer.pair.id,
      primer1.id=object@primer1.id,
      primer2.id=object@primer2.id,
      fragment1.id=object@fragment1.id,
      fragment2.id=object@fragment2.id,
      fragment12.id=object@fragment12.id,
      fragment12.primer.ids=object@fragment12.primer.ids,
      primer1.gc.content=object@primer1.gc.content,
      primer2.gc.content=object@primer2.gc.content,
      amplicon.gc.content=object@amplicon.gc.content,
      amplicon.genomic.gc.content=object@amplicon.genomic.gc.content,
      primer1.c2t.conversion=object@primer1.c2t.conversion,
      primer2.g2a.conversion=object@primer2.g2a.conversion,
      primer1.self.alignment=object@primer1.self.alignment,
      primer2.self.alignment=object@primer2.self.alignment,
      primer1.primer2.alignment=object@primer1.primer2.alignment,
      amplicon.sequence=object@amplicon.sequence,
      primer1.sequence.genomic=object@primer1.sequence.genomic,
      primer2.sequence.genomic=object@primer2.sequence.genomic,
      amplicon.sequence.genomic=object@amplicon.sequence.genomic,
      index=object@index))
}