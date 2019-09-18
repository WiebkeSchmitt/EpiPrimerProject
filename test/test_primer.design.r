library(testthat)
load("test_primer.design.bisulfite_GATA3_prom.rds")

test_that(
  "bisulfite design(GATA3_prom",
  {
    ref<-readRDS("test_primer.design.bisulfite_GATA3_prom.rds")
    sink("/dev/null")
    test<-bisulfite.primer.design(
      sequence="gttctttctgtccgtctacactgagcgtactcggggaatgagttagagccagtctcttcctcccctccccccttctcatccctcactgttgccactcaagtcaaaagcacacattgattacaaatattaggtctggaaagggcagctgcaacagctgaagcgtgttcactctgggggcttgagagcgcagaaggctcgggaaagaggtgacaatgacaacaaaattgacgcggacgctccagtcaaaggcatctcccctttatccgatgactcaccctcttaggaagtcggcccgagaggcaaatctcaaaataccttgacatgaaacattttgtttttctgatcaatttaacgcgcacgtttccccacatcgatgcgctctcccaaacaccctgcattagatcctaataatgatccatgcgtgcctattttttaaaagtctgaaaaagaaaattctgcccatcgaaatgaacttcatgaatggggcaggctggctgcaccgggacggaatcgtccacccgacccgaatgaattggcaggagccgcggccacatttaaagggccagagcgcgcgttccctcccgtccgcccccaagccccgcgggcctcgcccaccctgcccgccgcccctccgccggcggccgccctctgcggcgcccctttccggtcagtggaggggcgggaggaggggcgggggtgcgcggggcggggggagaagtcctggagcgggtttgggttgcagtttccttgtgccggggatcctgtcccctactcgccagcgccaggctcctcccccccggcgcggatgacactagaacctccttaagttgcgtcgcgccacagctgtctgcgaacactgagctgcctggcgccgtcttgatactttcagaaagaatgcattccctgtaaaaaaaaaaaaaaaatactgagagagggagagagagagagaagaagagagagagacggagggagagcgagacagagcgagcaacgcaatctgaccgagcaggtcgtacgccgccgcctcctcctcctctctgctcttcgctacccaggttggtactggtgactttttttttttttaagtttgattttttgcccccaaccacttgggaggacctaaatcaattttaaaaactcaactctcctcttttggaggttttctaggggctgagaggacggtcccgggaccggtgtccccgagggagggacttgccctccaagtcgtaacagtcagccctgggacttgccctccaagttgctcagccagccccggctcccgcgagccgggctgcagggacgtccccgagagccctgcgggctccgcggccgtgtccccgcgctcccgtgcgggtctcgggtgcgctgggcgggcgggcggcgcgaggggaggttgtgccactccagcaactcaggggctcatccaggtctcccattctctcccttgcaggtgacccgaggagggactccgcctccgagcggctgaggaccccggtgcagaggagcctggctcgcagaattgcagagtcgtcgcccctttttacaacctggtcccgttttattctgccgtacccagtttttggatttttgtcttccccttcttctctttgctaaacgacccctccaagataatttttaaaaaaccttctcctttgctcacctttgcttcccagccttcccatccccccaccgaaagcaaatcattcaacgacccccgaccctccgacggcaggagccccccgacctcccaggcggaccgccctccctccccgcgcgcgggttccgggcccggcgagagggcgcgagcacagccgaggccatggaggtgacggcggaccagccgcgctgggtgagccaccaccaccccgccgtgctcaacgggcagcacccggacacgcaccacccgggcctcagccactcctacatggacgcggcgcagtacccgctgccggaggaggtggatgtgctttttaacatcgacggtcaaggcaaccacgtcccgccctactacggaaactcggtcagggccacggtgcagaggtaccctccgacccaccacggtgagtgcgcccggggtgccggggctcccgccggccgcttcagccgtcccggctcggggaggtcgggagggacctgagggcggggagaggtcaagcgaaagcccccatctgccgttcctggttcatttacaaaaaaattgg",
      sequence.id="GATA3_prom",
      min.Tm.primer=48,
      max.Tm.primer=60,
      min.number.gc.amplicon=0,
      min.number.cg.amplicon=5,
      max.Tm.difference.primer=2,
      primer.align.binsize=12,
      min.length.primer=23,
      max.length.primer=34,
      low.complexity.primer.removal=TRUE,
      max.bins.low.complexity=7,
      remove.primers.with.n=TRUE,
      min.C2T.primer1=3,
      min.G2A.primer2=3,
      min.length.amplicon=150,
      max.length.amplicon=500,
      strand="top",
      mode="fast"
    )
    sink()
    expect_true(all.equal(ref,test))
  }
)

