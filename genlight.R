tuvitu1 <- gl.read.silicodart(
  filename="SilicoDArT_RAND.csv",
  ind.metafile="SNP_population.csv")


tuvitu2 <- gl.read.dart( 
  filename="SNP_RAND.csv", 
  ind.metafile="SNP_population.csv")

gl.save(tuvitu,file = "tmp.Rdata")


hetz <- gl.report.heterozygosity(tuvitu1)


write.csv(hetz,file = "./Results/tuvitu.csv")


# tidy Up!
rm(ibdReport)
