# load packages

library(devtools)
library(dartR) 


# manual explaining all the DArTR functions is here > https://cran.r-project.org/web/packages/dartR/dartR.pdf

# Convert data to genlight object

Scoreskig <- gl.read.silicodart(
  filename="SilicoDArT_selected.csv",
  ind.metafile="SNP_population.csv")
    

snpkig <- gl.read.dart( 
  filename="SNP_selected.csv", 
  ind.metafile="SNP_population.csv")

gl.save(snpkig,file = "tmp.Rdata")


######################## basic statistics ##################
# for each loci (Hs, Ho, Fis etc.)
# run ?basic.stats for details
basic.stats(snpkig)
basicStats <- gl.basic.stats(snpkig)

write.csv(basicStats,file = "./Results/basicStats.csv")

expectedPopHetz1 <- utils.het.pop(snpkig)


#Calculates the expected heterozygosities for each population in a genlight object
expectedPopHetz <- gl.test.heterozygosity(snpkig)


############### Expected Heterozygosity per loci #############
expectedHetz <- gl.He(snpkig)
print(expectedPopHetz)
write.csv(expectedHetz,file = "./Results/expectedHetz.csv")


############### observed Heterozygosity per loci #############

observedHetz <- gl.Ho(snpkig)
write.csv(observedHetz,file = "./Results/observedHetz.csv")

gl.report.heterozygosity(snpkig)
x <- snpkig
out <- utils.basic.stats(snpkig.gl, digits = 3)
out <- utils.het.pop(snpkig)

############################# Diversity Analysis ###########################

gl.report.maf(snpkig) # MAF for each locus for SNP data

gl.report.diversity(snpkig) #diversity indexes for SNPs
gl.report.heterozygosity(snpkig) #Estimates Heterozygosity

gl.diagnostics.hwe(snpkig)
gl.report.hwe(snpkig) # populations with less than 5 individuals are skipped


# The fixation index ( FST) is a measure of population differentiation due to genetic structure. 
gl.fst.pop(snpkig)

#### AMOVA ####
# AMOVA is used to detect whether or not there is significant population structure
AMOVA_snp <- gl.amova(snpkig)
AMOVA_snp

AMOVA_Scores <- gl.amova(Scoreskig)
AMOVA_Scores

# PCA Analysis
pcaScores <- gl.pcoa(
  Scoreskig,
  nfactors = 5,
  correction = NULL,
  mono.rm = TRUE,
  parallel = FALSE,
  n.cores = 16,
  plot.out = TRUE,
  save2tmp = FALSE,
  verbose = NULL
)

gl.pcoa.plot(pcaScores, Scoreskig)

gl.grm(snpkig) # mean probability of identity by state (IBS)


# Generates percentage allele frequencies by locus and population
gl.percent.freq(snpkig)


# Generate a geographical map
gl.map.interactive(snpkig)

###### Population STRUCTURE ##########################
gl.tree.nj(Scoreskig)
gl.tree.nj(snpkig)# nj tree to summarize genetic similarity among populations

## Run for SNP 
sr <- gl.run.structure(snpkig, k.range = 2:5, num.k.rep = 3, exec = './structure.exe')
ev <- gl.evanno(sr)
ev
qmat <- gl.plot.structure(sr, K=3)
qmat
head(qmat)
gl.map.structure(qmat, K=3, snpkig, scalex=1, scaley=0.5)

################### Hardy-Weinberg tests over loci and populations ############
gl.hwe.pop(snpkig)

################### Mantel test ############
mantelTest <- gl.ibd(snpkig)

print(mantelTest)
write.csv(mantelTest,file = "./Results/mantelTest.csv")


############## Report private alleles in one population compared with a second population ##############
gl.report.pa(snpkig)

gl.grm(snpkig)

########## Distances #########

gl.dist.pop(snpkig)
gl.dist.ind(snpkig)


########### Distance Matrices #############

gl.report.ld.map(snpkig)

snpdm <- gl.propShared(snpkig) #Calculates a similarity (distance) matrix
write.csv(snpdm,file = "./Results/snpdm.csv")


#### Reports ####


gl.report.rdepth(snpkig,save2tmp=TRUE)

gl.list.reports()

gl.print.reports(9)

write.csv()

utils.recalc.avgpic(snpkig, verbose= NULL)

gl.report.maf(snpkig)
gl.report.parent.offspring
gl.report.diversity(snpkig)
gl.report.hwe(snpkig)
