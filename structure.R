# For testing with sample data from DArTR website

str_snp <- gl.read.dart( 
  filename="str_snp.csv", 
  ind.metafile="str_pop.csv")


sr <- gl.run.structure(str_snp, k.range = 2:5, num.k.rep = 3, exec = './structure.exe')
ev <- gl.evanno(sr)
ev
qmat <- gl.plot.structure(sr, K=3)
head(qmat)
gl.map.structure(qmat, K=3, str_snp, scalex=1, scaley=0.5)
