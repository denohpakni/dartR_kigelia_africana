install.packages("devtools")


install.packages("BiocManager")
BiocManager::install(c("SNPRelate", "qvalue"))

install.packages("dartR") 

install.packages("tidyverse")

library(devtools)
library(dartR) 

# test installation
gl.smearplot(testset.gl)


# make sure to get the STRUCTURE.exe software for your computer from https://web.stanford.edu/group/pritchardlab/home.html

