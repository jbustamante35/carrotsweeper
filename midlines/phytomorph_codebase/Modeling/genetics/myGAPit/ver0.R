source("http://www.bioconductor.org/biocLite.R") 
#biocLite("multtest")
#install.packages("gplots")
#install.packages("LDheatmap")
#install.packages("genetics")
#install.packages("scatterplot3d") #added on 3/26/15
library(multtest)
library("gplots")
library("LDheatmap")
library("scatterplot3d") #added on 3/26/15
library("genetics")
library("compiler") #this library is already installed in R
source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
workPath <- "/mnt/snapper/Takeshi/GAPIT_CHTC/Schnable/";
workPathLocal <- "/mnt/scratch2/GAPIT/GAPIT_CHTC/Schnable/";
#Step 1: Set working directory and import data
setwd(workPathLocal)
myY <- read.table("diallele_GWAS_matched.txt", head = TRUE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
    Y=myY,
    PCA.total=3,
    file.G="reremake_geno_13M_chr",
    file.Ext.G="hmp.txt",
    file.from=1,
    file.to=10,
    SNP.fraction=.6,
    #file.fragment=128,
    file.path=workPathLocal
  )
