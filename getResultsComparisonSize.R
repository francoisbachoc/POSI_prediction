rm( list=ls() )   #equivalent de clear scilab
set.seed(1)
source("../functions.r")
source("functions_for_Taylor.R")

vn=c(20)
vp=c(10)

for (n in vn) {
  for (p in vp) {
    cat("################################################################### \n")
    cat("comparisonSize_n",n,"_p",p," \n",sep="")
    cat("################################################################### \n")
    load(paste("comparisonSize_n",n,"_p",p,sep=""))
    cat("median length with K1: ",median(vLengthK1[vContains==1])," \n",sep="")
    cat("median length with K3: ",median(vLengthK3[vContains==1])," \n",sep="")
    cat("median length with K4: ",median(vLengthK4[vContains==1])," \n",sep="")
    cat("median length with Taylor method: ",median(vLengthTaylor[vContains==1])," \n",sep="")
  }
}
  















