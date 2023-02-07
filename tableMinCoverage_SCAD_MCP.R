rm( list=ls() )   


#matrix the minimal coverage probabilities
#line number : a particular setting of n, design and model selection procedure
#column number: a particular target and constants
mMinCoverage = matrix(data=-1,nrow=12,ncol=8)  


cat("watershed_protected_intercept_p9_n20_SCAD_1 \n")
for (i in c(1,2,4,6,8)) {
  load(paste("watershed_protected_intercept_p9_n20_SCAD_1_withRes/res",i,".Rdata",sep=""))
  mMinCoverage[1,i] = minCovP
}

cat("watershed_protected_intercept_p9_n100_SCAD_1 \n")
for (i in c(1,2,4,6,8)) {
  load(paste("watershed_protected_intercept_p9_n100_SCAD_1_withRes/res",i,".Rdata",sep=""))
  mMinCoverage[2,i] = minCovP
}

cat("watershed_protected_intercept_p9_n20_MCP_1 \n")
for (i in c(2,4,5,6,8)) {
  load(paste("watershed_protected_intercept_p9_n20_MCP_1_withRes/res",i,".Rdata",sep=""))
  mMinCoverage[3,i] = minCovP
}

cat("watershed_protected_intercept_p9_n100_MCP_1 \n")
for (i in c(1,2,4,5,6,8)) {
  load(paste("watershed_protected_intercept_p9_n100_MCP_1_withRes/res",i,".Rdata",sep=""))
  mMinCoverage[4,i] = minCovP
}

cat("exchangeable_protected_intercept_p9_n20_SCAD_1 \n")
for (i in c(1,2,4,5,6,8)) {
  load(paste("exchangeable_protected_intercept_p9_n20_SCAD_1_withRes/res",i,".Rdata",sep=""))
  mMinCoverage[5,i] = minCovP
}

cat("exchangeable_protected_intercept_p9_n100_SCAD_1 \n")
for (i in c(1,2,4,5,6,8)) {
  load(paste("exchangeable_protected_intercept_p9_n100_SCAD_1_withRes/res",i,".Rdata",sep=""))
  mMinCoverage[6,i] = minCovP
}

cat("exchangeable_protected_intercept_p9_n20_MCP_1 \n")
for (i in c(1,2,4,5,6,8)) {
  load(paste("exchangeable_protected_intercept_p9_n20_MCP_1_withRes/res",i,".Rdata",sep=""))
  mMinCoverage[7,i] = minCovP
}

cat("exchangeable_protected_intercept_p9_n100_MCP_1 \n")
for (i in c(1,2,4,5,6,8)) {
  load(paste("exchangeable_protected_intercept_p9_n100_MCP_1_withRes/res",i,".Rdata",sep=""))
  mMinCoverage[8,i] = minCovP
}

cat("equicorrelated_protected_intercept_p9_n20_SCAD_1 \n")
for (i in c(1,2,4,5,6,8)) {
  load(paste("equicorrelated_protected_intercept_p9_n20_SCAD_1_withRes/res",i,".Rdata",sep=""))
  mMinCoverage[9,i] = minCovP
}

cat("equicorrelated_protected_intercept_p9_n100_SCAD_1 \n")
for (i in c(1,4,5,6)) {
  load(paste("equicorrelated_protected_intercept_p9_n100_SCAD_1_withRes/res",i,".Rdata",sep=""))
  mMinCoverage[10,i] = minCovP
}

cat("equicorrelated_protected_intercept_p9_n20_MCP_1 \n")
for (i in c(2,4,6,8)) {
  load(paste("equicorrelated_protected_intercept_p9_n20_MCP_1_withRes/res",i,".Rdata",sep=""))
  mMinCoverage[11,i] = minCovP
}

cat("equicorrelated_protected_intercept_p9_n100_MCP_1 \n")
for (i in c(4,5)) {
  load(paste("equicorrelated_protected_intercept_p9_n100_MCP_1_withRes/res",i,".Rdata",sep=""))
  mMinCoverage[12,i] = minCovP
}



#rounding to the nearest percentage
round(mMinCoverage,2)

