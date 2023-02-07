rm( list=ls() )   


vInfLambda = matrix(data=-1,nrow=12,ncol=1)  


cat("watershed_protected_intercept_p9_n20_SCAD_1 \n")
load("watershed_protected_intercept_p9_n20_SCAD_1_withRes/settings.Rdata")
vInfLambda[1] = min(svd(X)$d)

cat("watershed_protected_intercept_p9_n100_SCAD_1 \n")
load("watershed_protected_intercept_p9_n100_SCAD_1_withRes/settings.Rdata")
vInfLambda[2] = min(svd(X)$d)

cat("watershed_protected_intercept_p9_n20_SMCP_1 \n")
load("watershed_protected_intercept_p9_n20_MCP_1_withRes/settings.Rdata")
vInfLambda[3] = min(svd(X)$d)

cat("watershed_protected_intercept_p9_n100_MCP_1 \n")
load("watershed_protected_intercept_p9_n100_MCP_1_withRes/settings.Rdata")
vInfLambda[4] = min(svd(X)$d)

cat("exchangeable_protected_intercept_p9_n20_SCAD_1 \n")
load("exchangeable_protected_intercept_p9_n20_SCAD_1_withRes/settings.Rdata")
vInfLambda[5] = min(svd(X)$d)

cat("exchangeable_protected_intercept_p9_n100_SCAD_1 \n")
load("exchangeable_protected_intercept_p9_n100_SCAD_1_withRes/settings.Rdata")
vInfLambda[6] = min(svd(X)$d)

cat("exchangeable_protected_intercept_p9_n20_SMCP_1 \n")
load("exchangeable_protected_intercept_p9_n20_MCP_1_withRes/settings.Rdata")
vInfLambda[7] = min(svd(X)$d)

cat("exchangeable_protected_intercept_p9_n100_MCP_1 \n")
load("exchangeable_protected_intercept_p9_n100_MCP_1_withRes/settings.Rdata")
vInfLambda[8] = min(svd(X)$d)


cat("equicorrelated_protected_intercept_p9_n20_SCAD_1 \n")
load("equicorrelated_protected_intercept_p9_n20_SCAD_1_withRes/settings.Rdata")
vInfLambda[9] = min(svd(X)$d)

cat("equicorrelated_protected_intercept_p9_n100_SCAD_1 \n")
load("equicorrelated_protected_intercept_p9_n100_SCAD_1_withRes/settings.Rdata")
vInfLambda[10] = min(svd(X)$d)

cat("equicorrelated_protected_intercept_p9_n20_SMCP_1 \n")
load("equicorrelated_protected_intercept_p9_n20_MCP_1_withRes/settings.Rdata")
vInfLambda[11] = min(svd(X)$d)

cat("equicorrelated_protected_intercept_p9_n100_MCP_1 \n")
load("equicorrelated_protected_intercept_p9_n100_MCP_1_withRes/settings.Rdata")
vInfLambda[12] = min(svd(X)$d)
