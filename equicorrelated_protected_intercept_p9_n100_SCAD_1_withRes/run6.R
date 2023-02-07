rm( list=ls() )   #equivalent de clear scilab


####################################################
#loading of the settings
####################################################
load("settings.Rdata")
source("../functions.r")
set.seed(1)



######################################################
#execution
####################################################

#execution time
tStart=Sys.time()


#preparation of the constants
cat("########### \n")
cat("preparation of the constants \n")
cat("########### \n")
r=n-p
K = Kone(x0,Xc,r,alpha,I)

#choice of the target
target=2 #BetaM

#minCoverage probability
res=minCoverageProbability(nbeta1,nbeta2,nbatch3,nmc1,nmc2,nmc3,x0,X,Xc,sigma,Sigma,K,target,modelSelection,argModelSelection)
mbeta1 = res$mbeta1
mbeta2 = res$mbeta2
beta3 = res$beta3
vCov1 = res$vCov1
vCov2 = res$vCov2
vCov3 = res$vCov3
vHLCI1 = res$vHLCI1
vHLCI2 = res$vHLCI2
vHLCI3 = res$vHLCI3
vModSize1 = res$vModSize1
vModSize2 = res$vModSize2
vModSize3 = res$vModSize3
minCovP = res$minCovP


#the total execution time
tExec=Sys.time() - tStart



##############################################################
#Save of the result
##############################################################
save(tExec,mbeta1,mbeta2,beta3,vCov1,vCov2,vCov3,vHLCI1,vHLCI2,vHLCI3,vModSize1,vModSize2,vModSize3,minCovP, file = "res6.Rdata" )















