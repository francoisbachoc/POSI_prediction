rm( list=ls() )   #equivalent de clear scilab
set.seed(1)
source("../functions.r")

####################################
#Settings
####################################

vmModelSel = c(  "SCAD","MCP")
vTarget = c(1,2)  #1 BetaMn, 2 BetaM
vK = c(0,1,3,4) #Knaive, K1,K3,K4
vDesign = c("watershed","exchangeable","equicorrelated")

I=10^5
nbatch=10^3
nSample = 20
periodMessage=1
alpha=0.05


####################################
#Script generations
####################################

#save settings
save(I,nbatch,nSample,periodMessage,alpha,file="settings.Rdata")

dquote="\""
for (i in 1:length(vName)) {
  for (ir in 1:length(vr)) {
    vLine = c("rm( list=ls() )")
    vLine = c(vLine,"set.seed(1)")
    vLine = c(vLine,paste("source(",dquote,"../functions.r",dquote,")",sep=""))
    vLine = c(vLine,"load(\"settings.Rdata\")")
    vLine = c(vLine,paste("X = as.matrix(read.table(",dquote,"../",vName[i],".dat",dquote,"))",sep=""))
    vLine = c(vLine,paste("r = ",vr[ir],sep=""))
    vLine = c(vLine,"vK = sampleKoneUnstructured(X,r,alpha,I,nbatch,nSample,periodMessage)")
    vLine = c(vLine,paste("save(vK, file = ",dquote,vName[i],"_r",vr[ir],".Rdata",dquote,")",sep=""))
    fileConn<-file(paste("run",(2*(i-1)+ir),".R",sep=""))
    writeLines(vLine, fileConn)
    close(fileConn)
  }
}















