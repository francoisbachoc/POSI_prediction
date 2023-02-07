rm( list=ls() )   #equivalent de clear scilab
set.seed(1)
source("../functions.r")



#settings
sigma=1
n=100
p=10
alpha=0.05
I=1000   #I,J for numerical computatiom of the constants
J=1000
modelSelection=modelSelectionMCP
argModelSelection=0 # unused for MCP
nmc1=1000    #number of Monte Carlo steps and of simulated beta in the finding of the minimal coverage probability
nbeta1=1000
nmc2=10000
nbeta2=100
nbatch3=50 #nbatch MC computation of size nmc. The results are averaged
nmc3 = 10000


#Computation of the matrix Sigma from a textbook design
c <- sqrt( 0.8/(p-2) )
X0 <- diag(p-1)
X0[,1] <- c
X0[1,1] <- sqrt( 1-(p-2)*c^2)
Sigma = t(X0)%*%X0


#Sample of a random design matrix and predictor vector
cSigma=chol(Sigma)
X = matrix(nrow=n,ncol=p-1,data=rnorm(n=n*(p-1)))%*%cSigma
x0 = t(cSigma)%*%matrix(nrow=p-1,ncol=1,data=rnorm(n=p-1))

#add the intercept
Sigma = rbind(0,Sigma)
Sigma = cbind(0,Sigma)
Sigma[1,1]=1
X = as.matrix(cbind(1,X))
colnames(X) = c("1","2","3","4","5","6","7","8","9","10")
Xc = toCanonicalCoordinate(X)
x0 = c(1,x0)

#save of the settings and data
save(x0,X,Xc,Sigma,sigma,n,p,alpha,I,J,modelSelection,argModelSelection,nmc1,nbeta1,nmc2,nbeta2,nmc3,nbatch3,file="settings.Rdata")




















