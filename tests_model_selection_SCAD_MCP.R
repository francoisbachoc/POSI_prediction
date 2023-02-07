rm( list=ls() )   #equivalent de clear scilab
source("functions.r")


#test that the SCAD function model selection returns the restricted least square estimator
p=3
n=10
beta = matrix(nrow=p,ncol=1,data=c(0,0,0))
X <- matrix(nrow=n,ncol=p,data=c(seq(from=1,to=1,length=n),runif(n*(p-1))))
U <- as.matrix(  rnorm(n=n,sd=0.1,mean=0))
Y = X%*%beta + U
model=modelSelectionSCAD(X,Y)$model
hatbeta=modelSelectionSCAD(X,Y)$hatbeta
XM = as.matrix(X[,model[model!=0]])
hatbeta2 = solve(t(XM)%*%XM)%*%t(XM)%*%Y


#test of conservative SCAD model selection
nmc=100
n=1000
p=4
beta = c(0,3,3,0)
nonZeroCoeff = c(1,2,3)
X <- matrix(nrow=n,ncol=p,data=c(seq(from=1,to=1,length=n),runif(n*(p-1))))
vMcorrect = seq(length=nmc)
for (imc in 1:nmc) {
  U <- as.matrix(  rnorm(n=n,sd=1,mean=0))
  Y = X%*%beta + U
  model=modelSelectionSCAD(X,Y,k=2)$model
  vMcorrect[imc] = TRUE
  for (j in nonZeroCoeff) {
    vMcorrect[imc] = is.element(j,model)
  }
}
mean(vMcorrect)


#test of consistent post-model-selection estimation of beta with SCAD
nmc=100
n=1000
p=4
beta = c(2,0,0,3)
X <- matrix(nrow=n,ncol=p,data=c(seq(from=1,to=1,length=n),runif(n*(p-1))))
vSquareError = seq(length=nmc)
for (imc in 1:nmc) {
  U <- as.matrix(  rnorm(n=n,sd=1,mean=0))
  Y = X%*%beta + U
  resModSec = modelSelectionSCAD(X,Y,k=2)
  hatbeta = resModSec$hatbeta
  model = resModSec$model
  hatbetafull = c(0,0,0,0)
  hatbetafull[model[model!=0]] = hatbeta
  vSquareError[imc] = sum((beta - hatbetafull)^2) 
}
mean(vSquareError)


#test that the MCP function model selection returns the restricted least square estimator
p=3
n=10
beta = matrix(nrow=p,ncol=1,data=c(1,0,2))
X <- matrix(nrow=n,ncol=p,data=c(seq(from=1,to=1,length=n),runif(n*(p-1))))
U <- as.matrix(  rnorm(n=n,sd=0.1,mean=0))
Y = X%*%beta + U
model=modelSelectionMCP(X,Y)$model
hatbeta=modelSelectionMCP(X,Y)$hatbeta
XM = as.matrix(X[,model[model!=0]])
hatbeta2 = solve(t(XM)%*%XM)%*%t(XM)%*%Y


#test of conservative MCP model selection
nmc=100
n=1000
p=4
beta = c(0,3,3,0)
nonZeroCoeff = c(1,2,3)
X <- matrix(nrow=n,ncol=p,data=c(seq(from=1,to=1,length=n),runif(n*(p-1))))
vMcorrect = seq(length=nmc)
for (imc in 1:nmc) {
  U <- as.matrix(  rnorm(n=n,sd=1,mean=0))
  Y = X%*%beta + U
  model=modelSelectionMCP(X,Y,k=2)$model
  vMcorrect[imc] = TRUE
  for (j in nonZeroCoeff) {
    vMcorrect[imc] = is.element(j,model)
  }
}
mean(vMcorrect)


#test of consistent post-model-selection estimation of beta with MCP
nmc=100
n=1000
p=4
beta = c(2,0,0,3)
X <- matrix(nrow=n,ncol=p,data=c(seq(from=1,to=1,length=n),runif(n*(p-1))))
vSquareError = seq(length=nmc)
for (imc in 1:nmc) {
  U <- as.matrix(  rnorm(n=n,sd=1,mean=0))
  Y = X%*%beta + U
  resModSec = modelSelectionMCP(X,Y,k=2)
  hatbeta = resModSec$hatbeta
  model = resModSec$model
  hatbetafull = c(0,0,0,0)
  hatbetafull[model[model!=0]] = hatbeta
  vSquareError[imc] = sum((beta - hatbetafull)^2) 
}
mean(vSquareError)

















