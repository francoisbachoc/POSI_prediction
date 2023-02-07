rm( list=ls() )   #equivalent de clear scilab
source("functions.r")

maxOfInnerProducsbis = function(x_0,X,I) {
  #x_0........ : vector of size p
  #X...........: the pXk design matrix (full rank, k <=p)   (we can call this function when X is already an extracted matrix)
  #I...........: number of Monte Carlo sample
  #return......: vector of I realizations of max_M ( tM' V ). V unif on Sp
  #
  p = dim(X)[1]
  mtM = alltbM(x_0,X) #matrix of all the tM over all submodels
  U = matrix(nrow = p,ncol=I,data=rnorm(p*I))  #matrix of I realizations of N(0,I_p)
  U = U / (matrix(data=1,nrow=p,ncol=1) %*% matrix(data=sqrt(colSums(U^2)),nrow=1,ncol=I)) #renormalization (unif on Sp)
  vC = matrix(data=1,nrow=1,ncol=I)
  for (i in 1:I)
  {
    vC[i] = max(abs(mtM%*%U[,i]))
  }
  return( vC )  #the vectors of the I max
}

Konebis = function(x_0,X,r,alpha,I=1000) {
  #x_0........ : vector of size p
  #X...........: the pXp design matrix (full rank)  
  #r...........: number of degrees of freedom in variance estimation
  #alpha.......: confidence level
  #I...........: number of Monte Carlo sample for evaluating proba of excedeence of max of inner products
  #return......: Constant K1, estimated by Monte Carlo
  #
  p = dim(X)[1]
  vC = maxOfInnerProducsbis(x_0,X,I) #I realizations of maximum of inner product
  fconfidence = function(K){mean(tailN(K/vC,p,r))-alpha} #MC evaluation of confidence level for a constant K
  uniroot(fconfidence,interval=c(1,100))$root #empirical quantile 1-alpha ofNC
}

#test of adequation with HL function
p=20
r=50
alpha=0.05
I = 1000
x_0 = matrix(nrow=p,ncol=1,data=runif(p))
X = matrix(nrow=p,ncol=p,data=runif(p^2))
ptm1 <- proc.time()
K1 = Konebis(x_0,X,r,alpha,I)
ptm2 = proc.time() 
ptm2-ptm1



