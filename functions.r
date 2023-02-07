library(lars)
library(selectiveInference)
library(ncvreg)

subsets=function(n,r,v=1:n){
  # subsects() gets all subsets of size r from a set of size n
  if(r<=0) NULL else
    if(r>=n) v[1:n] else
      rbind(cbind(v[1],subsets(n-1,r-1,v[-1])),subsets(n-1,r,v[-1]))
}

alltbM = function(x_0,X) {
  #x_0........ : vector of size p
  #X...........: the pXk design matrix (full rank, k <=p)   (we can call this function when X is already an extracted matrix)
  #return......: the matrix of size (2^k-1) X p  of all the vectors bar{t}_M (normalized). 
  #
  p = dim(X)[1]
  k = dim(X)[2]
  mtM = matrix(data=0,nrow=2^k-1,ncol=p)
  index = 1 # index of the next tM to be calculated
  for (m.sz in 1:k) {   #loop on size of models
    if (m.sz<k){ #obtaining the matrix of all the models of a given size
      models.of.sz <- subsets(k,m.sz)
    }else{
      models.of.sz <- matrix(1:k,1,k)
    }
    for (i in 1:dim(models.of.sz)[1]) {   #loop on models of a given size
      model = models.of.sz[i,]
      x_0M = x_0[model]
      XM = as.matrix(X[,model])
      mtM[index,] = t(x_0M) %*% solve( t(XM)%*%XM ) %*% t(XM)  #computation of t_M for a given M
      index = index+1
    }
  }
  vNorm = sqrt(rowSums(mtM^2)) #normalization of the non-zero vectors
  mtM[vNorm!=0,] = mtM[vNorm!=0,] / (matrix(data=vNorm[vNorm!=0],nrow=length(vNorm[vNorm!=0]),ncol=1)%*%matrix(data=1,nrow=1,ncol=p)) 
  return(mtM)
}

maxOfInnerProducs = function(x_0,X,I) {
  #x_0........ : vector of size p
  #X...........: the pXk design matrix (full rank, k <=p)   (we can call this function when X is already an extracted matrix)
  #I...........: number of Monte Carlo sample
  #return......: vector of I realizations of max_M ( tM' V ). V unif on Sp
  #
  p = dim(X)[1]
  mtM = alltbM(x_0,X) #matrix of all the tM over all submodels
  U = matrix(nrow = p,ncol=I,data=rnorm(p*I))  #matrix of I realizations of N(0,I_p)
  U = U / (matrix(data=1,nrow=p,ncol=1) %*% matrix(data=sqrt(colSums(U^2)),nrow=1,ncol=I)) #renormalization (unif on Sp)
  return( apply(abs(mtM%*%U),2,max) )  #the vectors of the I max
}

Kone = function(x_0,X,r,alpha,I=1000) {
  #x_0........ : vector of size p
  #X...........: the pXp design matrix (full rank)  
  #r...........: number of degrees of freedom in variance estimation
  #alpha.......: confidence level
  #I...........: number of Monte Carlo sample for evaluating proba of excedeence of max of inner products
  #return......: Constant K1, estimated by Monte Carlo
  #
  p = dim(X)[1]
  vC = maxOfInnerProducs(x_0,X,I) #I realizations of maximum of inner product
  fconfidence = function(K){mean(tailN(K/vC,p,r))-alpha} #MC evaluation of confidence level for a constant K
  uniroot(fconfidence,interval=c(1,100))$root #empirical quantile 1-alpha ofNC
}

Ktwo = function(x0,M,Xc,r,alpha,nx01,nx02,nbatch3,I1,I2,I3) {
  #xo,M,Xc,r,alpha: cf function Kthree
  #maximize the constant K1 with the Monte Carlo 3 step approach of Leeb, Poetscher, Ewald 2013.
  #nx01,nx02,nx03,I1,I2,I3: similar to function minCoverageProbability
  p = dim(X)[2]
  compM = seq(from=1,to=p,length=p)
  compM=compM[-M]
  XccompM = Xc[,compM]
  
  #matrices of beta and vector of corresponding coverage probabilities
  x01=matrix(nrow=1,ncol=p,data=x0)
  mx02 = matrix(nrow=nx02,ncol=1,data=1)%*%matrix(nrow=1,ncol=p,data=x0)
  vK2 = seq(from=-1,to=-1,length=nx02)
  
  #step 1: evaluate nxo1 K1 for x01 and keep ONLINE (memory issues) the best x01 to store then in mx02
  cat("########### \n")
  cat("Step 1 \n")
  cat("########### \n")
  for (ix01 in 1:nx01) {
    x01[,compM] = mRandomx0M(XccompM,1)
    cat("x0 number ",ix01,"\n")
    cat(Sys.time(),"\n")
    K1 = Kone(x01,Xc,r,alpha,I1)
    if (K1 > min(vK2))  {#if we improve on the nxo2 current values
      mx02[which.min(vK2),] = x01
      vK2[which.min(vK2)] = K1
    }
  }
  
  
  #step 2
  cat("########### \n")
  cat("Step 2 \n")
  cat("########### \n")
  for (ix02 in 1:nx02) {
    cat("x0 number ",ix02,"\n")
    cat(Sys.time(),"\n")
    x02 = mx02[ix02,1:p]
    vK2[ix02] = Kone(x02,Xc,r,alpha,I2)
  }
  
  #pick the best beta for step 3
  ind2 = sort(vK2,decreasing=TRUE,index.return=TRUE)$ix
  x03 = mx02[ind2[1],]
  
  #step 3
  cat("########### \n")
  cat("Step 3 \n")
  cat("########### \n")
  vK3 = seq(length=nbatch3)
  for (ix03 in 1:nbatch3) {
    cat("Monte Carlo batch number ",ix03,"\n")
    cat(Sys.time(),"\n")
    vK3[ix03] = Kone(x03,Xc,r,alpha,I3)
  }
  
  #the max K1
  maxK1 = mean(vK3)
  
  #return of the results
  list(mx02=mx02,x03=x03,vK2=vK2,vK3=vK3,maxK1=maxK1)
}

mRandomx0M = function(XM,m) {
  #Sample m random vectors x0Mnew (matrix m times M) that are Gaussian with covariance matrix t(XM)*XM
  XM = as.matrix(XM)
  pM = dim(XM)[2]
  Sigma = (1/pM)*(t(XM)%*%XM)
  cSigma=chol(Sigma)
  matrix(nrow=m,ncol=pM,data=rnorm(n=m*pM))%*%cSigma
}

Kfour = function(p,r,alpha,I=1000) {
  #p...........: number of variables
  #r...........: number of degrees of freedom in variance estimation
  #alpha.......: confidence level
  #I...........:number of quantiles for numerical evaluation of the mean value
  #return......: Constant K4, estimated by Monte Carlo
  #
  vC = sqrt(qbeta(p=seq(from=0,to=1/(2^p-1),length=I),shape1=1/2,shape2=(p-1)/2, lower.tail=FALSE)) #vector of quantiles of Beta distribution
  fconfidence = function(K){mean(tailN(K/vC,p,r))-alpha} #MC evaluation of confidence level for a constant K
  uniroot(fconfidence,interval=c(1,100))$root 
}

KthreeInf = function(x_0,M,X,r,alpha,I=1000,J=1000)  {
  #x_0.........: vector of size p
  #M...........: vector of selected variable indices, identifying the submodel
  #X...........: the pXp design matrix (full rank)  
  #r...........: number of degrees of freedom in variance estimation
  #alpha.......: confidence level
  #I...........: number of Monte Carlo samples for evaluating proba of exceedence of max of inner products
  #J...........: number of quantiles for the numerical evaluation of the mean value
  #return......: Constant K3, estimated by Monte Carlo
  #
  p = dim(X)[1]
  k = length(M)
  x0M = x_0[M]
  XM = as.matrix(X[,M])
  #
  #computations of I realizations of max of inner products
  vC = maxOfInnerProducs(x0M,XM,I) 
  #
  #binary search for finding mstar
  fUnionBound = function(m){mean(vC>m) + upperTailB(m,p,k) } #union bound giving tail function of C(x_{0M})
  mstar = uniroot(function(m){fUnionBound(m)-1},interval=c(0,1))$root
  
  #J binary searches for finding quantiles of of the beta part
  vM = seq(length=J-1)
  for (j in 1:(J-1)) {
    vM[j] = uniroot(function(m){upperTailB(mstar,p,k)*(j/J) - upperTailB(m,p,k)},interval=c(0,1))$root   
  }
  #
  #Binary search for calculating K3  
  fconfidence = function(K) {
    res = 1-tailN(K/mstar,p,r)
    res = res - (1/I) *sum( (vC > mstar)*(tailN(K/vC,p,r)-tailN(K/mstar,p,r)*seq(from=1,to=1,length=I) )  )
    res = res + upperTailB(mstar,p,k) *(1/J)*sum(tailN(K/mstar,p,r)*seq(from=1,to=1,length=J-1)-tailN(K/vM,p,r))
    res - (1-alpha)
  } #MC evaluation of confidence level for a constant K
  uniroot(fconfidence,interval=c(1,100))$root 
}

KthreeSup = function(x_0,M,X,r,alpha,I=1000,J=1000)  {
  #x_0........ : vector of size p
  #M...........: vector of selected variable indices, identifying the submodel
  #X...........: the pXp design matrix (full rank)  
  #r...........: number of degrees of freedom in variance estimation
  #alpha.......: confidence level
  #I..........: number of Monte Carlo samples for evaluating proba of exceedence of max of inner products
  #J..........: number of quantiles for the numerical evaluation of the mean value
  #return......: Constant K3, estimated by Monte Carlo
  #
  p = dim(X)[1]
  k = length(M)
  x0M = x_0[M]
  XM = as.matrix(X[,M])
  #
  #computations of I realizations of max of inner products
  vC = maxOfInnerProducs(x0M,XM,I) 
  #
  #binary search for finding mstar
  fUnionBound = function(m){mean(vC>m) + upperTailB(m,p,k)} #union bound giving tail function of C(x_{0M})
  mstar = uniroot(function(m){fUnionBound(m)-1},interval=c(0,1))$root
  
  #J binary searches for finding quantiles of of the beta part
  vM = seq(length=J)
  vM[1] = 1
  for (j in 2:J) {
    vM[j] = uniroot(function(m){upperTailB(mstar,p,k)*((j-1)/J) - upperTailB(m,p,k)},interval=c(0,1))$root   
  }
  #
  #Binary search for calculating K3  
  fconfidence = function(K) {
    res = 1-tailN(K/mstar,p,r)
    res = res - (1/I) *sum( (vC > mstar)*(tailN(K/vC,p,r)-tailN(K/mstar,p,r)*seq(from=1,to=1,length=I) )  )
    res = res + upperTailB(mstar,p,k) *(1/J)*sum(tailN(K/mstar,p,r)*seq(from=1,to=1,length=J)-tailN(K/vM,p,r))
    res - (1-alpha)
  } #MC evaluation of confidence level for a constant K
  uniroot(fconfidence,interval=c(1,100))$root 
}

KthreeOld = function(x_0,M,X,r,alpha,I1=1000,I2=1000)  {
  #x_0........ : vector of size p
  #M...........: vector of selected variable indices, identifying the submodel
  #X...........: the pXp design matrix (full rank)  
  #r...........: number of degrees of freedom in variance estimation
  #alpha.......: confidence level
  #I1..........: number of Monte Carlo samples for evaluating proba of exceedence of max of inner products
  #I2..........: number of quantiles for the numerical evaluation of the mean value
  #return......: Constant K3, estimated by Monte Carlo
  #
  p = dim(X)[1]
  k = length(M)
  x0M = x_0[M]
  XM = as.matrix(X[,M])
  #
  #computations of I1 realizations of max of inner products
  vM = maxOfInnerProducs(x0M,XM,I1) 
  #
  #I2 binary searches for finding quantiles of C
  vC = seq(length=I2)
  fUnionBound = function(c){mean(vM>c) + (2^p-2^k)*pbeta(q=c^2,shape1=1/2,shape2=(p-1)/2,lower.tail=FALSE)} #union bound giving tail function of C(x_{0M})
  for (i2 in 1:I2) {
    vC[i2] = uniroot(function(c){fUnionBound(c)-(i2/I2)},interval=c(0,1))$root   #computation of quantile i2/I2 of C(x_{0M})
  }
  #
  #Binary search for calculating K3  
  fconfidence = function(K){mean(tailN(K/vC,p,r))-alpha} #MC evaluation of confidence level for a constant K
  uniroot(fconfidence,interval=c(1,100))$root 
}

Kfive = function(p,r,alpha) {
  #p...........: number of variables
  #r...........: number of degrees of freedom in variance estimation
  #alpha.......: confidence level
  #return......: Constant K5 (scheffe constant)
  #
  sqrt(p*qf(alpha,df1=p,df2=r,lower.tail=FALSE))
}
  
tailN = function(q,p,r)  {
  #q...........: a quantile
  #p...........: parameter of distribution N in the paper
  #r...........: parameter of distribution N in the paper
  #return......: tail function of N at q (cf paper)
  #
  1-pf(q^2/p,p,r)
}

upperTailB = function(q,p,k)  {
  #q...........: a quantile
  #p...........: dimension on sphere
  #k...........: in the (2^p - 2^k) (cf paper)
  #return......: upper bound of the tail function of the beta part in the paper
  #
  (2^p - 2^k) * pbeta(q^2,shape1=1/2,shape2=(p-1)/2,lower.tail=FALSE)
}
  
modelSelectionAIC = function(X,Y,k=2) {
  #X...........: the design matrix (NOT in canonical coordinate)
  #Y...........: observation vector
  #k...........: penalty factor (k=2 default value for AIC)
  #return......: list composed of:
  #model.......: vector defining the model of the form [index_1,...,index_M,0...,0] (total size p) (can be the empty model)
  #hatbeta.....: the post model selection least square estimator (matrix of dim 0 if empty model)
  #
  #preparative
  n = dim(X)[1]
  p = dim(X)[2]
  X  = data.frame(X)
  for (i in 1:p) {
    colnames(X)[i] = paste("V",i,sep="") #variable indexed by numbers only
  }
  #
  #AIC procedure with protected intercept
  f0<-lm(Y~.-1, data=X)
  f1 <- step(f0, trace=0,scope=list(lower="Y~V1"),k=k)
  #
  #getting selected model and restricted least square estimator
  model = colnames(f1$qr$qr)  # indices of selected variables
  for (i in 1:length(model)) {  #loop for removing the quotation marks
    model[i] = substr(model[i],start=2,stop=nchar(model[i]))
  }
  model = as.numeric(c(model,seq(from=0,to=0,length=p-length(model)))) #completion with 0s for having length p
  hatbeta = as.vector(f1$coefficients) #least square estimator
  list(model=model,hatbeta=hatbeta)
}

modelSelectionLASSO = function(X,Y,k=0) {
  #preparative
  n = dim(X)[1]
  p = dim(X)[2]
  X  = data.frame(X)
  #
  #removing the 1st column and taking residuals
  Xr = X[,-1]
  yr = residuals(lm(Y~X[,1]-1))      
  #
  #LASSO procedure
  lasso.cv <- cv.lars(Xr,yr, type="lasso",plot.it=F, intercept=F, se=F)
  lasso <- lars(Xr,yr,type="lasso",intercept=F)
  betahat.lasso <- coef(lasso, mode="fraction",s=lasso.cv$index[which.min(lasso.cv$cv)])
  #
  #model corresponds to non-zero coefficients
  M <- which(betahat.lasso!=0)
  M <- c(1,M+1)
  M = c(M,seq(from=0,to=0,length=p-length(M)))
  #
  #returning least square estimator
  XM = as.matrix(X[,M[M!=0]])
  hatbeta = solve(t(XM)%*%XM)%*%t(XM)%*%Y
  list(model=M,hatbeta=hatbeta)
}

modelSelectionSCAD = function(X,Y,k=0) {
  #preparative
  n = dim(X)[1]
  p = dim(X)[2]
  X  = data.frame(X)
  #
  #removing the 1st column 
  Xr = X[,-1]
  yr = Y
  #
  #SCAD procedure (always add and protect an intercept)
  cvfit <- cv.ncvreg(Xr, yr,family="gaussian",
                     penalty = "SCAD")
  fit <- cvfit$fit
  betahat.SCAD <- fit$beta[,cvfit$min]
  #
  #model corresponds to non-zero coefficients
  M <- which(betahat.SCAD!=0)
  M = c(M,seq(from=0,to=0,length=p-length(M)))
  #
  #returning least square estimator
  XM = as.matrix(X[,M[M!=0]])
  hatbeta = solve(t(XM)%*%XM)%*%t(XM)%*%Y
  list(model=M,hatbeta=hatbeta)
}

modelSelectionMCP = function(X,Y,k=0) {
  #preparative
  n = dim(X)[1]
  p = dim(X)[2]
  X  = data.frame(X)
  #
  #removing the 1st column 
  Xr = X[,-1]
  yr = Y
  #
  #MCP procedure (always add and protect an intercept)
  cvfit <- cv.ncvreg(Xr, yr,family="gaussian",
                     penalty = "MCP")
  fit <- cvfit$fit
  betahat.MCP <- fit$beta[,cvfit$min]
  #
  #model corresponds to non-zero coefficients
  M <- which(betahat.MCP!=0)
  M = c(M,seq(from=0,to=0,length=p-length(M)))
  #
  #returning least square estimator
  XM = as.matrix(X[,M[M!=0]])
  hatbeta = solve(t(XM)%*%XM)%*%t(XM)%*%Y
  list(model=M,hatbeta=hatbeta)
}

matrixKthreeOld = function(x_0,X,r,alpha,I=1000,J=1000) {
  #input arguments........:cf function KthreeInf
  #return.................:matrix (2^p-1)X(p+1). each line is a model M and the corresponding K3(M)
  p = dim(X)[1]
  mK3 = matrix(data=0,nrow=2^p-1,ncol=p+1)
  index = 1 # index of the next K3 to be calculated
  for (m.sz in 1:p) {   #loop on size of models
    if (m.sz<p){ #obtaining the matrix of all the models of a given size
      models.of.sz <- subsets(p,m.sz)
    }else{
      models.of.sz <- matrix(1:p,1,p)
    }
    for (i in 1:dim(models.of.sz)[1]) {   #loop on models of a given size
      model = models.of.sz[i,]
      mK3[index,1:m.sz] = model   #fill in the model in the corresponding line
      x_0M = x_0[model]
      XM = as.matrix(X[,model])
      mK3[index,p+1] = KthreeOld(x_0,model,X,r,alpha,I,J) #computing K3 for the current model
      index = index+1
    }
  }
  mK3
}
 
matrixKthreeInf = function(x_0,X,r,alpha,I=1000,J=1000) {
  #input arguments........:cf function KthreeInf
  #return.................:matrix (2^p-1)X(p+1). each line is a model M and the corresponding K3(M)
  p = dim(X)[1]
  mK3 = matrix(data=0,nrow=2^p-1,ncol=p+1)
  index = 1 # index of the next K3 to be calculated
  for (m.sz in 1:p) {   #loop on size of models
    if (m.sz<p){ #obtaining the matrix of all the models of a given size
      models.of.sz <- subsets(p,m.sz)
    }else{
      models.of.sz <- matrix(1:p,1,p)
    }
    for (i in 1:dim(models.of.sz)[1]) {   #loop on models of a given size
      model = models.of.sz[i,]
      mK3[index,1:m.sz] = model   #fill in the model in the corresponding line
      x_0M = x_0[model]
      XM = as.matrix(X[,model])
      mK3[index,p+1] = KthreeInf(x_0,model,X,r,alpha,I,J) #computing K3 for the current model
      index = index+1
    }
  }
  mK3
}

matrixKthreeSup = function(x_0,X,r,alpha,I=1000,J=1000) {
  #input arguments........:cf function KthreeSup
  #return.................:matrix (2^p-1)X(p+1). each line is a model M and the corresponding K3(M)
  p = dim(X)[1]
  mK3 = matrix(data=0,nrow=2^p-1,ncol=p+1)
  index = 1 # index of the next K3 to be calculated
  for (m.sz in 1:p) {   #loop on size of models
    if (m.sz<p){ #obtaining the matrix of all the models of a given size
      models.of.sz <- subsets(p,m.sz)
    }else{
      models.of.sz <- matrix(1:p,1,p)
    }
    for (i in 1:dim(models.of.sz)[1]) {   #loop on models of a given size
      model = models.of.sz[i,]
      mK3[index,1:m.sz] = model   #fill in the model in the corresponding line
      x_0M = x_0[model]
      XM = as.matrix(X[,model])
      mK3[index,p+1] = KthreeSup(x_0,model,X,r,alpha,I,J) #computing K3 for the current model
      index = index+1
    }
  }
  mK3
}
  
coverageProbability = function(x0,X,Sigma,beta,sigma,I,fK,modelSelection,argModelSelection) {
  #X..................: the design matrix (NOT in canonical coordinate)
  #Sigma..............: covariance matrix of the predictors
  #beta...............: true vector of coefficients
  #sigma..............: true standard deviation of Gaussian errors
  #I..................: number of Monte Carlo samples
  #Knaive.............: naive constant for confidence intervals
  #fK.................: either a constant of the form K1, or a matrix of constant for all submodels of the form mK3
  #modelSelection.....: the model selection function (AIC, BIC, LASSO)
  #argModelSelection..: For LASSO a conventional 0, for AIC, BIC, the factor k
  #return.............: 3 averages for:
  #coverageBetaMn......: coverage of design-dependent coverage target
  #coverageBetaM.......: coverage of design-independent coverage target
  #HLCI................: half length of CI
  #modSize.............: size of the selected model
  #
  #
  #
  #preparative
  n = dim(X)[1]
  p = dim(X)[2]
  Xc = toCanonicalCoordinate(X)
  mU = matrix(nrow=n,ncol=I,data=rnorm(n=n*I,mean=0,sd=sigma))  #sample of the I noise vectors
  vCoverageBetaM = seq(length=I,from=1,to=1)  #column index: number of constant K. row index: number of MC simulation. For betaM
  vCoverageBetaMn = seq(length=I,from=1,to=1)  #For betaMn   IMPORTANT: we fill with 1 and it is not changed if empty model is selected
  vHLCI = seq(length=I,from=1,to=1)
  vModSize = seq(length=I,from=1,to=1)
  #
  #Monte Carlo procedure
  for (imc in 1:I) {  #loop of the Monte Carlo sample
    Y = as.matrix(X)%*%matrix(nrow=p,ncol=1,data=beta) + mU[,imc]   #realization of the observation vector
    modelSelectionOutcome = modelSelection(X,Y,argModelSelection) # the selected model and the estimated beta
    model = modelSelectionOutcome$model
    hatbeta = modelSelectionOutcome$hatbeta 
    if (length(hatbeta)>0) {   #if the selected model is not empty, then there is something to cover
      if ( length(fK)==1 ) {  #if we have a constant constant (e.g. K1)
        K=fK
      } else { #if the constant depends on the selected model
        K = selectFromModel(fK,model)   #choice of the right K3, as a function of the selected model
      }
      errorBetaM = sum(x0[model]*hatbeta) - sum(x0[model]*betaM(Sigma,beta,model)) #difference of inner products
      errorBetaMn = sum(x0[model]*hatbeta) - sum(x0[model]*betaMn(Xc,beta,model))  #difference of inner products
      hatSigma = sqrt((1/(n-p))*sum( (Y - X%*%solve((t(X)%*%X))%*%t(X)%*%Y)^2 ))   #estimation of the variance in full linear model
      normtM = sqrt( sum(tM(x0,Xc,model[model!=0])^2) )  #computation of norm of t_M for the selected model
      HLCI_ = hatSigma*normtM*K
      vCoverageBetaMn[imc] = abs(errorBetaMn) < HLCI_
      vCoverageBetaM[imc] = abs(errorBetaM) < HLCI_
      vHLCI[imc] = HLCI_
      vModSize[imc] = length(hatbeta)
    }
  }
  list(coverageBetaMn=mean(vCoverageBetaMn),coverageBetaM=mean(vCoverageBetaM),HLCI = mean(vHLCI),modSize=mean(vModSize))
}

tM = function(x0,X,model) {
  #x0..........: vector of predictors
  #X...........: design matrix
  #model.......: vector of size length(M) composed by the indices in the model
  #return......: vector tM of the paper
  #
  n = dim(X)[1]
  p = dim(X)[2]
  x0M = matrix(x0[model],nrow=length(model),ncol=1)  
  XM = X[,model]
  XM%*% solve( t(XM)%*%XM ) %*% x0M
}

betaMn = function(X,beta,model) {
  #X...........: design matrix
  #beta........: the true regression vector
  #model.......: vector of size length(M) composed by the indices in the model
  #return......: beta_Mn, the vector giving the design-dependant coverage target
  #
  model=model[model!=0]
  n = dim(X)[1]
  p = dim(X)[2]
  XM = X[,model]
  XcompM = as.matrix(X[,-model])
  betaM_ = matrix(nrow=length(model),ncol=1,data=beta[model])
  betaCompM = matrix(nrow=(p-length(model)),ncol=1,data=beta[-model])
  betaM_ + solve(t(XM)%*%XM)%*%t(XM)%*%XcompM%*%betaCompM #formula of the paper
}
  
betaM = function(Sigma,beta,model) {
  #Sigma.......: the covariance matrix of the predictors
  #beta........: the true regression vector
  #model.......: vector of size length(M) composed by the indices in the model
  #return......: beta_M, the vector giving the distribution-dependant coverage target
  #
  model=model[model!=0]
  p = dim(Sigma)[1]
  SigmaMM = Sigma[model,model]
  SigmaMcompM = Sigma[model,-model]
  betaM_ = matrix(nrow=length(model),ncol=1,data=beta[model])
  betaCompM = matrix(nrow=(p-length(model)),ncol=1,data=beta[-model])
  betaM_ + solve(SigmaMM)%*%SigmaMcompM%*%betaCompM   #formula of the paper
}


selectFromModel = function(M,model) {
  #M..........: matrix of size lX(p+1), corresponding to l models among 1...p and l associated numbers
  #model......: vector of size p defining a model
  #return.....: the number in the matrix corresponding to model. 0 if the model is not in the matrix
  #
  p = dim(M)[2] -1
  l = dim(M)[1]
  index = which.min( rowSums((M[,1:p]- (matrix(data=1,nrow=l,ncol=1)%*%matrix(data=model,nrow=1,ncol=p)))^2) )
  res = 0
  if (norm(M[index,1:p]-matrix(nrow=1,ncol=p,data=model)) < 0.001  ) {
    res = M[index,p+1]
  }
  res
}



mRandomBeta <- function(X,scale,m) {
  #   randomBeta(X, scale, m)
  #  
  #	X .......a matrix 
  #	scale ... a positive number
  #	m ....... a positive integer
  #
  #	For a matrix X of dimension n by p, this returns an p by m matrix 
  #	whose column-vectors are i.i.d. Each column vector v is such that
  #	X v is a vector of independent standard Gaussians within the column
  #	space of X, multiplied by scale.
  #	
  decomp <- svd(X);
  decomp$v %*% diag(1/decomp$d) %*% matrix(rnorm(m*dim(X)[2]), ncol=m) *scale
}
  
  
toCanonicalCoordinate = function(X) {
  #return a canonical coordinate version of X
  Q = svd(X)$u
  t(Q)%*%X
}


Knaive = function(r,alpha) {
  #constant for the naive student confidence intervals
  qt(p=1-alpha/2,df=r)
}
  

minCoverageProbability = function(nbeta1,nbeta2,nbatch3,nmc1,nmc2,nmc3,x0,X,Xc,sigma,Sigma,K,target,modelSelection,argModelSelection) {
  #calculate the min coverage probability with the Monte Carlo 3 step approach of Leeb, Poetscher, Ewald 2013.
  #nbeta1,nbeta2,nbatch3,nmc1,nmc2,nmc3: 2 on nbeta MC computation of size nmc. Step 3 on the best beta with nbatch*nmc MC evaluation, with nbatch time messages 
  #K: constant, or table of constants for the CI
  #modelSelection: a model selection procedure. See e.g. modelSelectionAIC
  #argModelSelection: see modelSelectionAIC, modelSelection LASSO
  
  
  #matrices of beta and vector of corresponding coverage probabilities
  mbeta1 = t(mRandomBeta(X,scale=1,m=nbeta1))
  mbeta2 = t(mRandomBeta(X,scale=1,m=nbeta2))
  vCov1 = seq(length=nbeta1)
  vCov2 = seq(length=nbeta2)
  vCov3 = seq(length=nbatch3)
  vHLCI1 = seq(length=nbeta1)
  vHLCI2 = seq(length=nbeta2)
  vHLCI3 = seq(length=nbatch3)
  vModSize1 = seq(length=nbeta1)
  vModSize2 = seq(length=nbeta2)
  vModSize3 = seq(length=nbatch3)
  #
  #step 1
  cat("########### \n")
  cat("Step 1 \n")
  cat("########### \n")
  for (ibeta1 in 1:nbeta1) {
    cat("beta number ",ibeta1,"\n")
    cat(Sys.time(),"\n")
    beta1 = mbeta1[ibeta1,1:p]
    lres1 = coverageProbability(x0,X,Sigma,beta1,sigma,nmc1,K,modelSelection,argModelSelection)
    if (target == 1) {
      vCov1[ibeta1] =lres1$coverageBetaMn
    }
    if (target == 2) {
      vCov1[ibeta1] =lres1$coverageBetaM
    }
    vHLCI1[ibeta1] =lres1$HLCI
    vModSize1[ibeta1] = lres1$modSize
  }
  
  #pick the best betas for step 2
  ind1 = sort(vCov1,decreasing=FALSE,index.return=TRUE)$ix
  mbeta2 = mbeta1[ind1[1:nbeta2],]
  
  
  #step 2
  cat("########### \n")
  cat("Step 2 \n")
  cat("########### \n")
  for (ibeta2 in 1:nbeta2) {
    cat("beta number ",ibeta2,"\n")
    cat(Sys.time(),"\n")
    beta2 = mbeta2[ibeta2,1:p]
    lres2 = coverageProbability(x0,X,Sigma,beta2,sigma,nmc2,K,modelSelection,argModelSelection)
    if (target == 1) {
      vCov2[ibeta2] =lres2$coverageBetaMn
    }
    if (target == 2) {
      vCov2[ibeta2] =lres2$coverageBetaM
    }
    vHLCI2[ibeta2] =lres2$HLCI
    vModSize2[ibeta2] = lres2$modSize
  }
  
  #pick the best beta for step 3
  ind2 = sort(vCov2,decreasing=FALSE,index.return=TRUE)$ix
  beta3 = mbeta2[ind2[1],]
  
  #step 3
  cat("########### \n")
  cat("Step 3 \n")
  cat("########### \n")
  for (ibeta3 in 1:nbatch3) {
    cat("Monte Carlo batch number ",ibeta3,"\n")
    cat(Sys.time(),"\n")
    lres3 = coverageProbability(x0,X,Sigma,beta3,sigma,nmc3,K,modelSelection,argModelSelection)
    if (target == 1) {
      vCov3[ibeta3] =lres3$coverageBetaMn
    }
    if (target == 2) {
      vCov3[ibeta3] =lres3$coverageBetaM
    }
    vHLCI3[ibeta3] = lres3$HLCI
    vModSize3[ibeta3] = lres3$modSize
  }
  
  #the min coverage probability
  minCovP = mean(vCov3)
  
  #return of the results
  list(mbeta1=mbeta1,mbeta2=mbeta2,beta3=beta3,vCov1=vCov1,vCov2=vCov2,vCov3=vCov3,vHLCI1=vHLCI1,vHLCI2=vHLCI2,vHLCI3=vHLCI3,vModSize1=vModSize1,vModSize2=vModSize2,vModSize3=vModSize3,minCovP=minCovP,target=target)
}
  


diffTargets = function(x0,X,Sigma,beta,sigma,I,modelSelection,argModelSelection) {
  #x0.................: the vector of new regressors
  #X..................: the design matrix (NOT in canonical coordinate)
  #Sigma..............: covariance matrix of the predictors
  #beta...............: true vector of coefficients
  #sigma..............: true standard deviation of Gaussian errors
  #I..................: number of Monte Carlo samples
  #modelSelection.....: the model selection function (AIC, BIC, LASSO)
  #argModelSelection..: For LASSO a conventional 0, for AIC, BIC, the factor k
  #return.............: vector of sample (size I) for the two targets
  #vDesignDep.........: design-dependent coverage target
  #vDesignInd.........: design-independent coverage target
  #
  #
  #
  #preparative
  n = dim(X)[1]
  p = dim(X)[2]
  Xc = toCanonicalCoordinate(X)
  mU = matrix(nrow=n,ncol=I,data=rnorm(n=n*I,mean=0,sd=sigma))  #sample of the I noise vectors
  vDesignDep = seq(length=I,from=1,to=1)  
  vDesignInd = seq(length=I,from=1,to=1)  
  #
  #Monte Carlo procedure
  for (imc in 1:I) {  #loop on the Monte Carlo sample
    Y = as.matrix(X)%*%matrix(nrow=p,ncol=1,data=beta) + mU[,imc]   #realization of the observation vector
    modelSelectionOutcome = modelSelection(X,Y,argModelSelection) # the selected model and the estimated beta
    model = modelSelectionOutcome$model
    hatbeta = modelSelectionOutcome$hatbeta 
    if (length(hatbeta)>0) {   #if the selected model is not empty, then there are no targets
      vDesignDep[imc] = sum(x0[model]*betaMn(Xc,beta,model))  
      vDesignInd[imc] = sum(x0[model]*betaM(Sigma,beta,model))  
    }
  }
  list(vDesignDep=vDesignDep,vDesignInd=vDesignInd)
}


allTarForMod = function(x0,X,Sigma,beta_) {
  #Give all the target values for all possible submodels containing the 1st regressor
  # return: (2^(p-1) -1)X(p+3) matrix. Line: a submodel. Column: 1-p: the submodel. p+1:design-dependant target. p+2: design-independent target. p+3: norm tM
  #
  n = dim(X)[1]
  p = dim(X)[2]
  #
  mTar = matrix(nrow = (2^(p-1)),ncol=(p+3),data=0)
  #
  #submodel with only intercept
  model = seq(from=0,to=0,length=p)
  model[1] = 1
  mTar[1,1:p] = model
  mTar[1,p+1] = sum(x0[model]*betaMn(Xc,beta_,model))  
  mTar[1,p+2] = sum(x0[model]*betaM(Sigma,beta_,model))  
  mTar[1,p+3] = sqrt( sum(tM(x0,X,model[model!=0])^2) )
  #
  index=2
  for (m.sz in 1:(p-1)) {   #loop on size of models (EXCLUDING the intercept)
    if (m.sz<p-1){ #obtaining the matrix of all the models of a given size
      models.of.sz <- subsets(p-1,m.sz)
    }else{
      models.of.sz <- matrix(1:(p-1),1,(p-1))
    }
    for (i in 1:dim(models.of.sz)[1]) {   #loop on models of a given size
      model = models.of.sz[i,]
      model = c(1,model+1) #adding the 1st regressor
      model = c(model,seq(from=0,to=0,length=p-length(model)))
      mTar[index,1:p] = model
      mTar[index,p+1] = sum(x0[model]*betaMn(Xc,beta_,model))  
      mTar[index,p+2] = sum(x0[model]*betaM(Sigma,beta_,model))  
      mTar[index,p+3] = sqrt( sum(tM(x0,X,model[model!=0])^2) )
      index = index+1
    }
  }
  mTar
}

mMod = function(x0,X,Sigma,beta,sigma,I,modelSelection,argModelSelection) {
  #x0.................: the vector of new regressors
  #X..................: the design matrix (NOT in canonical coordinate)
  #Sigma..............: covariance matrix of the predictors
  #beta...............: true vector of coefficients
  #sigma..............: true standard deviation of Gaussian errors
  #I..................: number of Monte Carlo samples
  #modelSelection.....: the model selection function (AIC, BIC, LASSO)
  #argModelSelection..: For LASSO a conventional 0, for AIC, BIC, the factor k
  #return.............: matrix of sample (size IXp) for the selected model  #
  #
  #
  #preparative
  n = dim(X)[1]
  p = dim(X)[2]
  Xc = toCanonicalCoordinate(X)
  mU = matrix(nrow=n,ncol=I,data=rnorm(n=n*I,mean=0,sd=sigma))  #sample of the I noise vectors
  mMod_ = matrix(nrow=I,ncol=p,data=0)
  #
  #Monte Carlo procedure
  for (imc in 1:I) {  #loop on the Monte Carlo sample
    Y = as.matrix(X)%*%matrix(nrow=p,ncol=1,data=beta) + mU[,imc]   #realization of the observation vector
    modelSelectionOutcome = modelSelection(X,Y,argModelSelection) # the selected model and the estimated beta
    model = modelSelectionOutcome$model
    mMod_[imc,] = model
  }
  mMod_
}



  