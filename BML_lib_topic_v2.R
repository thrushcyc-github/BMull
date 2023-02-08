BML_binary_same_alpha = function(XX,yy,z,XtX,niter = 5000,p=5,id_train,nu=5*3+2,psi=10,cutoff = 4){
  # BML 1: suppose alpha is the same (model 1)
  ## now start the sampling:
   
  Beta = matrix(0,niter,p*3)
  Alpha = matrix(0,niter,4)
  Ys = matrix(0,niter,length(id_train))
  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
     
  Psi = diag(rep(psi,p*3))
  Sigma_list = {}
     
  m <- polr(as.factor(y)~X, Hess=TRUE, method="probit")
  beta = m$coefficients
  alpha = m$zeta
  for(i in 1:niter){
    if(i%%5000==0) print(i)
    ## samlpe y_star
    low = c(-Inf,alpha[1:4])
    up = c(alpha[1:4],Inf)
    y_star = rtruncnorm(length(y), a=low[y], b=up[y], mean = X%*%beta, sd = 1)

    ## sample Sigma ## should this Sigma be block diagonal?
    #Sigma = riwish(nu+.5, Psi+beta%*%t(beta))
    A = Psi+beta%*%t(beta)
    Sigma = matrix(0,p*3,p*3) 
    for(k in 1:p) Sigma[1:3+3*(k-1),1:3+3*(k-1)] = riwish(nu+1, A[1:3+3*(k-1),1:3+3*(k-1)])
    
    # sample alpha
    alpha1 = rtruncnorm(1, a=max(y_star[y==1]), b=min(y_star[y==2]), mean = 0, sd = sqrt(10))
    alpha2 = rtruncnorm(1, a=max(y_star[y==2]), b=min(y_star[y==3]), mean = 0, sd = sqrt(10))
    alpha3 = rtruncnorm(1, a=max(y_star[y==3]), b=min(y_star[y==4]), mean = 0, sd = sqrt(10))
    alpha4 = rtruncnorm(1, a=max(y_star[y==4]), b=min(y_star[y==5]), mean = 0, sd = sqrt(10))
    alpha = c(alpha1,alpha2,alpha3,alpha4)
    ## sample beta
    A = solve(XtX+solve(Sigma))
    beta = mvrnorm(n = 1, mu = A%*%(t(X)%*%y_star), Sigma = A)
    ## store the results
    Beta[i,] = beta
    Alpha[i,]=alpha
    Ys[i,]=y_star
    Sigma_list[[i]] = Sigma
  }
  ## compare the performance
  alpha_post = apply(Alpha[-(1:1000),],2,mean)
  beta_post = apply(Beta[-(1:1000),],2,mean)
  ## for model 1
  mu = Xt%*%beta_post
  p11 = pnorm(alpha_post[1]-mu)
  p12 = pnorm(alpha_post[2]-mu) - p11
  p13 = pnorm(alpha_post[3]-mu) - pnorm(alpha_post[2]-mu)
  p14 = pnorm(alpha_post[4]-mu) - pnorm(alpha_post[3]-mu)
  p15 = 1 - pnorm(alpha_post[4]-mu)
  p1 = cbind(p11,p12,p13,p14,p15)
  #y1 = apply(p1,1,which.max)
  out = list()
  y1 = apply(p1[,3:5],1,sum)
  y2 = apply(p1[,4:5],1,sum)
  out$y1 = y1
  out$y2 = y2
  return(out)
}

BML_binary_diff_alpha = function(XX,yy,zz,XtX,niter = 5000,p=5,id_train,nu=5*3+2,psi=10,cutoff=4){
  Beta = matrix(0,niter,p*3)
  Alpha = matrix(0,niter,4*3)
  Ys = matrix(0,niter,length(id_train))
  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
  z = zz[id_train]
  zt = zz[-id_train]
  
  Psi = diag(rep(psi,p*3))
  Sigma_list = {}
  
  m <- polr(as.factor(y)~X, Hess=TRUE, method="probit")
  beta = m$coefficients
  alpha = rep(m$zeta,3)
  for(i in 1:niter){
    if(i%%5000==0) print(i)
    ## samlpe y_star
    low = c(-Inf,alpha[1:4])
    up = c(alpha[1:4],Inf)
    y_star = rtruncnorm(length(y), a=low[y], b=up[y], mean = X%*%beta, sd = 1)
    
    low = c(-Inf,alpha[1:4+4])
    up = c(alpha[1:4+4],Inf)
    y_star[z==2] = rtruncnorm(length(y[z==2]), a=low[y[z==2]], b=up[y[z==2]], mean = X[which(z==2),]%*%beta, sd = 1)
    
    low = c(-Inf,alpha[1:4+8])
    up = c(alpha[1:4+8],Inf)
    y_star[z==3] = rtruncnorm(length(y[z==3]), a=low[y[z==3]], b=up[y[z==3]], mean = X[which(z==3),]%*%beta, sd = 1)
    
    ## sample Sigma ## should this Sigma be block diagonal?
    #Sigma = riwish(nu+.5, Psi+beta%*%t(beta))
    A = Psi+beta%*%t(beta)
    Sigma = matrix(0,p*3,p*3) 
    for(k in 1:p) Sigma[1:3+3*(k-1),1:3+3*(k-1)] = riwish(nu+1, A[1:3+3*(k-1),1:3+3*(k-1)])
    
    # sample alpha
    for(j in 1:3){
      alpha[1+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==1&z==j]), b=min(y_star[y==2&z==j]), mean = 0, sd = sqrt(10))
      alpha[2+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==2&z==j]), b=min(y_star[y==3&z==j]), mean = 0, sd = sqrt(10))
      alpha[3+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==3&z==j]), b=min(y_star[y==4&z==j]), mean = 0, sd = sqrt(10))
      alpha[4+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==4&z==j]), b=min(y_star[y==5&z==j]), mean = 0, sd = sqrt(10))
    }

    ## sample beta
    A = solve(XtX+solve(Sigma))
    beta = mvrnorm(n = 1, mu = A%*%(t(X)%*%y_star), Sigma = A)
    
    ## store the results
    Beta[i,] = beta
    Alpha[i,]=alpha
    Ys[i,]=y_star
    Sigma_list[[i]] = Sigma
  }
  ## compare the performance
  alpha_post = apply(Alpha[-(1:1000),],2,mean)
  beta_post = apply(Beta[-(1:1000),],2,mean)
  
  ## for model 1
  mu = Xt%*%beta_post
  p11 = pnorm(alpha_post[1+(zt-1)*4]-mu)
  p12 = pnorm(alpha_post[2+(zt-1)*4]-mu) - p11
  p13 = pnorm(alpha_post[3+(zt-1)*4]-mu) - pnorm(alpha_post[2+(zt-1)*4]-mu)
  p14 = pnorm(alpha_post[4+(zt-1)*4]-mu) - pnorm(alpha_post[3+(zt-1)*4]-mu)
  p15 = 1 - pnorm(alpha_post[4+(zt-1)*4]-mu)
  
  p1 = cbind(p11,p12,p13,p14,p15)
  #y1 = apply(p1,1,which.max)
  out = list()
  y1 = apply(p1[,3:5],1,sum)
  y2 = apply(p1[,4:5],1,sum)
  out$y1 = y1
  out$y2 = y2
  return(out)
}

BML_binary_same_alpha_ind_sig = function(XX,yy,z,XtX,niter = 5000,p=5,id_train,nu=5*3+2,psi=10,cutoff = 4){
  # BML 1: suppose alpha is the same (model 1)
  ## now start the sampling:
  
  Beta = matrix(0,niter,p*3)
  Alpha = matrix(0,niter,4)
  Ys = matrix(0,niter,length(id_train))
  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
  
  Psi = diag(rep(psi,p*3))
  Sigma_list = {}
  
  m <- polr(as.factor(y)~X, Hess=TRUE, method="probit")
  beta = m$coefficients
  alpha = m$zeta
  for(i in 1:niter){
    if(i%%5000==0) print(i)
    ## samlpe y_star
    low = c(-Inf,alpha[1:4])
    up = c(alpha[1:4],Inf)
    y_star = rtruncnorm(length(y), a=low[y], b=up[y], mean = X%*%beta, sd = 1)
    
    ## sample Sigma ## should this Sigma be block diagonal?
    #Sigma = riwish(nu+.5, Psi+beta%*%t(beta))
    A = Psi+beta%*%t(beta)
    Sigma = matrix(0,p*3,p*3) 
    diag(Sigma)=1
    #for(k in 1:p) Sigma[1:3+3*(k-1),1:3+3*(k-1)] = riwish(nu+1, A[1:3+3*(k-1),1:3+3*(k-1)])
    
    # sample alpha
    alpha1 = rtruncnorm(1, a=max(y_star[y==1]), b=min(y_star[y==2]), mean = 0, sd = sqrt(10))
    alpha2 = rtruncnorm(1, a=max(y_star[y==2]), b=min(y_star[y==3]), mean = 0, sd = sqrt(10))
    alpha3 = rtruncnorm(1, a=max(y_star[y==3]), b=min(y_star[y==4]), mean = 0, sd = sqrt(10))
    alpha4 = rtruncnorm(1, a=max(y_star[y==4]), b=min(y_star[y==5]), mean = 0, sd = sqrt(10))
    alpha = c(alpha1,alpha2,alpha3,alpha4)
    ## sample beta
    A = solve(XtX+solve(Sigma))
    beta = mvrnorm(n = 1, mu = A%*%(t(X)%*%y_star), Sigma = A)
    ## store the results
    Beta[i,] = beta
    Alpha[i,]=alpha
    Ys[i,]=y_star
    Sigma_list[[i]] = Sigma
  }
  ## compare the performance
  alpha_post = apply(Alpha[-(1:1000),],2,mean)
  beta_post = apply(Beta[-(1:1000),],2,mean)
  ## for model 1
  mu = Xt%*%beta_post
  p11 = pnorm(alpha_post[1]-mu)
  p12 = pnorm(alpha_post[2]-mu) - p11
  p13 = pnorm(alpha_post[3]-mu) - pnorm(alpha_post[2]-mu)
  p14 = pnorm(alpha_post[4]-mu) - pnorm(alpha_post[3]-mu)
  p15 = 1 - pnorm(alpha_post[4]-mu)
  p1 = cbind(p11,p12,p13,p14,p15)
  #y1 = apply(p1,1,which.max)
  out = list()
  y1 = apply(p1[,3:5],1,sum)
  y2 = apply(p1[,4:5],1,sum)
  out$y1 = y1
  out$y2 = y2
  return(out)
}

BML_binary_diff_alpha_ind_sig = function(XX,yy,zz,XtX,niter = 5000,p=5,id_train,nu=5*3+2,psi=10,cutoff=4){
  Beta = matrix(0,niter,p*3)
  Alpha = matrix(0,niter,4*3)
  Ys = matrix(0,niter,length(id_train))
  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
  z = zz[id_train]
  zt = zz[-id_train]
  
  Psi = diag(rep(psi,p*3))
  Sigma_list = {}
  
  m <- polr(as.factor(y)~X, Hess=TRUE, method="probit")
  beta = m$coefficients
  alpha = rep(m$zeta,3)
  for(i in 1:niter){
    if(i%%5000==0) print(i)
    ## samlpe y_star
    low = c(-Inf,alpha[1:4])
    up = c(alpha[1:4],Inf)
    y_star = rtruncnorm(length(y), a=low[y], b=up[y], mean = X%*%beta, sd = 1)
    
    low = c(-Inf,alpha[1:4+4])
    up = c(alpha[1:4+4],Inf)
    y_star[z==2] = rtruncnorm(length(y[z==2]), a=low[y[z==2]], b=up[y[z==2]], mean = X[which(z==2),]%*%beta, sd = 1)
    
    low = c(-Inf,alpha[1:4+8])
    up = c(alpha[1:4+8],Inf)
    y_star[z==3] = rtruncnorm(length(y[z==3]), a=low[y[z==3]], b=up[y[z==3]], mean = X[which(z==3),]%*%beta, sd = 1)
    
    ## sample Sigma ## should this Sigma be block diagonal?
    #Sigma = riwish(nu+.5, Psi+beta%*%t(beta))
    A = Psi+beta%*%t(beta)
    Sigma = matrix(0,p*3,p*3)
    diag(Sigma)=1
    #for(k in 1:p) Sigma[1:3+3*(k-1),1:3+3*(k-1)] = riwish(nu+1, A[1:3+3*(k-1),1:3+3*(k-1)])
    
    # sample alpha
    for(j in 1:3){
      alpha[1+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==1&z==j]), b=min(y_star[y==2&z==j]), mean = 0, sd = sqrt(10))
      alpha[2+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==2&z==j]), b=min(y_star[y==3&z==j]), mean = 0, sd = sqrt(10))
      alpha[3+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==3&z==j]), b=min(y_star[y==4&z==j]), mean = 0, sd = sqrt(10))
      alpha[4+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==4&z==j]), b=min(y_star[y==5&z==j]), mean = 0, sd = sqrt(10))
    }
    
    ## sample beta
    A = solve(XtX+solve(Sigma))
    beta = mvrnorm(n = 1, mu = A%*%(t(X)%*%y_star), Sigma = A)
    
    ## store the results
    Beta[i,] = beta
    Alpha[i,]=alpha
    Ys[i,]=y_star
    Sigma_list[[i]] = Sigma
  }
  ## compare the performance
  alpha_post = apply(Alpha[-(1:1000),],2,mean)
  beta_post = apply(Beta[-(1:1000),],2,mean)
  
  ## for model 1
  mu = Xt%*%beta_post
  p11 = pnorm(alpha_post[1+(zt-1)*4]-mu)
  p12 = pnorm(alpha_post[2+(zt-1)*4]-mu) - p11
  p13 = pnorm(alpha_post[3+(zt-1)*4]-mu) - pnorm(alpha_post[2+(zt-1)*4]-mu)
  p14 = pnorm(alpha_post[4+(zt-1)*4]-mu) - pnorm(alpha_post[3+(zt-1)*4]-mu)
  p15 = 1 - pnorm(alpha_post[4+(zt-1)*4]-mu)
  
  p1 = cbind(p11,p12,p13,p14,p15)
  #y1 = apply(p1,1,which.max)
  out = list()
  y1 = apply(p1[,3:5],1,sum)
  y2 = apply(p1[,4:5],1,sum)
  out$y1 = y1
  out$y2 = y2
  return(out)
}


BML_binary_diff_alpha_v2 = function(XX,yy,zz,XtX,niter = 5000,p=5,id_train,nu=5*3+2,psi=10,cutoff=4,data,var_str=var_str){
  ## different from v1, the initial points are obtained from three regression
  nn = nrow(data)
  id_train1 = id_train[which(id_train<=nn)]
  id_train2 = id_train[which(id_train<=nn*2&id_train>nn)]-nn
  id_train3 = id_train[which(id_train>nn*2)]-nn*2

  m1 <- polr(paste("satisfaction.fac~",var_str,sep=""),data=data[id_train1,], Hess=TRUE, method="probit")
  m2 <- polr(paste("ease.fac~",var_str,sep=""),data=data[id_train2,], Hess=TRUE, method="probit")
  m3 <- polr(paste("effect.fac~",var_str,sep=""),data=data[id_train3,], Hess=TRUE, method="probit")
  
  beta = c(m1$coefficients,m2$coefficients,m3$coefficients)
  alpha = c(m1$zeta,m2$zeta,m3$zeta)
  
  ## different from v1, the initial points are obtained from three regression
  Beta = matrix(0,niter,p*3)
  Alpha = matrix(0,niter,4*3)
  Ys = matrix(0,niter,length(id_train))
  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
  z = zz[id_train]
  zt = zz[-id_train]
  
  Psi = diag(rep(psi,p*3))
  Sigma_list = {}
  
  for(i in 1:niter){
    if(i%%5000==0) print(i)
    ## samlpe y_star
    low = c(-Inf,alpha[1:4])
    up = c(alpha[1:4],Inf)
    y_star = rtruncnorm(length(y), a=low[y], b=up[y], mean = X%*%beta, sd = 1)
    
    low = c(-Inf,alpha[1:4+4])
    up = c(alpha[1:4+4],Inf)
    y_star[z==2] = rtruncnorm(length(y[z==2]), a=low[y[z==2]], b=up[y[z==2]], mean = X[which(z==2),]%*%beta, sd = 1)
    
    low = c(-Inf,alpha[1:4+8])
    up = c(alpha[1:4+8],Inf)
    y_star[z==3] = rtruncnorm(length(y[z==3]), a=low[y[z==3]], b=up[y[z==3]], mean = X[which(z==3),]%*%beta, sd = 1)
    
    ## sample Sigma ## should this Sigma be block diagonal?
    #Sigma = riwish(nu+.5, Psi+beta%*%t(beta))
    A = Psi+beta%*%t(beta)
    Sigma = matrix(0,p*3,p*3) 
    for(k in 1:p) Sigma[1:3+3*(k-1),1:3+3*(k-1)] = riwish(nu+1, A[1:3+3*(k-1),1:3+3*(k-1)])
    
    # sample alpha
    for(j in 1:3){
      alpha[1+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==1&z==j]), b=min(y_star[y==2&z==j]), mean = 0, sd = sqrt(10))
      alpha[2+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==2&z==j]), b=min(y_star[y==3&z==j]), mean = 0, sd = sqrt(10))
      alpha[3+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==3&z==j]), b=min(y_star[y==4&z==j]), mean = 0, sd = sqrt(10))
      alpha[4+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==4&z==j]), b=min(y_star[y==5&z==j]), mean = 0, sd = sqrt(10))
    }
    
    ## sample beta
    A = solve(XtX+solve(Sigma))
    beta = mvrnorm(n = 1, mu = A%*%(t(X)%*%y_star), Sigma = A)
    
    ## store the results
    Beta[i,] = beta
    Alpha[i,]=alpha
    Ys[i,]=y_star
    Sigma_list[[i]] = Sigma
  }
  ## compare the performance
  alpha_post = apply(Alpha[-(1:1000),],2,mean)
  beta_post = apply(Beta[-(1:1000),],2,mean)
  
  ## for model 1
  mu = Xt%*%beta_post
  p11 = pnorm(alpha_post[1+(zt-1)*4]-mu)
  p12 = pnorm(alpha_post[2+(zt-1)*4]-mu) - p11
  p13 = pnorm(alpha_post[3+(zt-1)*4]-mu) - pnorm(alpha_post[2+(zt-1)*4]-mu)
  p14 = pnorm(alpha_post[4+(zt-1)*4]-mu) - pnorm(alpha_post[3+(zt-1)*4]-mu)
  p15 = 1 - pnorm(alpha_post[4+(zt-1)*4]-mu)
  
  p1 = cbind(p11,p12,p13,p14,p15)
  #y1 = apply(p1,1,which.max)
  out = list()
  y1 = apply(p1[,3:5],1,sum)
  y2 = apply(p1[,4:5],1,sum)
  out$y1 = y1
  out$y2 = y2
  return(out)
  
}


BML_binary_same_alpha_collapse_y = function(XX,yy,z,XtX,niter = 5000,p=5,id_train,nu=5*3+2,psi=10,cutoff=4){
  # BML 1: suppose alpha is the same (model 1)
  ## now start the sampling:
  yy[yy<cutoff] = 1
  yy[yy>=cutoff] = 2
  
  Beta = matrix(0,niter,p*3)
  Alpha = matrix(0,niter,1)

  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
  
  Psi = diag(rep(psi,p*3))
  Sigma_list = {}
  
  #m <- polr(as.factor(y)~X, Hess=TRUE, method="probit")
  m = glm((y-1)~X,family = binomial(link="probit"))
  beta = m$coefficients[-1]
  alpha = m$coefficients[1]
  for(i in 1:niter){
    if(i%%5000==0) print(i)
    ## samlpe y_star
    low = c(-Inf,alpha)
    up = c(alpha,Inf)
    y_star = rtruncnorm(length(y), a=low[y], b=up[y], mean = X%*%beta, sd = 1)

    ## sample Sigma ## should this Sigma be block diagonal?
    #Sigma = riwish(nu+.5, Psi+beta%*%t(beta))
    A = Psi+beta%*%t(beta)
    Sigma = matrix(0,p*3,p*3) 
    for(k in 1:p) Sigma[1:3+3*(k-1),1:3+3*(k-1)] = riwish(nu+1, A[1:3+3*(k-1),1:3+3*(k-1)])
    
    # sample alpha
    alpha = rtruncnorm(1, a=max(y_star[y==1]), b=min(y_star[y==2]), mean = 0, sd = sqrt(10))

    ## sample beta
    A = solve(XtX+solve(Sigma))
    beta = mvrnorm(n = 1, mu = A%*%(t(X)%*%y_star), Sigma = A)
    ## store the results
    Beta[i,] = beta
    Alpha[i,]=alpha
    #Ys[i,]=y_star
    Sigma_list[[i]] = Sigma
  }
  ## compare the performance
  #alpha_post = apply(Alpha[-(1:1000),],2,mean)
  alpha_post = mean(Alpha[-(1:1000),])
  beta_post = apply(Beta[-(1:1000),],2,mean)
  ## for model 1
  mu = Xt%*%beta_post
  p1 = 1-pnorm(alpha_post[1]-mu)
  return(p1)
}

BML_binary_diff_alpha_collapse_y = function(XX,yy,zz,XtX,niter = 5000,p=5,id_train,nu=5*3+2,psi=10,cutoff=4){
  Beta = matrix(0,niter,p*3)
  Alpha = matrix(0,niter,3)
  
  yy[yy<cutoff] = 1
  yy[yy>=cutoff] = 2
  
  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
  z = zz[id_train]
  zt = zz[-id_train]
  
  Psi = diag(rep(psi,p*3))
  Sigma_list = {}
  
  #m <- polr(as.factor(y)~X, Hess=TRUE, method="probit")
  m = glm((y-1)~X,family = binomial(link="probit"))
  
  beta = m$coefficients[-1]
  alpha = rep(m$coefficients[1],3)
  for(i in 1:niter){
    if(i%%5000==0) print(i)
    ## samlpe y_star
    low = c(-Inf,alpha[1])
    up = c(alpha[1],Inf)
    y_star = rtruncnorm(length(y), a=low[y], b=up[y], mean = X%*%beta, sd = 1)
    
    low = c(-Inf,alpha[2])
    up = c(alpha[2],Inf)
    y_star[z==2] = rtruncnorm(length(y[z==2]), a=low[y[z==2]], b=up[y[z==2]], mean = X[which(z==2),]%*%beta, sd = 1)
    
    low = c(-Inf,alpha[3])
    up = c(alpha[3],Inf)
    y_star[z==3] = rtruncnorm(length(y[z==3]), a=low[y[z==3]], b=up[y[z==3]], mean = X[which(z==3),]%*%beta, sd = 1)
    
    ## sample Sigma ## should this Sigma be block diagonal?
    #Sigma = riwish(nu+.5, Psi+beta%*%t(beta))
    A = Psi+beta%*%t(beta)
    Sigma = matrix(0,p*3,p*3) 
    for(k in 1:p) Sigma[1:3+3*(k-1),1:3+3*(k-1)] = riwish(nu+1, A[1:3+3*(k-1),1:3+3*(k-1)])
    
    # sample alpha
    for(j in 1:3){
      alpha[j] = rtruncnorm(1, a=max(y_star[y==1&z==j]), b=min(y_star[y==2&z==j]), mean = 0, sd = sqrt(10))
    }
    
    
    ## sample beta
    A = solve(XtX+solve(Sigma))
    beta = mvrnorm(n = 1, mu = A%*%(t(X)%*%y_star), Sigma = A)
    
    ## store the results
    Beta[i,] = beta
    Alpha[i,]=alpha
    #Ys[i,]=y_star
    Sigma_list[[i]] = Sigma
  }
  ## compare the performance
  alpha_post = apply(Alpha[-(1:1000),],2,mean)
  beta_post = apply(Beta[-(1:1000),],2,mean)
  
  ## for model 1
  mu = Xt%*%beta_post
  p1 = 1 - pnorm(alpha_post[zt]-mu)

  return(p1)
}

BML_binary_diff_alpha_v2_collapse_y = function(XX,yy,zz,XtX,niter = 5000,p=5,id_train,nu=5*3+2,psi=10,cutoff=4){
  yy[yy<cutoff] = 1
  yy[yy>=cutoff] = 2
  
  nn = nrow(data)
  id_train1 = id_train[which(id_train<=nn)]
  id_train2 = id_train[which(id_train<=nn*2&id_train>nn)]
  id_train3 = id_train[which(id_train>nn*2)]
  
  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
  z = zz[id_train]
  zt = zz[-id_train]
  
  m1 = glm((yy[id_train1]-1)~XX[id_train1,1:p],family = binomial(link="probit"))
  m2 = glm((yy[id_train2]-1)~XX[id_train2,1:p+p],family = binomial(link="probit"))
  m3 = glm((yy[id_train3]-1)~XX[id_train3,1:p+2*p],family = binomial(link="probit"))
  
  beta = c(m1$coefficients[-1],m2$coefficients[-1],m3$coefficients[-1])
  alpha = c(m1$coefficients[1],m2$coefficients[1],m3$coefficients[1])
  
  ## different from v1, the initial points are obtained from three regression
  Beta = matrix(0,niter,p*3)
  Alpha = matrix(0,niter,3)
  Ys = matrix(0,niter,length(id_train))
  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
  z = zz[id_train]
  zt = zz[-id_train]
  
  Psi = diag(rep(psi,p*3))
  Sigma_list = {}
  
  for(i in 1:niter){
    if(i%%5000==0) print(i)
    ## samlpe y_star
    low = c(-Inf,alpha[1])
    up = c(alpha[1],Inf)
    y_star = rtruncnorm(length(y), a=low[y], b=up[y], mean = X%*%beta, sd = 1)
    
    low = c(-Inf,alpha[2])
    up = c(alpha[2],Inf)
    y_star[z==2] = rtruncnorm(length(y[z==2]), a=low[y[z==2]], b=up[y[z==2]], mean = X[which(z==2),]%*%beta, sd = 1)
    
    low = c(-Inf,alpha[3])
    up = c(alpha[3],Inf)
    y_star[z==3] = rtruncnorm(length(y[z==3]), a=low[y[z==3]], b=up[y[z==3]], mean = X[which(z==3),]%*%beta, sd = 1)
    
    ## sample Sigma ## should this Sigma be block diagonal?
    #Sigma = riwish(nu+.5, Psi+beta%*%t(beta))
    A = Psi+beta%*%t(beta)
    Sigma = matrix(0,p*3,p*3) 
    for(k in 1:p) Sigma[1:3+3*(k-1),1:3+3*(k-1)] = riwish(nu+1, A[1:3+3*(k-1),1:3+3*(k-1)])
    
    # sample alpha
    for(j in 1:3){
      alpha[j] = rtruncnorm(1, a=max(y_star[y==1&z==j]), b=min(y_star[y==2&z==j]), mean = 0, sd = sqrt(10))
    }
    
    ## sample beta
    A = solve(XtX+solve(Sigma))
    beta = mvrnorm(n = 1, mu = A%*%(t(X)%*%y_star), Sigma = A)
    
    ## store the results
    Beta[i,] = beta
    Alpha[i,]=alpha
    Ys[i,]=y_star
    Sigma_list[[i]] = Sigma
  }
  ## compare the performance
  alpha_post = apply(Alpha[-(1:1000),],2,mean)
  beta_post = apply(Beta[-(1:1000),],2,mean)
  
  ## for model 1
  mu = Xt%*%beta_post
  p1 = 1 - pnorm(alpha_post[zt]-mu)
  
  return(p1)
}

BML_binary_same_alpha_ind_sig_collapse_y = function(XX,yy,z,XtX,niter = 5000,p=5,id_train,nu=5*3+2,psi=10,cutoff=4){
  # BML 1: suppose alpha is the same (model 1)
  ## now start the sampling:
  yy[yy<cutoff] = 1
  yy[yy>=cutoff] = 2
  
  Beta = matrix(0,niter,p*3)
  Alpha = matrix(0,niter,1)
  
  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
  
  Psi = diag(rep(psi,p*3))
  Sigma_list = {}
  
  #m <- polr(as.factor(y)~X, Hess=TRUE, method="probit")
  m = glm((y-1)~X,family = binomial(link="probit"))
  beta = m$coefficients[-1]
  alpha = m$coefficients[1]
  for(i in 1:niter){
    if(i%%5000==0) print(i)
    ## samlpe y_star
    low = c(-Inf,alpha)
    up = c(alpha,Inf)
    y_star = rtruncnorm(length(y), a=low[y], b=up[y], mean = X%*%beta, sd = 1)
    
    ## sample Sigma ## should this Sigma be block diagonal?
    #Sigma = riwish(nu+.5, Psi+beta%*%t(beta))
    A = Psi+beta%*%t(beta)
    Sigma = matrix(0,p*3,p*3) 
    diag(Sigma) = 1
    #for(k in 1:p) Sigma[1:3+3*(k-1),1:3+3*(k-1)] = riwish(nu+1, A[1:3+3*(k-1),1:3+3*(k-1)])
    
    # sample alpha
    alpha = rtruncnorm(1, a=max(y_star[y==1]), b=min(y_star[y==2]), mean = 0, sd = sqrt(10))
    
    ## sample beta
    A = solve(XtX+solve(Sigma))
    beta = mvrnorm(n = 1, mu = A%*%(t(X)%*%y_star), Sigma = A)
    ## store the results
    Beta[i,] = beta
    Alpha[i,]=alpha
    #Ys[i,]=y_star
    Sigma_list[[i]] = Sigma
  }
  ## compare the performance
  #alpha_post = apply(Alpha[-(1:1000),],2,mean)
  alpha_post = mean(Alpha[-(1:1000),])
  beta_post = apply(Beta[-(1:1000),],2,mean)
  ## for model 1
  mu = Xt%*%beta_post
  p1 = 1-pnorm(alpha_post[1]-mu)
  return(p1)
}

BML_binary_diff_alpha_ind_sig_collapse_y = function(XX,yy,zz,XtX,niter = 5000,p=5,id_train,nu=5*3+2,psi=10,cutoff=4){
  Beta = matrix(0,niter,p*3)
  Alpha = matrix(0,niter,3)
  
  yy[yy<cutoff] = 1
  yy[yy>=cutoff] = 2
  
  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
  z = zz[id_train]
  zt = zz[-id_train]
  
  Psi = diag(rep(psi,p*3))
  Sigma_list = {}
  
  #m <- polr(as.factor(y)~X, Hess=TRUE, method="probit")
  m = glm((y-1)~X,family = binomial(link="probit"))
  
  beta = m$coefficients[-1]
  alpha = rep(m$coefficients[1],3)
  for(i in 1:niter){
    if(i%%5000==0) print(i)
    ## samlpe y_star
    low = c(-Inf,alpha[1])
    up = c(alpha[1],Inf)
    y_star = rtruncnorm(length(y), a=low[y], b=up[y], mean = X%*%beta, sd = 1)
    
    low = c(-Inf,alpha[2])
    up = c(alpha[2],Inf)
    y_star[z==2] = rtruncnorm(length(y[z==2]), a=low[y[z==2]], b=up[y[z==2]], mean = X[which(z==2),]%*%beta, sd = 1)
    
    low = c(-Inf,alpha[3])
    up = c(alpha[3],Inf)
    y_star[z==3] = rtruncnorm(length(y[z==3]), a=low[y[z==3]], b=up[y[z==3]], mean = X[which(z==3),]%*%beta, sd = 1)
    
    ## sample Sigma ## should this Sigma be block diagonal?
    #Sigma = riwish(nu+.5, Psi+beta%*%t(beta))
    A = Psi+beta%*%t(beta)
    Sigma = matrix(0,p*3,p*3)
    diag(Sigma) = 1
    #for(k in 1:p) Sigma[1:3+3*(k-1),1:3+3*(k-1)] = riwish(nu+1, A[1:3+3*(k-1),1:3+3*(k-1)])
    
    # sample alpha
    for(j in 1:3){
      alpha[j] = rtruncnorm(1, a=max(y_star[y==1&z==j]), b=min(y_star[y==2&z==j]), mean = 0, sd = sqrt(10))
    }
    
    
    ## sample beta
    A = solve(XtX+solve(Sigma))
    beta = mvrnorm(n = 1, mu = A%*%(t(X)%*%y_star), Sigma = A)
    
    ## store the results
    Beta[i,] = beta
    Alpha[i,]=alpha
    #Ys[i,]=y_star
    Sigma_list[[i]] = Sigma
  }
  ## compare the performance
  alpha_post = apply(Alpha[-(1:1000),],2,mean)
  beta_post = apply(Beta[-(1:1000),],2,mean)
  
  ## for model 1
  mu = Xt%*%beta_post
  p1 = 1 - pnorm(alpha_post[zt]-mu)
  
  return(p1)
}

BML_lasso_binary_same_alpha = function(XX,yy,z,XtX,niter = 5000,p=5,id_train,nu=5*3+2,psi=10,cutoff = 4){
  # BML 1: suppose alpha is the same (model 1)
  ## now start the sampling:
  
  Beta = matrix(0,niter,p*3)
  Alpha = matrix(0,niter,4)
  Tau = matrix(0,niter,p) ## here the tau is tau square
  gam = c() ## here gam is lambda (copied from BEL)
  
  Ys = matrix(0,niter,length(id_train))
  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
  
  Psi = diag(rep(psi,p*3))
  Sigma_list = {}
  
  m <- polr(as.factor(y)~X, Hess=TRUE, method="probit")
  beta = m$coefficients
  alpha = m$zeta
  
  tau = rep(1,p)
  a = diag(rep(tau,each=3))
  Tau[1,] = tau
  
  for(i in 1:niter){
    if(i%%5000==0) print(i)
    
    if(i>1){
      ## sample tau
      gam[i] = sqrt(2*p/sum(Tau[1:(i-1),])*(i-1))
      for(k in 1:p){
        temp = Beta[i-1,1:3+3*(k-1)]
        m_temp = (Sigma_list[[i-1]])[1:3+3*(k-1),1:3+3*(k-1)]
        d = t(temp)%*%m_temp%*%temp
        Tau[i,k] = 1/rinvgauss(n=1,m=pmax(sqrt(gam[i]^2/d),.1),s=gam[i]^2)
      }
      tau = Tau[i,]
      a = diag(rep(tau,each=3))
    }
    
    ## samlpe y_star
    low = c(-Inf,alpha[1:4])
    up = c(alpha[1:4],Inf)
    y_star = rtruncnorm(length(y), a=low[y], b=up[y], mean = X%*%beta, sd = 1)
    
    ## sample Sigma ## should this Sigma be block diagonal?
    #Sigma = riwish(nu+.5, Psi+beta%*%t(beta))
    A = Psi+solve(a)%*%beta%*%t(beta)
    Sigma = matrix(0,p*3,p*3) 
    for(k in 1:p) Sigma[1:3+3*(k-1),1:3+3*(k-1)] = riwish(nu+1, A[1:3+3*(k-1),1:3+3*(k-1)])
    
    # sample alpha
    alpha1 = rtruncnorm(1, a=max(y_star[y==1]), b=min(y_star[y==2]), mean = 0, sd = sqrt(10))
    alpha2 = rtruncnorm(1, a=max(y_star[y==2]), b=min(y_star[y==3]), mean = 0, sd = sqrt(10))
    alpha3 = rtruncnorm(1, a=max(y_star[y==3]), b=min(y_star[y==4]), mean = 0, sd = sqrt(10))
    alpha4 = rtruncnorm(1, a=max(y_star[y==4]), b=min(y_star[y==5]), mean = 0, sd = sqrt(10))
    alpha = c(alpha1,alpha2,alpha3,alpha4)
    ## sample beta
    A = solve(XtX+solve(a)%*%solve(Sigma))
    beta = mvrnorm(n = 1, mu = A%*%(t(X)%*%y_star), Sigma = A)
    ## store the results
    Beta[i,] = beta
    Alpha[i,]=alpha
    Ys[i,]=y_star
    Sigma_list[[i]] = Sigma
  }
  ## compare the performance
  #alpha_post = apply(Alpha[-(1:1000),],2,mean)
  #beta_post = apply(Beta[-(1:1000),],2,mean)
  ## for model 1
  #mu = Xt%*%beta_post
  #p11 = pnorm(alpha_post[1]-mu)
  #p12 = pnorm(alpha_post[2]-mu) - p11
  #p13 = pnorm(alpha_post[3]-mu) - pnorm(alpha_post[2]-mu)
  #p14 = pnorm(alpha_post[4]-mu) - pnorm(alpha_post[3]-mu)
  #p15 = 1 - pnorm(alpha_post[4]-mu)
  #p1 = cbind(p11,p12,p13,p14,p15)
  #y1 = apply(p1,1,which.max)
  #out = list()
  #y1 = apply(p1[,3:5],1,sum)
  #y2 = apply(p1[,4:5],1,sum)
  #out$y1 = y1
  #out$y2 = y2
  
  ## now select the important id
  #ind = unique(c(1:5,ceiling(which(abs(apply(Beta[-(1:1000),],2,mean)/apply(Beta[-(1:1000),],2,sd))>1)/3)))
  ind = unique(c(1,2,4,ceiling(which(abs(apply(Beta[-(1:1000),],2,mean)/apply(Beta[-(1:1000),],2,sd))>1)/3)))
  p2 = length(ind)
  ind = sort(c(ind*3-2,ind*3-1,ind*3))
  out = BML_binary_same_alpha(XX[,ind],yy,zz,t(XX[id_train,ind])%*%XX[id_train,ind],
                              niter = 5000,p=p2,id_train,nu=p2*3+2,psi=10,cutoff = cutoff)
  #out$ind = ind
  return(out)
}


BML_lasso_binary_diff_alpha = function(XX,yy,zz,XtX,niter = 5000,p=5,id_train,nu=3+2,psi=10,cutoff=4){
  Beta = matrix(0,niter,p*3)
  Alpha = matrix(0,niter,4*3)
  Tau = matrix(0,niter,p) ## here the tau is tau square
  gam = c() ## here gam is lambda (copied from BEL)
  Ys = matrix(0,niter,length(id_train))
  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
  z = zz[id_train]
  zt = zz[-id_train]
  
  Psi = diag(rep(psi,p*3))
  Sigma_list = {}
  
  m <- polr(as.factor(y)~X, Hess=TRUE, method="probit")
  beta = m$coefficients
  alpha = rep(m$zeta,3)
  tau = rep(1,p)
  a = diag(rep(tau,each=3))
  Tau[1,] = tau
  for(i in 1:niter){
    if(i%%5000==0) print(i)
    
    if(i>1){
      ## sample tau
      gam[i] = sqrt(2*p/sum(Tau[1:(i-1),])*(i-1))
      for(k in 1:p){
        temp = Beta[i-1,1:3+3*(k-1)]
        m_temp = (Sigma_list[[i-1]])[1:3+3*(k-1),1:3+3*(k-1)]
        d = t(temp)%*%m_temp%*%temp
        Tau[i,k] = 1/rinvgauss(n=1,m=pmax(sqrt(gam[i]^2/d),.1),s=gam[i]^2)
      }
      tau = Tau[i,]
      a = diag(rep(tau,each=3))
    }
    ## samlpe y_star
    low = c(-Inf,alpha[1:4])
    up = c(alpha[1:4],Inf)
    y_star = rtruncnorm(length(y), a=low[y], b=up[y], mean = X%*%beta, sd = 1)
    
    low = c(-Inf,alpha[1:4+4])
    up = c(alpha[1:4+4],Inf)
    y_star[z==2] = rtruncnorm(length(y[z==2]), a=low[y[z==2]], b=up[y[z==2]], mean = X[which(z==2),]%*%beta, sd = 1)
    
    low = c(-Inf,alpha[1:4+8])
    up = c(alpha[1:4+8],Inf)
    y_star[z==3] = rtruncnorm(length(y[z==3]), a=low[y[z==3]], b=up[y[z==3]], mean = X[which(z==3),]%*%beta, sd = 1)
    
    ## sample Sigma ## should this Sigma be block diagonal?
    #Sigma = riwish(nu+.5, Psi+beta%*%t(beta))
    A = Psi+solve(a)%*%beta%*%t(beta)
    Sigma = matrix(0,p*3,p*3) 
    for(k in 1:p) Sigma[1:3+3*(k-1),1:3+3*(k-1)] = riwish(nu+1, A[1:3+3*(k-1),1:3+3*(k-1)])
                  
    # sample alpha
    for(j in 1:3){
      alpha[1+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==1&z==j]), b=min(y_star[y==2&z==j]), mean = 0, sd = sqrt(10))
      alpha[2+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==2&z==j]), b=min(y_star[y==3&z==j]), mean = 0, sd = sqrt(10))
      alpha[3+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==3&z==j]), b=min(y_star[y==4&z==j]), mean = 0, sd = sqrt(10))
      alpha[4+(j-1)*4] = rtruncnorm(1, a=max(y_star[y==4&z==j]), b=min(y_star[y==5&z==j]), mean = 0, sd = sqrt(10))
    }
    
    ## sample beta
    A = solve(XtX+solve(a)%*%solve(Sigma))
    beta = mvrnorm(n = 1, mu = A%*%(t(X)%*%y_star), Sigma = A)
    
    
    ## store the results
    Beta[i,] = beta
    Alpha[i,]=alpha
    Ys[i,]=y_star
    Sigma_list[[i]] = Sigma
  }
  ## compare the performance
  alpha_post = apply(Alpha[-(1:1000),],2,mean)
  beta_post = apply(Beta[-(1:1000),],2,mean)
  
  ## for model 1
  mu = Xt%*%beta_post
  p11 = pnorm(alpha_post[1+(zt-1)*4]-mu)
  p12 = pnorm(alpha_post[2+(zt-1)*4]-mu) - p11
  p13 = pnorm(alpha_post[3+(zt-1)*4]-mu) - pnorm(alpha_post[2+(zt-1)*4]-mu)
  p14 = pnorm(alpha_post[4+(zt-1)*4]-mu) - pnorm(alpha_post[3+(zt-1)*4]-mu)
  p15 = 1 - pnorm(alpha_post[4+(zt-1)*4]-mu)
  
  p1 = cbind(p11,p12,p13,p14,p15)
  #y1 = apply(p1,1,which.max)
  out = list()
  y1 = apply(p1[,3:5],1,sum)
  y2 = apply(p1[,4:5],1,sum)
  out$y1 = y1
  out$y2 = y2
  return(out)
}

BML_lasso_binary_same_alpha_collapse_y = function(XX,yy,z,XtX,niter = 5000,p=5,id_train,nu=5*3+2,psi=10,cutoff=4){
  # BML 1: suppose alpha is the same (model 1)
  ## now start the sampling:
  yy_ori = yy
  
  yy[yy<cutoff] = 1
  yy[yy>=cutoff] = 2
  
  Beta = matrix(0,niter,p*3)
  Alpha = matrix(0,niter,1)
  Tau = matrix(0,niter,p) ## here the tau is tau square
  gam = c() ## here gam is lambda (copied from BEL)
  
  X = XX[id_train,]
  y = yy[id_train]
  yt = yy[-id_train]
  Xt = XX[-id_train,]
  
  Psi = diag(rep(psi,p*3))
  Sigma_list = {}
  
  #m <- polr(as.factor(y)~X, Hess=TRUE, method="probit")
  m = glm((y-1)~X,family = binomial(link="probit"))
  beta = m$coefficients[-1]
  alpha = m$coefficients[1]
  
  tau = rep(1,p)
  a = diag(rep(tau,each=3))
  Tau[1,] = tau
  
  for(i in 1:niter){
    if(i%%5000==0) print(i)
    
    if(i>1){
      ## sample tau
      gam[i] = sqrt(2*p/sum(Tau[1:(i-1),])*(i-1))
      for(k in 1:p){
        temp = Beta[i-1,1:3+3*(k-1)]
        m_temp = (Sigma_list[[i-1]])[1:3+3*(k-1),1:3+3*(k-1)]
        d = t(temp)%*%m_temp%*%temp
        Tau[i,k] = 1/rinvgauss(n=1,m=pmax(sqrt(gam[i]^2/d),.1),s=gam[i]^2)
      }
      tau = Tau[i,]
      a = diag(rep(tau,each=3))
    }
    
    ## samlpe y_star
    low = c(-Inf,alpha)
    up = c(alpha,Inf)
    y_star = rtruncnorm(length(y), a=low[y], b=up[y], mean = X%*%beta, sd = 1)
    
    ## sample Sigma ## should this Sigma be block diagonal?
    #Sigma = riwish(nu+.5, Psi+beta%*%t(beta))
    A = Psi+solve(a)%*%beta%*%t(beta)
    Sigma = matrix(0,p*3,p*3) 
    for(k in 1:p) Sigma[1:3+3*(k-1),1:3+3*(k-1)] = riwish(nu+1, A[1:3+3*(k-1),1:3+3*(k-1)])
    
    # sample alpha
    alpha = rtruncnorm(1, a=max(y_star[y==1]), b=min(y_star[y==2]), mean = 0, sd = sqrt(10))
    
    ## sample beta
    A = solve(XtX+solve(a)%*%solve(Sigma))
    beta = mvrnorm(n = 1, mu = A%*%(t(X)%*%y_star), Sigma = A)
    ## store the results
    Beta[i,] = beta
    Alpha[i,]=alpha
    #Ys[i,]=y_star
    Sigma_list[[i]] = Sigma
  }
  ## compare the performance
  #alpha_post = apply(Alpha[-(1:1000),],2,mean)
  #alpha_post = mean(Alpha[-(1:1000),])
  #beta_post = apply(Beta[-(1:1000),],2,mean)
  ## for model 1
  #mu = Xt%*%beta_post
  #p1 = 1-pnorm(alpha_post[1]-mu)
  
#  ind = unique(c(1:5,ceiling(which(abs(apply(Beta[-(1:1000),],2,mean)/apply(Beta[-(1:1000),],2,sd))>1)/3)))
  ind = unique(c(1,2,4,ceiling(which(abs(apply(Beta[-(1:1000),],2,mean)/apply(Beta[-(1:1000),],2,sd))>1)/3)))
  p2 = length(ind)
  ind = sort(c(ind*3-2,ind*3-1,ind*3))
  out = BML_binary_same_alpha_collapse_y(XX[,ind],yy_ori,zz,t(XX[id_train,ind])%*%XX[id_train,ind],
                              niter = 5000,p=p2,id_train,nu=p2*3+2,psi=10,cutoff = cutoff)
  #out$ind = ind
  return(out)
  #return(p1)
}


check_accu = function(y,yt,cutoff,cutoff2){
  accu = (sum(y>=cutoff&yt>=cutoff2)+sum(y<cutoff&yt<cutoff2))/length(yt)
  return(accu)
}

max_accu = function(y,yt,cutoff2){
  a = c()
  for(i in 1:100){
    a[i] = check_accu(y,yt,cutoff=i/100,cutoff2=cutoff2)
  }
  return(max(a))
}

get_auc = function(y,yt,cutoff2){
  roc <- roc(yt>=cutoff2, y)
  # Get the full AUC
  return(auc(roc))
}

