polr_compare = function(y,X,Xt){
  m <- polr(as.factor(y)~X, Hess=TRUE, method="probit")
  beta = m$coefficients
  alpha = m$zeta
  mu4 = Xt%*%beta
  p41 = pnorm(alpha[1]-mu4)
  p42 = pnorm(alpha[2]-mu4) - p41
  p43 = pnorm(alpha[3]-mu4) - pnorm(alpha[2]-mu4)
  p44 = pnorm(alpha[4]-mu4) - pnorm(alpha[3]-mu4)
  p45 = 1 - pnorm(alpha[4]-mu4)
  
  pp4 = cbind(p41,p42,p43,p44,p45)
  
  out = list()
  y1 = apply(pp4[,3:5],1,sum)
  y2 = apply(pp4[,4:5],1,sum)
  out$y1 = y1
  out$y2 = y2
  return(out)
}

polr_compare_diff_alpha = function(y,X,Xt,z_train,zt,p){
  out = list()
  n = nrow(Xt)
  out$y1 = rep(0,n)
  out$y2 = rep(0,n)
  for(i in 1:3){
    m <- polr(as.factor(y[z_train==i])~X[z_train==i,1:p*3-3+i], Hess=TRUE, method="probit")
    beta = m$coefficients
    alpha = m$zeta
    mu4 = Xt[zt==i,1:p*3-3+i]%*%beta
    p41 = pnorm(alpha[1]-mu4)
    p42 = pnorm(alpha[2]-mu4) - p41
    p43 = pnorm(alpha[3]-mu4) - pnorm(alpha[2]-mu4)
    p44 = pnorm(alpha[4]-mu4) - pnorm(alpha[3]-mu4)
    p45 = 1 - pnorm(alpha[4]-mu4)
    
    pp4 = cbind(p41,p42,p43,p44,p45)
    y1 = apply(pp4[,3:5],1,sum)
    y2 = apply(pp4[,4:5],1,sum)
    out$y1[zt==i] = y1
    out$y2[zt==i] = y2
  }
  return(out)
}


compare_accu_v3 = function(n,zz,yy,XX,VV,WW,n_train=300){
  #now only keep the first part. 
  mm = 13
  #res = matrix(0,mm*4,12)
  res = matrix(0,mm,3)
  str = c("BML_topic","BML_non","probit_topic","probit_non","svm_topic","svm_non","rf_topic","rf_non",
          "nn_topic","nn_non","ada_topic","ada_non","BML_lasso_topic")
  #rownames(res) = c(str, paste(str,"_sat",sep=""), paste(str,"_ease",sep=""), paste(str,"_eff",sep=""))
  rownames(res) = str
  id_train = sample(1:n,n_train)
  id_test = setdiff(1:n,id_train)
  
  #id_train = sample(1:(n/3),n_train/3)
  #id_train = c(id_train,id_train+(n/3),id_train+(n/3*2))
  #id_test = setdiff(1:n,id_train)
  
  z_train = zz[id_train]
  zt = zz[-id_train]
  y = yy[id_train]
  yt = yy[id_test]
  
  X_train = XX[id_train,]
  X = X_train
  Xt = XX[id_test,]
  XtX = t(X_train)%*%X_train
  
  V_train = VV[id_train,]
  V = V_train
  Vt = VV[id_test,]
  VtV = t(V_train)%*%V_train
  
  W_train = WW[id_train,]
  W = W_train
  Wt = WW[id_test,]
  WtW = t(W_train)%*%W_train
  
  data1_train = data.frame(y,X)
  data1_test = data.frame(yt,Xt)
  data3_train = data.frame(y,W)
  data3_test = data.frame(yt,Wt)
  
  p1 = ncol(XX)/3
  p3 = ncol(WW)/3
  
  ###############################
  ##
  ## first part: 5 class
  ##
  ###############################
  
  ## Bayesian
  out = BML_binary_same_alpha(XX,yy,zz,XtX,niter = 5000,p=p1,id_train,nu=p1*3+2,psi=10,cutoff = cutoff)
  y1 = out$y1;y21 = out$y2
  out = BML_binary_same_alpha(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = cutoff)
  y2 = out$y1;y22 = out$y2
  
  ## Frequentist
  # for model 4 
  out = polr_compare(y,X,Xt)
  y3 = out$y1;y23 = out$y2
  # for model 6 
  out = polr_compare(y,W,Wt)
  y4 = out$y1;y24 = out$y2
  
  ## competitor: svm, rf, nn
  svm1 <- svm(as.factor(y)~., data=data1_train, probability = TRUE)
  pred = predict(svm1,data1_test,probability = TRUE)
  pred = attr(pred, "probabilities")
  y5 = 1- pred[,"1"]-pred[,"2"]
  #y25 = 1- pred[,"1"]-pred[,"2"]-pred[,"3"]

  svm1 <- svm(as.factor(y)~., data=data3_train, probability = TRUE)
  pred = predict(svm1,data3_test,probability = TRUE)
  pred = attr(pred, "probabilities")
  y6 = 1- pred[,"1"]-pred[,"2"]
  #y26 = 1- pred[,"1"]-pred[,"2"]-pred[,"3"]
  
  rf <- randomForest(as.factor(y) ~ .,data=data1_train)
  y7 = apply(predict(rf, newdata=data1_test,type="prob")[,3:5],1,sum)
  #y27 = apply(predict(rf, newdata=data1_test,type="prob")[,4:5],1,sum)

  rf <- randomForest(as.factor(y) ~ .,data=data3_train)
  y8 = apply(predict(rf, newdata=data3_test,type="prob")[,3:5],1,sum)
  #y28 = apply(predict(rf, newdata=data3_test,type="prob")[,4:5],1,sum)

  print("nn")
  nn <- neuralnet(as.factor(y) ~ .,data=data1_train, hidden=c(1), linear.output=FALSE, threshold=0.01)
  pred <- compute(nn, data1_test)
  y9 = apply(pred$net.result[,3:5],1,sum)
  #y29 = apply(pred$net.result[,4:5],1,sum)

  nn <- neuralnet(as.factor(y) ~ .,data=data3_train, hidden=c(1), linear.output=FALSE, threshold=0.01)
  pred <- compute(nn, data3_test)
  y10 = apply(pred$net.result[,3:5],1,sum)
  #y30 = apply(pred$net.result[,4:5],1,sum)
  
  print("ada")
  data1_train$y = as.factor(data1_train$y)
  data1_test$y = as.factor(data1_test$y)
  ada = boosting(y ~ .,data=data1_train, boos=TRUE, mfinal=50)
  y11 = predict(ada, data1_test)
  y11 = apply(y11$prob[,3:5],1,sum)
  
  data3_train$y = as.factor(data3_train$y)
  data3_test$y = as.factor(data3_test$y)
  ada = boosting(y ~ .,data=data3_train, boos=TRUE, mfinal=50)
  y12 = predict(ada, data3_test)
  y12 = apply(y12$prob[,3:5],1,sum)
  
  print("BMul")
  out = BML_lasso_binary_same_alpha(XX,yy,zz,XtX,niter = 5000,p=p1,id_train,nu=p1*3+2,psi=10,cutoff = cutoff)
  y13 = out$y1;#y29 = out$y2
  #out = BML_lasso_binary_same_alpha(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = cutoff)
  #y10 = out$y1;y30 = out$y2
  
  # removed columns from 10, only keep the first 3
  ###############################
  ##
  ## second part: 2 class
  ##
  ###############################
  
  ## Bayesian
#  out = BML_binary_same_alpha_collapse_y(XX,yy,zz,XtX,niter = 5000,p=p1,id_train,nu=p1*3+2,psi=10,cutoff = 3)
#  y11 = out
#  out = BML_binary_same_alpha_collapse_y(XX,yy,zz,XtX,niter = 5000,p=p1,id_train,nu=p1*3+2,psi=10,cutoff = 4)
#  y31 = out
#  out = BML_binary_same_alpha_collapse_y(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = 3)
#  y12 = out
#  out = BML_binary_same_alpha_collapse_y(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = 4)
#  y32 = out
  
  ## Frequentist probit
  # for model 4 
#  y_3 = I(y>=3)
#  y_4 = I(y>=4)
#  A = data.frame(X,y_3)
#  B = data.frame(X,y_4)
#  C = data.frame(Xt,yt)
#  m = glm(y_3~.,data = A,family = binomial(link="probit"))
#  y13 = predict.glm(m,C,type="response")
#  m = glm(y_4~.,data = B,family = binomial(link="probit"))
#  y33 = predict.glm(m,C,type="response")
  
  # for model 6 
#  A = data.frame(W,y_3)
#  B = data.frame(W,y_4)
#  C = data.frame(Wt,yt)
#  m = glm(y_3~.,data = A,family = binomial(link="probit"))
#  y14 = predict.glm(m,C,type="response")
#  m = glm(y_4~.,data = B,family = binomial(link="probit"))
#  y34 = predict.glm(m,C,type="response")
  
  ## competitor: svm, rf, nn
#  svm1 <- svm(I((y>=3)*1)~., data=data1_train,method="C-classification", kernal="radial")
#  y15 = predict(svm1,data1_test,probability = TRUE)
#  svm1 <- svm(I((y>=4)*1)~., data=data1_train,method="C-classification", kernal="radial")
#  y35 = predict(svm1,data1_test,probability = TRUE)
#  
#  svm1 <- svm(I((y>=3)*1)~., data=data3_train,method="C-classification", kernal="radial")
#  y16 = predict(svm1,data3_test,probability = TRUE)
#  svm1 <- svm(I((y>=4)*1)~., data=data3_train,method="C-classification", kernal="radial")
#  y36 = predict(svm1,data3_test,probability = TRUE)
  
  
#  rf <- randomForest(as.factor(y>=3) ~ .,data=data1_train)
#  y17 = predict(rf, newdata=data1_test,type="prob")[,2]
#  rf <- randomForest(as.factor(y>=4) ~ .,data=data1_train)
#  y37 = predict(rf, newdata=data1_test,type="prob")[,2]
  
#  rf <- randomForest(as.factor(y>=3) ~ .,data=data3_train)
#  y18 = predict(rf, newdata=data3_test,type="prob")[,2]
#  rf <- randomForest(as.factor(y>=4) ~ .,data=data3_train)
#  y38 = predict(rf, newdata=data3_test,type="prob")[,2]
  
#  nn <- neuralnet(as.factor(y>=3) ~ .,data=data1_train, hidden=c(5), linear.output=FALSE, threshold=0.01)
#  pred <- compute(nn, data1_test)
#  y19 = pred$net.result[,2]
#  nn <- neuralnet(as.factor(y>=4) ~ .,data=data1_train, hidden=c(5), linear.output=FALSE, threshold=0.01)
#  pred <- compute(nn, data1_test)
#  y39 = pred$net.result[,2]

#  nn <- neuralnet(as.factor(y>=3) ~ .,data=data3_train, hidden=c(5), linear.output=FALSE, threshold=0.01)
#  pred <- compute(nn, data3_test)
#  y20 = pred$net.result[,2]
#  nn <- neuralnet(as.factor(y>=4) ~ .,data=data3_train, hidden=c(5), linear.output=FALSE, threshold=0.01)
#  pred <- compute(nn, data3_test)
#  y40 = pred$net.result[,2]

  ## Bayesian lasso
#  out = BML_lasso_binary_same_alpha_collapse_y(XX,yy,zz,XtX,niter = 5000,p=p1,id_train,nu=p1*3+2,psi=10,cutoff = 3)
#  y19 = out
#  out = BML_lasso_binary_same_alpha_collapse_y(XX,yy,zz,XtX,niter = 5000,p=p1,id_train,nu=p1*3+2,psi=10,cutoff = 4)
#  y39 = out
  #out = BML_lasso_binary_same_alpha_collapse_y(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = 3)
  #y20 = out
  #out = BML_lasso_binary_same_alpha_collapse_y(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = 4)
  #y40 = out
  
  y_res = data.frame(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13)#y10)
                     #y11,y12,y13,y14,y15,y16,y17,y18,y19,#y20,
                     #y21,y22,y23,y24,y25,y26,y27,y28,y29,#y30,
                     #y31,y32,y33,y34,y35,y36,y37,y38,y39)
  print("done")
  ## overall
  for(k in 1:mm){  
    
    name = paste("y",k,sep="")
    res[k,1:3] = c(check_accu(y_res[,name],yt,cutoff=.5,cutoff2=3),
                   max_accu(y_res[,name],yt,cutoff2=3),get_auc(y_res[,name],yt,cutoff2=3))
  }
  #for(k in 1:mm+10){  
  #  name = paste("y",k,sep="")
  #  res[k-10,1:3+3] = c(check_accu(y_res[,name],yt,cutoff=.5,cutoff2=3),
  #                     max_accu(y_res[,name],yt,cutoff2=3),get_auc(y_res[,name],yt,cutoff2=3))
  #}
  #for(k in 1:mm+20){  
  #  name = paste("y",k,sep="")
  #  res[k-20,1:3+3*2] = c(check_accu(y_res[,name],yt,cutoff=.5,cutoff2=4),
  #                         max_accu(y_res[,name],yt,cutoff2=4),get_auc(y_res[,name],yt,cutoff2=4))
  #}
  #for(k in 1:mm+30){  
  #  name = paste("y",k,sep="")
  #  res[k-30,1:3+3*3] = c(check_accu(y_res[,name],yt,cutoff=.5,cutoff2=4),
  #                         max_accu(y_res[,name],yt,cutoff2=4),get_auc(y_res[,name],yt,cutoff2=4))
  #}
  
  #For evaluation of each task 
  #for(j in 1:3){
  #  for(k in 1:mm){  
  #    name = paste("y",k,sep="")
  #    res[k+mm*j,1:3] = c(check_accu(y_res[zt==j,name],yt[zt==j],cutoff=.5,cutoff2=3),
  #                   max_accu(y_res[zt==j,name],yt[zt==j],cutoff2=3),get_auc(y_res[zt==j,name],yt[zt==j],cutoff2=3))
  #  }
    #for(k in 1:mm+10){  
    #  name = paste("y",k,sep="")
    #  res[k-10+mm*j,1:3+3] = c(check_accu(y_res[zt==j,name],yt[zt==j],cutoff=.5,cutoff2=3),
    #                      max_accu(y_res[zt==j,name],yt[zt==j],cutoff2=3),get_auc(y_res[zt==j,name],yt[zt==j],cutoff2=3))
    #}
    #for(k in 1:mm+20){  
    #  name = paste("y",k,sep="")
    #  res[k-20+mm*j,1:3+3*2] = c(check_accu(y_res[zt==j,name],yt[zt==j],cutoff=.5,cutoff2=4),
    #                         max_accu(y_res[zt==j,name],yt[zt==j],cutoff2=4),get_auc(y_res[zt==j,name],yt[zt==j],cutoff2=4))
    #}
    #for(k in 1:mm+30){  
    #  name = paste("y",k,sep="")
    #  res[k-30+mm*j,1:3+3*3] = c(check_accu(y_res[zt==j,name],yt[zt==j],cutoff=.5,cutoff2=4),
    #                        max_accu(y_res[zt==j,name],yt[zt==j],cutoff2=4),get_auc(y_res[zt==j,name],yt[zt==j],cutoff2=4))
    #}
  #}
  return(res)
} 


compare_accu_diff_alpha_beta = function(n,zz,yy,WW,n_train=300,p3){
  res = matrix(0,4*4,6)
  str = c("BML_same_alpha","BML_diff_alpha","probit_same_alpha","probit_diff_alpha")
  rownames(res) = c(str, paste(str,"_sat",sep=""), paste(str,"_ease",sep=""), paste(str,"_eff",sep=""))
  id_train = sample(1:n,n_train)
  id_test = setdiff(1:n,id_train)
  
  z_train = zz[id_train]
  zt = zz[-id_train]
  y = yy[id_train]
  yt = yy[id_test]

  W_train = WW[id_train,]
  W = W_train
  Wt = WW[id_test,]
  WtW = t(W_train)%*%W_train
  
  data3_train = data.frame(y,W)
  data3_test = data.frame(yt,Wt)
  
  ###############################
  ##
  ## first part: 5 class
  ##
  ###############################
  
  ## Bayesian
  out = BML_binary_same_alpha(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = cutoff)
  y1 = out$y1
  out = BML_binary_diff_alpha(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = cutoff)
  y2 = out$y1
  

  out = BML_binary_same_alpha_ind_sig(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = cutoff)
  y3 = out$y1
  out = BML_binary_diff_alpha_ind_sig(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = cutoff)
  y4 = out$y1
  
  ###############################
  ##
  ## second part: 2 class
  ##
  ###############################
  
  ## Bayesian
  out = BML_binary_same_alpha_collapse_y(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = 3)
  y5 = out
  out = BML_binary_diff_alpha_collapse_y(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = 3)
  y6 = out
  
  
  out = BML_binary_same_alpha_ind_sig_collapse_y(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = 3)
  y7 = out
  out = BML_binary_diff_alpha_ind_sig_collapse_y(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = 3)
  y8 = out
  

  y_res = data.frame(y1,y2,y3,y4,y5,y6,y7,y8)
  ## overall
  for(k in 1:4){  
    name = paste("y",k,sep="")
    res[k,1:3] = c(check_accu(y_res[,name],yt,cutoff=.5,cutoff2=3),
                   max_accu(y_res[,name],yt,cutoff2=3),get_auc(y_res[,name],yt,cutoff2=3))
  }
  for(k in 1:4+4){  
    name = paste("y",k,sep="")
    res[k-4,1:3+3] = c(check_accu(y_res[,name],yt,cutoff=.5,cutoff2=3),
                   max_accu(y_res[,name],yt,cutoff2=3),get_auc(y_res[,name],yt,cutoff2=3))
  }
  

  for(j in 1:3){
    for(k in 1:4){  
      name = paste("y",k,sep="")
      res[k+4*j,1:3] = c(check_accu(y_res[zt==j,name],yt[zt==j],cutoff=.5,cutoff2=3),
                         max_accu(y_res[zt==j,name],yt[zt==j],cutoff2=3),get_auc(y_res[zt==j,name],yt[zt==j],cutoff2=3))
    }
    for(k in 1:4+4){  
      name = paste("y",k,sep="")
      res[(k-4)+4*j,1:3+3] = c(check_accu(y_res[zt==j,name],yt[zt==j],cutoff=.5,cutoff2=3),
                         max_accu(y_res[zt==j,name],yt[zt==j],cutoff2=3),get_auc(y_res[zt==j,name],yt[zt==j],cutoff2=3))
    }
  }
  return(res)
} 


compare_accu_diff_alpha = function(n,zz,yy,WW,n_train=300,p3){
  res = matrix(0,4*4,6)
  str = c("BML_same_alpha","BML_diff_alpha","probit_same_alpha","probit_diff_alpha")
  rownames(res) = c(str, paste(str,"_sat",sep=""), paste(str,"_ease",sep=""), paste(str,"_eff",sep=""))
  id_train = sample(1:n,n_train)
  id_test = setdiff(1:n,id_train)
  
  z_train = zz[id_train]
  zt = zz[-id_train]
  y = yy[id_train]
  yt = yy[id_test]
  
  W_train = WW[id_train,]
  W = W_train
  Wt = WW[id_test,]
  WtW = t(W_train)%*%W_train
  
  data3_train = data.frame(y,W)
  data3_test = data.frame(yt,Wt)
  
  ###############################
  ##
  ## first part: 5 class
  ##
  ###############################
  
  ## Bayesian
  out = BML_binary_same_alpha(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = cutoff)
  y1 = out$y1
  out = BML_binary_diff_alpha(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = cutoff)
  y2 = out$y1
  
  ## Frequentist
  # for model 6 
  out = polr_compare(y,W,Wt)
  y3 = out$y1;
  # for model 6 
  out = polr_compare_diff_alpha(y,W,Wt,z_train,zt,p=p3)
  y4 = out$y1;
  
  ###############################
  ##
  ## second part: 2 class
  ##
  ###############################
  
  ## Bayesian
  out = BML_binary_same_alpha_collapse_y(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = 3)
  y5 = out
  out = BML_binary_diff_alpha_collapse_y(WW,yy,zz,WtW,niter = 5000,p=p3,id_train,nu=p3*3+2,psi=10,cutoff = 3)
  y6 = out
  
  ## Frequentist probit
  
  y_3 = I(y>=3)
  
  # for model 6 
  A = data.frame(W,y_3)
  C = data.frame(Wt,yt)
  y8 = rep(0,length(yt))
  
  m = glm(y_3~.,data = A,family = binomial(link="probit"))
  y7 = predict.glm(m,C,type="response")
  
  for(i in 1:3){
    m = glm(y_3~.,data = A[z_train==i,c(1:p3*3-3+i,p3*3+1)],family = binomial(link="probit"))
    y8[zt==i] = predict.glm(m,C[zt==i,c(1:p3*3-3+i,p3*3+1)],type="response")
    
  }
  
  y_res = data.frame(y1,y2,y3,y4,y5,y6,y7,y8)
  ## overall
  for(k in 1:4){  
    name = paste("y",k,sep="")
    res[k,1:3] = c(check_accu(y_res[,name],yt,cutoff=.5,cutoff2=3),
                   max_accu(y_res[,name],yt,cutoff2=3),get_auc(y_res[,name],yt,cutoff2=3))
  }
  for(k in 1:4+4){  
    name = paste("y",k,sep="")
    res[k-4,1:3+3] = c(check_accu(y_res[,name],yt,cutoff=.5,cutoff2=3),
                       max_accu(y_res[,name],yt,cutoff2=3),get_auc(y_res[,name],yt,cutoff2=3))
  }
  
  
  for(j in 1:3){
    for(k in 1:4){  
      name = paste("y",k,sep="")
      res[k+4*j,1:3] = c(check_accu(y_res[zt==j,name],yt[zt==j],cutoff=.5,cutoff2=3),
                         max_accu(y_res[zt==j,name],yt[zt==j],cutoff2=3),get_auc(y_res[zt==j,name],yt[zt==j],cutoff2=3))
    }
    for(k in 1:4+4){  
      name = paste("y",k,sep="")
      res[(k-4)+4*j,1:3+3] = c(check_accu(y_res[zt==j,name],yt[zt==j],cutoff=.5,cutoff2=3),
                               max_accu(y_res[zt==j,name],yt[zt==j],cutoff2=3),get_auc(y_res[zt==j,name],yt[zt==j],cutoff2=3))
    }
  }
  return(res)
} 
