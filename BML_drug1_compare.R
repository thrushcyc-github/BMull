require("truncnorm")
require("MCMCpack")
require("MASS")
require("statmod")
require("pROC")
library(e1071)   
library(randomForest)
library(neuralnet)
library(adabag)
library(caret)

source("./Rcode/revision/BML_lib_topic_v2.R")
source("./Rcode/revision/compare_lib.R")
drug2 = read.csv("data/depression_drug1.csv")

## recode treatment to numbers
drug2$Treat_num = -1
drug2$Treat_num[drug2$Treat=="on Treatment for less than 1 month"] = 0
drug2$Treat_num[drug2$Treat=="on Treatment for 1 to 6 months"] = 1
drug2$Treat_num[drug2$Treat=="on Treatment for 6 months to less than 1 year"] = 2
drug2$Treat_num[drug2$Treat=="on Treatment for 1 to less than 2 years"] = 3
drug2$Treat_num[drug2$Treat=="on Treatment for 2 to less than 5 years"] = 4
drug2$Treat_num[drug2$Treat=="on Treatment for 5 to less than 10 years"] = 5
drug2$Treat_num[drug2$Treat=="on Treatment for 10 years or more"] = 6

## recode age to numbers
drug2$age_num = 0
drug2$age_num[drug2$Age=="3 6"] = 1
drug2$age_num[drug2$Age=="7 12"] = 2
drug2$age_num[drug2$Age=="13 18"] = 3
drug2$age_num[drug2$Age=="19 24"] = 4
drug2$age_num[drug2$Age=="25 34"] = 5
drug2$age_num[drug2$Age=="35 44"] = 6
drug2$age_num[drug2$Age=="45 54"] = 7
drug2$age_num[drug2$Age=="55 64"] = 8
drug2$age_num[drug2$Age=="65 74"] = 9
drug2$age_num[drug2$Age=="75 or over"] = 10

data = drug2[drug2$Treat_num>=0&drug2$age_num>0&drug2$Gender!="none",]
data = subset(data,select=-c(X,condition,drug,reviewer_x,text,time,useful,reviewer_y,Taker,Age,Treat,User))

data$satisfaction.fac = factor(data$satisfaction)
data$ease.fac = factor(data$ease)
data$effect.fac = factor(data$effectiveness)
data$Female = as.numeric(data$Gender=="Female")
data$age_num2 = data$age_num^2
data$Treat_num2 = data$Treat_num^2


data$topic0 = as.numeric(data$topic0_lda>.1)
data$topic1 = as.numeric(data$topic1_lda>.1)
data$topic2 = as.numeric(data$topic2_lda>.1)
data$topic3 = as.numeric(data$topic3_lda>.1)
data$topic4 = as.numeric(data$topic4_lda>.1)
data$topic5 = as.numeric(data$topic5_lda>.1)
data$topic6 = as.numeric(data$topic6_lda>.1)
data$topic7 = as.numeric(data$topic7_lda>.1)
data$topic8 = as.numeric(data$topic8_lda>.1)
data$topic9 = as.numeric(data$topic9_lda>.1)


## get the dominant topic
data$topic = apply(data[,32:41],1,which.max)

data$prior_sentiment = rep(NA,nrow(data))
for(i in 1:10){
  id = which(data$topic==i)
  temp = cumsum(data$sentiment[id])
  l = length(temp)
  data$prior_sentiment[id[-1]] = temp[1:(l-1)]/(1:(l-1))
}

data = data[!is.na(data$prior_sentiment),]

yy = c(data$satisfaction,data$ease,data$effectiveness)
n = length(yy)
n1 = nrow(data)
n2 = n1
n3 = n1
n = n1+n2+n3

p1 = 5+11
XX = matrix(0,n,p1*3)
var_list1 = c("Female","age_num","age_num2","Treat_num","Treat_num2","topic0","topic1","topic2","topic3","topic4","topic5","topic6","topic7","topic8","topic9","prior_sentiment")
var_str1 = paste(var_list1, collapse='+' )
temp = as.matrix(data[,var_list1])
#temp[,1:5] = scale(temp[,1:5])
XX[1:n1,1:p1*3-2] = temp
XX[1:n2+n1,1:p1*3-1] = temp
XX[1:n3+n1+n2,1:p1*3] = temp
colnames(XX) = c(var_list1,paste(var_list1,"2",sep=""),paste(var_list1,"3",sep=""))


p2 = 5+10
VV = matrix(0,n,p2*3)
var_list2 = c("Female","age_num","age_num2","Treat_num","Treat_num2","topic0","topic1","topic2","topic3","topic4","topic5","topic6","topic7","topic8","topic9")
var_str2 = paste(var_list2, collapse='+' )
temp = as.matrix(data[,var_list2])
VV[1:n1,1:p2*3-2] = temp
VV[1:n2+n1,1:p2*3-1] = temp
VV[1:n3+n1+n2,1:p2*3] = temp
colnames(VV) = c(var_list2,paste(var_list2,"2",sep=""),paste(var_list2,"3",sep=""))


p3 = 5
WW = matrix(0,n,p3*3)
var_list3 = c("Female","age_num","age_num2","Treat_num","Treat_num2")
var_str3 = paste(var_list3, collapse='+' )
temp = as.matrix(data[,var_list3])
#temp[,1:5] = scale(temp[,1:5])
WW[1:n1,1:p3*3-2] = temp
WW[1:n2+n1,1:p3*3-1] = temp
WW[1:n3+n1+n2,1:p3*3] = temp
colnames(WW) = c(var_list3,paste(var_list3,"2",sep=""),paste(var_list3,"3",sep=""))


zz = c(rep(1,n1),rep(2,n2),rep(3,n3))

library(foreach)
library(doParallel)

no_cores <- detectCores()
registerDoParallel(makeCluster(no_cores-1))


a = Sys.time()

out <- foreach(i=1:400,.packages = c('foreach','doParallel','MASS','statmod','truncnorm','MCMCpack','pROC',
                                     'e1071','randomForest','neuralnet','adabag'),
               .errorhandling = "remove") %dopar% 
  compare_accu_v3(n,zz,yy,XX,VV,WW,n_train=150)

Sys.time() - a

stopCluster(makeCluster(no_cores))

temp = matrix(0,13,3)
iter = length(out)
for(i in 1:iter){
  temp = temp+out[[i]]
}
print(round(temp[1:13,]/iter,3)*100)
