library(glmnet)

library(ncvreg)

library(SIS)


# setwd("C:/Users/Administrator/Desktop/LA_research/L0Regression/R_code/new_simulation")

rm( list = ls ( all = TRUE))

# source('logisticL0.R')  

source('evaluate.R')
source('dev.R')

#--------initializing-------------------

set.seed(321)
seeds=round(runif(100)*100000)

# beta=c(2,0,0,-4,5,rep(0,5)) 

beta=c(0.5,0,0,-1,1.2,rep(0,5)) 

p=length(beta)

var_no=seq(1,p)

n=300

iter=100

# lambda=2

e=1e-6

e0=1e-10


#-------------saved variables---------------------------


beta_SCAD_ncv=matrix(1000,iter,p)

beta_SCAD_SIS=matrix(0,iter,p)


beta_lasso=matrix(1000,iter,p)

beta_alasso=matrix(1000,iter,p)



#-------------------------------------------------------------

for (t in 1:iter){
  
  
  #--------------data generating-------------------------------
  
  
  set.seed(seeds[t])
  
  x_all = rnorm(n*p)   
  
  X=matrix(x_all,n,p) 
  
  #simulation design: Y~Bernoulli{pr(t(X)*Beta) where pr(u)=exp(u)/(1-exp(u)),      
  
  z = X%*%beta       
  
  pr = 1/(1+exp(-z))         
  
  Y = rbinom(n,1,pr) 
  
  #-----initial beta-------------------------------
  
  cvob0=cv.glmnet(X,Y,family="binomial", alpha = 0, nfolds=10)
  
  
  # print(cbind(t,cvob0$lambda.min))
  
  initial=glmnet(X, Y, family="binomial", alpha = 0,lambda=cvob0$lambda.min)
  
  
  
  
  #---------------------SCAD,ncv------------------------------------------------
  cv.scad=cv.ncvreg(X, Y, family="binomial",penalty="SCAD", alpha=1,nfolds=10)
  
  l.min=cv.scad$lambda.min
  
  scad=ncvreg(X, Y, family="binomial",penalty="SCAD", alpha=1,lambda=l.min)
  
  beta_SCAD_ncv[t,]=scad$beta[2:(p+1)]
  
  #---------------------SCAD,SIS------------------------------------------------
  
  
  scad.sis=SIS(X, Y, family="binomial",penalty="SCAD", standardize=FALSE)
  
  coef=scad.sis$coef.est
  
  coef_1=coef[-1]
  
  b=as.numeric(substring(names(coef_1),2))
  
  beta_SCAD_SIS[t,b]=coef_1
  
  #--------------------MCP, ncv------------------------------------------------
  cv.mcp=cv.ncvreg(X, Y, family="binomial",penalty="MCP",nfolds=10)
  
  l.min=cv.mcp$lambda.min
  
  mcp=ncvreg(X, Y, family="binomial",penalty="MCP", lambda=l.min)
  
  beta_MCP_ncv[t,]=mcp$beta[2:(p+1)]
  
  #--------------------Lasso------------------------------------------------
  
  cvob1=cv.glmnet(X,Y,family="binomial", alpha = 1, nfolds=10)
  
  l1=glmnet(X, Y, family="binomial", alpha = 1,lambda=cvob1$lambda.min)
  
  b_1=l1$beta
  
  # no_1=which(abs(b_1)>e0)
  
  beta_lasso[t,]=t(as.matrix(b_1))
  
  #-----------------alasso--------------------------------------------------
  
  alasso.cv = cv.glmnet(X,Y,family="binomial",alpha=1,penalty.factor=1/abs(initial$beta))
  
  # alcoef<-coef(alasso.cv,s="lambda.min")
  
  alasso = glmnet(X,Y,family="binomial",lambda=alasso.cv$lambda.min,alpha=1,penalty.factor=1/abs(initial$beta))
  
  beta_alasso[t,]=t(as.matrix(alasso$beta))
  
  
}




write.csv(beta_SCAD_ncv,file="beta_SCAD_n100p10.csv")

evaluate(beta_SCAD_ncv,beta)

write.csv(beta_SCAD_SIS,file="beta_SCAD_SIS_n100p10.csv")

evaluate(beta_SCAD_SIS,beta)

write.csv(beta_MCP_ncv,file="beta_MCP_n100p10.csv")

evaluate(beta_MCP_ncv,beta)

write.csv(beta_lasso,file="beta_lasso_n100p10.csv")

evaluate(beta_lasso,beta)

write.csv(beta_alasso,file="beta_alasso_n100p10.csv")

evaluate(beta_alasso,beta)




