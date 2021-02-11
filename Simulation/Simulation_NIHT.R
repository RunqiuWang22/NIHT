library(glmnet)

#############------------------------------------------IHT_L function-----------------------------------------------------------------------------------

IHT_L=function(x,y, Beta0, K){
  
  n=dim(x)[1]
  
  p=dim(x)[2]
  
  G=rep(0,p)
  
  h=rep(0,p)
  
  Beta=rep(0,p)
  
  NewBeta=rep(0,p)
  
  OddBeta=rep(0,p)
  
  New_Beta=rep(0,p)
  
  mu=rep(0,1000)
  
  w=rep(0,1000)
  
  c=0.001
  
  k=3
  
  thres=10^(-6)
  
  j=0
  
  h[order(t(x)%*%y,decreasing=TRUE)[1:K]]=1
  
  T0=order(h*(t(x)%*%y),decreasing=TRUE)[1:K]
  
  #--------------------------main algorithm--------------------------
  
  Beta=Beta0
  
  NT=T0
  
  OddBeta=Beta
  
  while(sum(abs(New_Beta-Beta))>thres&j<1000){
    
    
    
    j=j+1#iteration time
    
    T=order(Beta,decreasing=TRUE)[1:K]
    
    Beta=OddBeta
    
    mu[j]=sum(Beta^2)/sum((x%*%Beta)^2)
    
    z=x%*%Beta
    
    d=rep(1,n)/(1+exp(-z))
    
    G=t(x)%*%(d-y) #gradient
    
    NewBeta=Beta-mu[j]*G
    
    h[order(abs(NewBeta),decreasing=TRUE)[1:K]]=1
    
    h[-order(abs(NewBeta),decreasing=TRUE)[1:K]]=0
    
    NewBeta=h*NewBeta
    
    if (setequal(NT,T)){New_Beta=NewBeta}
    
    
    
    else{
      
      w[j]=(1-c)*sum((NewBeta-Beta)^2)/sum((x%*%(NewBeta-Beta))^2)
      
      if (is.na(w[j])==TRUE) { New_Beta=NewBeta &break}
      
      
      
      ####whether this algorithm convergence
      
      while (mu[j]>w[j]){
        
        mu[j]=mu[j]/(k*(1-c))
        
        NewBeta=Beta-mu[j]*G
        
        h[order(abs(NewBeta),decreasing=TRUE)[1:K]]=1
        
        h[-order(abs(NewBeta),decreasing=TRUE)[1:K]]=0
        
        NewBeta=h*NewBeta
        
        w[j]=(1-c)*sum((NewBeta-Beta)^2)/sum((x%*%(NewBeta-Beta))^2)
        
        if (is.na(w[j])==TRUE) {break}
        
      }
      
      New_Beta=NewBeta
      
    }
    
    NT=order(New_Beta,decreasing=TRUE)[1:K]
    
    OddBeta=New_Beta
    
  }
  
  return (New_Beta)
  
}



##################------------------

dev=function(x,y,beta)
  
{
  
  # beta is a vertical vector of dim p*1;x is the design matrix of dim n*p
  
  # p is the probabilty for y=1 in logistic regression
  
  # return pLogit is a vector of dim n*1
  
  p=exp(x%*%beta)/(1+exp(x%*%beta))
  
  lh=y*log(p)+(1-y)*log(1-p)
  
  lh[which(is.na(lh)==TRUE)]=0
  
  dev=-2*sum(lh)
  
  return(dev)
  
}

set.seed(321)
seeds=round(runif(100)*100000)
beta=c(0.5,0,0,-1,1.2,rep(0, 995))

p=length(beta)

n=300

iter=100

K=rep(0,iter)

mab=rep(0,p)

msb=rep(0,p)

#----------------------saved Beta---------------

Beta_IHT=matrix(0,iter,p)

##----------------

for(j in 1:1000) {
  
  for (t in 1:iter){
    
    ###-------------generate data-------------------
    
    set.seed(seeds[t])
    
    x_all = rnorm(n*p)  
    
    X=matrix(x_all,n,p)
    
    #simulation design: Y~Bernoulli{pr(t(X)*Beta) where pr(u)=exp(u)/(1-exp(u)),     
    
    z = X%*%beta      
    
    pr = 1/(1+exp(-z))        
    
    Y = rbinom(n,1,pr)
    
    ####----initial beta------------------------------------------------------
    
    cvob0=cv.glmnet(X,Y,family="binomial", alpha = 0, nfolds=10)
    
    # print(cbind(t,cvob0$lambda.min))
    
    initial=glmnet(X, Y, family="binomial", alpha = 0,lambda=cvob0$lambda.min)
    
    #---------------K IHT------------------------------
    
    Beta_IHT[t,]= IHT_L(X,Y, as.vector(initial$beta),j)}
  
  mab[j]=evaluate(Beta_IHT,beta)[1]
  
  msb[j]=evaluate(Beta_IHT,beta)[2]
  
  
}
