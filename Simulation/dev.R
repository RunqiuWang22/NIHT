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