
evaluate=function(beta_cv,beta) {

t=100

# beta=c(0.5,0,0,-1,1.2,rep(0,5)) 

p=length(beta)


beta_nonzero=c(1,4,5)

beta_zero=c(1:p)[-beta_nonzero]


pno=3



beta_matrix=t(matrix(rep(beta,t),p,t))

 r=beta_cv


dif=abs(r-beta_matrix) 

mab=mean(colMeans(dif))

sd_mab=sd(colMeans(dif))

msb=colMeans(dif*dif)

msb=sqrt(sum(msb))/p

############################################################

ifzero=(r==0)

nonzero=mean(p-apply(ifzero,1,sum))  #求行不为零的系数的个数，再平均

sd_nonzero=sd(p-apply(ifzero,1,sum))

c_nonzero=mean(3-apply(ifzero[,beta_nonzero],1,sum))

sd_c_nonzero=sd(3-apply(ifzero[,beta_nonzero],1,sum))

c_zero=mean(apply(ifzero[,beta_zero],1,sum))

sd_c_zero=sd(apply(ifzero[,beta_zero],1,sum))

number_c=length(which(apply(ifzero[,beta_nonzero],1,sum)==0 & apply(ifzero[,beta_zero],1,sum)==(p-3)))

percent_c=number_c/t

assess=cbind(mab,sd_mab,msb,nonzero,sd_nonzero,c_nonzero,sd_c_nonzero,c_zero,sd_c_zero,percent_c)

return(assess) }