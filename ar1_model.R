## Simulation study R code for AR1 model in Table 2

#######################
#	Library call 
#######################

library(devtools)
library(tryCatchLog)
library(mvtnorm)
library(mclust)
library(clustMD)
library(QUIC)
library(MomTrunc)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(TruncatedNormal)

########################################################
#	source code for main function
########################################################

source("solver1.54.R")
sourceCpp("cal_prob.cpp") 

## sub-function for AR1 model
ar1_cor <- function(n, rho) {
exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
    (1:n - 1))
rho^exponent
}

########################################################
#	AR1 model parameter settings
########################################################

n=100 ## the number of sample size
p=c(0.5,0.5) ## mixing probability phi
grp=c(rep(1,n*p[1]),rep(2,n*p[2])) ## group membership vector

pp=11 ## the number of parameters + 1 for nominal

mu1 = rep(0, pp) ## mean vector for cluster 1
mu2 = rep(1, pp) ## mean vector for cluster 2

rho=0.0 ## from 0 to 0.8 by 0.2: as the number is higher, the correlation is higher.

## covariance matrix
sigma1=ar1_cor(pp,rho)
sigma1=sigma2=cov2cor(solve(sigma1))

########################################################
#
#		Bootstrap (100 replications)
#
#	1. clustMD method
#
#	2. regclustMD method (newly proposed)
#
########################################################

accuracy=matrix(0,100,2) 

for(B in 1:100){ # start of loop for B

#####################
## data generation
#####################

## pp-q continuous variables, q categorical variables (1 binary, 1 ordinal with 3 levels, 1 nominal)

set.seed(B*10)

new.B=B
nom.length=0

while(nom.length<3){

dat=rbind(rmvnorm(n*p[1],mean=mu1,sigma=sigma1),
rmvnorm(n*p[2],mean=mu2,sigma=sigma2))


	ordi.z1=dat[,8]
	ordi1<-NULL
	cut1=median(ordi.z1)
	ordi1[ordi.z1<=cut1]=1;ordi1[ordi.z1>cut1]=2

	ordi.z2=dat[,9]
	ordi2=NULL
	cut2=quantile(ordi.z2,c(0.25,0.75))
	ordi2[ordi.z2<=cut2[1]]=1
	ordi2[ordi.z2>cut2[1] & ordi.z2<=cut2[2]]=2
	ordi2[ordi.z2>cut2[2]]=3

	nomi<-NULL
	nomi1 = apply(dat[,10:11],1,which.max) + 1
	nomi1[apply(dat[,10:11],1,max)<=0] = 1     

	ordi=cbind(ordi1,ordi2,nomi1)

set.seed((new.B+100)*10)
new.B=new.B+1

nom.length=length(table(nomi1))
}

data.y=cbind(dat[,1:7],ordi)
head(data.y)

#####################
## clustMD method
#####################

md.res=try(clustMD(X=data.y, G=2, CnsIndx=pp-4, OrdIndx=pp-2, Nnorms=20000, 
	MaxIter=100, model="VVI", startCL="kmeans", stop.tol=0.001))

if(class(md.res)!= "try-error"){
	md.tab=table(md.res$cl,grp)
	if(sum(md.tab)==0){
		accuracy[B,1]=NA
	}else{
		accuracy[B,1]=max(sum(diag(md.tab))/n, 1-sum(diag(md.tab))/n) # accuracy
	}
}else{
	accuracy[B,1]=NA
}

#####################
## regclustMD method
#####################

is.debug=TRUE
reg.res=tryLog(regClustMD.BIC(X=data.y, G=2, CnsIndx=pp-4, OrdIndx=pp-2, Nnorms=20000,
	MaxIter=100, store.params = TRUE, scale = TRUE, stop.tol=0.001,
	verbose.outer=TRUE, verbose.inner=FALSE), execution.context.msg = B, 	write.error.dump.file = TRUE)

if(class(reg.res)!= "try-error"){
	reg.tab= table(factor(reg.res$ind, levels=c(1,2)), grp) 
	if(sum(reg.tab)==0){
		accuracy[B,2]=NA
	}else{
		accuracy[B,2]=max(sum(diag(reg.tab))/n, 1-sum(diag(reg.tab))/n) # accuracy
	}
}else{
	accuracy[B,2]=NA
}

} # end of loop for B


#####################
## 	result
#####################

apply(accuracy,2,mean,na.rm=T) # accuracy mean
apply(accuracy,2,sd,na.rm=T)/sqrt(100) # accuracy standard error
apply(apply(accuracy,2,is.na),2,mean)

