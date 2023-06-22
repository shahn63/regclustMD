#############################################################################################################
# PROGRAM 
# Code for real data example : Prostate Cancer data
#############################################################################################################
rm(list=ls())
setwd("")

require(mvtnorm);require(mclust);require(clustMD);require(QUIC);require(MomTrunc)
require(Matrix);require(Rcpp);require(RcppArmadillo);require(TruncatedNormal); require(psych)
source("solver1.54.R")
sourceCpp("cal_prob.cpp") 
source("cluster_validation_final.R")

data(Byar)

## data transformation referring to clustMD
Byar$Size.of.primary.tumour = sqrt(Byar$Size.of.primary.tumour)
Byar$Serum.prostatic.acid.phosphatase = log(Byar$Serum.prostatic.acid.phosphatase)
Y = as.matrix(Byar[, c(1, 2, 5, 6, 8, 9, 10, 11, 3, 4, 12, 7)])
Y[, 9:12] = Y[, 9:12] + 1
Y[, 1:8] = scale(Y[, 1:8])
Yekg = rep(NA, nrow(Y))
Yekg[Y[,12]==1] = 1
Yekg[(Y[,12]==2)|(Y[,12]==3)|(Y[,12]==4)] = 2
Yekg[(Y[,12]==5)|(Y[,12]==6)|(Y[,12]==7)] = 3
Y[, 12] = Yekg

Y1 = Y
Y1 = data.frame(Y1)
Y1[,9] = as.factor(Y1[,9])
Y1[,10] = as.factor(Y1[,10])
Y1[,11] = as.factor(Y1[,11])
Y1[,12] = as.factor(Y1[,12])

###############################################################################
## Correlation matrix estimation
corr1 = mixedCor(Y, c=1:9, p=12, d=10:11)
cor_res1 = round(corr1$rho,3)

Byar2 = as.matrix(Byar[, c(1, 2, 5, 6, 8, 9, 10, 11, 3, 4, 12, 7)])
corr2 = mixedCor(Byar2, c=1:9, p=12, d=10:11)
cor_res2 = round(corr2$rho,3)

###############################################################################
## Inverse matrix estimation
inv_cov_Byar_raw = solve(cov(Byar2))
round(inv_cov_Byar_raw, 3)
inv_cov_Byar_trasformed = solve(cov(Y))
round(inv_cov_Byar_trasformed, 3)

lambda = lambdaest(Y1) ;lambda       

tval.md = tval.reg = c()
G = 2:5

for (i in 1:length(G)){
set.seed(123)

### code for clustMD
is.debug = TRUE
res = clustMD(X = Y, G = G[i], CnsIndx = 8, OrdIndx = 11, Nnorms = 20000,
MaxIter = 500, model = "EVI", store.params = FALSE, scale = TRUE,
startCL = "kmeans", autoStop = TRUE, ma.band = 30, stop.tol = 0.0001)
cluster = res$cl
size = table(cluster)

r1 = c_index(data = dat4, k = G[i], lambda = lambda) 
r2 = dunn_index(data = dat4, k = G[i], lambda = lambda)
r3 = gamma_index(data = dat4, k = G[i], lambda = lambda) 
r4 = gplus_index(data = dat4, k = G[i], lambda = lambda)
r5 = mcclain_index(data = dat4, k = G[i], lambda = lambda) 
r6 = ptbiserial_index(data = dat4, k = G[i], lambda = lambda)
r7 = sil_index(data = dat4, k = G[i], lambda = lambda)
r8 = tau_index(data = dat4, k = G[i], lambda = lambda)      

val.md = c(G[i],res$BIChat, r1,r2,r3,r4,r5,r6,r7,r8)
tval.md = rbind(tval.md,val.md)  


### code for regclustMD
is.debug = TRUE
reg.res = try(regClustMD.BIC(X = Y, G = G[i], CnsIndx = 8, OrdIndx = 11, Nnorms = 20000,
MaxIter = 500, store.params = TRUE, scale = TRUE, stop.tol=0.0001, verbose.outer=TRUE, verbose.inner=FALSE))
cluster = reg.res$ind
size = table(cluster)

r1 = c_index(data = Y1, k = G[i], lambda = lambda) 
r2 = dunn_index(data = Y1, k = G[i], lambda = lambda)
r3 = gamma_index(data = Y1, k = G[i], lambda = lambda) 
r4 = gplus_index(data = Y1, k = G[i], lambda = lambda)
r5 = mcclain_index(data = Y1, k = G[i], lambda = lambda) 
r6 = ptbiserial_index(data = Y1, k = G[i], lambda = lambda)
r7 = sil_index(data = Y1, k = G[i], lambda = lambda)
r8 = tau_index(data = Y1, k = G[i], lambda = lambda)      

val.reg = c(G[i],reg.res$BIC, r1,r2,r3,r4,r5,r6,r7,r8)
tval.reg = rbind(tval.reg,val.reg)  
}   

rownames(tval.md) = rownames(tval.reg) = paste("G",G,sep="=") 
colnames(tval.md) = colnames(tval.reg) =c("G","BIC","c","dunn","gamma","gplus","mcclain","ptbiserial","silhouette","tau")


