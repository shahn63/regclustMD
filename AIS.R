#############################################################################################################
# PROGRAM 
# Code for real data example : Australian Institute of Sports data
#############################################################################################################
rm(list=ls())

require(mvtnorm);require(mclust);require(clustMD);require(QUIC);require(MomTrunc)
require(Matrix);require(Rcpp);require(RcppArmadillo);require(TruncatedNormal); require(psych)
source("solver1.54.R")
sourceCpp("cal_prob.cpp") 
source("cluster_validation_final.R")

library(DAAG)
data(ais)

dat1 = ais[,c(12,13,1:11)]
dat1$sex = ifelse(dat1$sex=="f",2,1)         
dat1$sport = ifelse(dat1$sport=="B_Ball",1,
             ifelse(dat1$sport=="Netball",1,
             ifelse(dat1$sport=="Tennis",1,
             ifelse(dat1$sport=="W_Polo",2,
             ifelse(dat1$sport=="Row",2,
             ifelse(dat1$sport=="T_Sprnt",3,
             ifelse(dat1$sport=="T_400m",3,
             ifelse(dat1$sport=="Field",3,
             ifelse(dat1$sport=="Gym",3,NA)))))))))
             
dat2 = dat1[,c(3:11,        # Order variables (Continuous,
               1,           #                  binary, 
               2)]          #                  nominal)
dat3 = data.frame(dat2)   
            
dat3[,1:9] = scale(dat3[,1:9])
dat3 = dat3[complete.cases(dat3),]
dat2 = dat2[complete.cases(dat2),]

###############################################################################
## Correlation matrix estimation
corr1 = mixedCor(dat2, c=c(1:9), p=11,d=10)
cor_res1 = round(corr1$rho,3)

###############################################################################
## Inverse matrix estimation               
inv_cov_raw = solve(cov(dat2))
round(inv_cov_raw,3) 
inv_cov_trasformed = solve(cov(dat3))
round(inv_cov_trasformed,3)

dat4 = dat3
dat4[,10] = as.factor(dat4[,10])
dat4[,11] = as.factor(dat4[,11])

lambda = lambdaest(dat4) ;lambda       

tval.md = tval.reg = c()
G = 2:5

for (i in 1:length(G)){
set.seed(123)

### code for clustMD
is.debug = TRUE
res_1 = clustMD(X = dat3, G = G[i], CnsIndx = 9, OrdIndx = 10, Nnorms = 20000,
MaxIter = 500, model = "EVI", store.params = FALSE, scale = TRUE,
startCL = "kmeans", autoStop = TRUE, ma.band=30, stop.tol=0.0001)
cluster = res_1$cl
size = table(cluster)

r1 = c_index(data = dat4, k = G[i], lambda = lambda) 
r2 = dunn_index(data = dat4, k = G[i], lambda = lambda)
r3 = gamma_index(data = dat4, k = G[i], lambda = lambda) 
r4 = gplus_index(data = dat4, k = G[i], lambda = lambda)
r5 = mcclain_index(data = dat4, k = G[i], lambda = lambda) 
r6 = ptbiserial_index(data = dat4, k = G[i], lambda = lambda)
r7 = sil_index(data = dat4, k = G[i], lambda = lambda)
r8 = tau_index(data = dat4, k = G[i], lambda = lambda)      

val.md = c(G[i],res_1$BIChat, r1,r2,r3,r4,r5,r6,r7,r8)
tval.md = rbind(tval.md,val.md)  


### code for regClustMD
reg.res = try(regClustMD(X = dat3, G = G[i], CnsIndx = 9, OrdIndx = 10, Nnorms = 20000, lambda=0,
MaxIter = 500, store.params = TRUE, scale = TRUE, stop.tol=0.0001) )
cluster = reg.res$ind
size = table(cluster)

r1 = c_index(data = dat4, k = G[i], lambda = lambda) 
r2 = dunn_index(data = dat4, k = G[i], lambda = lambda)
r3 = gamma_index(data = dat4, k = G[i], lambda = lambda) 
r4 = gplus_index(data = dat4, k = G[i], lambda = lambda)
r5 = mcclain_index(data = dat4, k = G[i], lambda = lambda) 
r6 = ptbiserial_index(data = dat4, k = G[i], lambda = lambda)
r7 = sil_index(data = dat4, k = G[i], lambda = lambda)
r8 = tau_index(data = dat4, k = G[i], lambda = lambda)      

val.reg = c(G[i],reg.res$BIC, r1,r2,r3,r4,r5,r6,r7,r8)
tval.reg = rbind(tval.reg,val.reg)  
}   

rownames(tval.md) = rownames(tval.reg) = paste("G",G,sep="=") 
colnames(tval.md) = colnames(tval.reg) =c("G","BIC","c","dunn","gamma","gplus","mcclain","ptbiserial","silhouette","tau")
