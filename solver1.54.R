cal_acc = function(avec, bvec) {
  n = length(avec)
  tbl = table(avec, bvec)
  acc = max( sum(diag(tbl))/n,  1 - sum(diag(tbl))/n )
  return(acc)
}

regClustMD = function(X, G, CnsIndx, OrdIndx, lambda = 10,
  Nnorms = 1e+05, MaxIter = 100, store.params = TRUE, scale = TRUE,
  stop.tol=0.001, verbose=TRUE,
  seed=NULL, ind.init=NULL, Omega.stack.init=NULL, Sigma.stack.init=NULL,
  grp=NULL) {
  
  ## @brief Runs a proposed procedure with a fixed tuning parameter. Usage is similar to clsutmd.
  
  ## @dependency  TruncatedNormal, tmvtnorm, glasso
  
  ## @param X, data matrix of size N x D
  ## @param G, integer, the number of groups
  ## @param CnsIndx, integer, the index of the last continuous variable
  ## @param OrdIndx, integer, the index of the last ordinal variable
  ## @param lambda, numeric, tuning parameter for the graphical lasso in the M step
  ## @param Nnorms, integer, the size for the MC samples for TruncatedNormal::mvNcdfin the E step, default 1e+5
  ## @param MaxIter, integer, maximum number of outer-loop iterations, default 500
  ## @param store.params, logical, whether or not store all the updates, default TRUE
  ## @param scale, logical, whether or not scale first, default TRUE
  ## @param stop.tol, tolerence for convergence
  
  ## @return ind, N by 1 vector of the label
  ## @return probs, N by G vector for tau_ig E(l_{ig})
  ## @return pi.vec, G by 1 vector of class probability (pi_g)
  ## @return mu, D by G array for mu_g
  ## @return Omega, D by D by G array for Omega_g,
  ## @return Sigma, D by D by G array for Sigma_g,
  ## @return loglik, final value of negative normalized loglikelihood (objFtn without penalty)
  ## @return objFtn, final value of neg normalized likelihood plus penalty
  ## @return BIC, BIC value
  ## @return objFtn.store, history of objFtn per iteration
  ## @return loglik.store, history of loglik per iteration
  ## @return stored, list of other stored value. stored[[t]] indicates the list of the stored values for the $t$-th iteration
  # return( list( ind=ind, probs=probs, pi.vec = pi.vec,
  # mu=mu, Omega=Omega, Sigma=Sigma, loglik=loglik, objFtn=objFtn, 
  # BIC=BIC, est.time=est.time, 
  # objFtn.store=objFtn.store, loglik.store=loglik.store, stored=stored ))
  
  require(TruncatedNormal)
  #require(tmvtnorm) no longer use
  #require(glasso)
  library(QUIC)
  
  # for truncated normal integral
  require(Rcpp)
  require(RcppArmadillo)
  # new package for truncated normal moments
  require(MomTrunc)
  
  is.zero.diag = FALSE 
  # organize inputs
  ###### convert to hidden variable #####
  
  if (OrdIndx < ncol(X)) {
    X1 = X[ ,1:OrdIndx]
    X2 = data.frame(X[ ,-(1:OrdIndx), drop=FALSE])
    for (j in 1:ncol(X2)) X2[ ,j] = as.factor(X2[ ,j])
    Y = cbind( X1, model.matrix( ~ . , data = X2)[ ,-1] + 1) 
    Y = as.matrix(Y)
  } else {
    Y <- as.matrix(X)  
  }
  
  N <- nrow(Y) 
  J <- ncol(Y) 
  D <- J  
  # standardize the continous part
  if (scale) {
    if (CnsIndx > 0) 
      Y[, 1:CnsIndx] <- scale(Y[, 1:CnsIndx])
  }
  # K : the number of categories 
  # K_j is NA for indices with continuous values
  K <- apply(Y, 2, max)
  if (CnsIndx > 0) { K[1:CnsIndx] <- NA }
  
  Ez <- array(NA, c(N, D, G))

  for (g in 1:G) Ez[, 1:D, g] <- Y

  perc.cutoffs2<-function (CnsIndx, Y, N) 
  {
    perc.cut <- list()
    for (j in (CnsIndx + 1):ncol(Y)) {
      perc.cut[[j]] <- qnorm(c(0, cumsum(table(Y[, j])/N)))
    }
    return(perc.cut)
  }  
  if (D > CnsIndx) {  
    perc.cut <- perc.cutoffs2(CnsIndx, Y, N)
    zlimits <- array(NA, c(N, D, 2))
    zlimits[, 1:CnsIndx, 1] <- -Inf
    zlimits[, 1:CnsIndx, 2] <- Inf
    for (j in (CnsIndx + 1):D) { 
      for (k in 1:K[j]) { 
        zlimits[Y[, j] == k, j, 1] <- perc.cut[[j]][k]
        zlimits[Y[, j] == k, j, 2] <- perc.cut[[j]][k + 1]
      }
    }
  } else { 
    perc.cut <- list()
    zlimits <- array(NA, c(N, D, 2))
  }
  
  Zinit <- matrix(NA, N, D)
  Zinit[, 1:D] <- Y[, 1:D]   
  
  if (D > CnsIndx) {
    Y1 = Y[ ,(1:CnsIndx)]
    Y2 = scale(Y[ ,(CnsIndx+1):D]) + matrix(rnorm(N*(D-CnsIndx), mean=0, sd=0.1), nrow=N)
    Ystd = cbind(Y1, Y2)
  }

  if (is.null(ind.init)) {
    res.mc = Mclust(X, G=G)
    ind = res.mc$classification
  } else {
    ind = ind.init
  }

  pi.vec <- table(ind)/N

  mu <- matrix(NA, D, G)
  for (g in 1:G) mu[, g] <- colMeans(matrix(Ystd[ind == g, ], sum(ind == g), D))

  Sigma.stack.old <- array(NA, c(D, D, G))
  Omega.stack.old <- array(NA, c(D, D, G))
  for (g in 1:G) {
    Sigma.stack.old[, , g] = cov(Ystd[ind == g, ])
    Omega.stack.old[, , g] = solve(Sigma.stack.old[, , g] + 0.001 * diag(D))
  }  
  acc=NA
  if (!is.null(grp)) { acc = cal_acc(ind, grp) }
  list_inits = list(Ystd = Ystd, acc = acc, ind = ind, pi.vec = pi.vec, mu = mu,
    Sigma = Sigma.stack.old, Omega = Omega.stack.old)
  
  
  objFtn.old = 0
  if (!is.null(Omega.stack.init) ) Omega.stack.old = Omega.stack.init
  if (!is.null(Sigma.stack.init) ) Sigma.stack.old = Sigma.stack.init
  cov.glasso = Sigma.stack.old
  
  
  acc.store = rep(NA, MaxIter)
  if (store.params == TRUE) {
    ind.store <- matrix(NA, N, MaxIter)
    tau.store <- array(NA, c(N, G, MaxIter))
    mu.store <- array(NA, c(D, G, MaxIter))
    Sigma.store <- array(NA, c(D, D, G, MaxIter))
    Omega.store <- array(NA, c(D, D, G, MaxIter))
  }
  objFtn.store = array(NA, MaxIter)
  loglik.store = array(NA, MaxIter)
  
  
  time1 = Sys.time()
  
  ############## Iteration start #################
  ############## Iteration start #################
  ############## Iteration start #################
  for (t in 1:MaxIter) {
    if (verbose) print(Sys.time())
    if (verbose) cat(sprintf("================ Iteration %d ================\n", t))
    if (verbose) cat(sprintf("================ Iteration %d ================\n", t))
    if (verbose) cat(sprintf("================ Iteration %d ================\n", t))
    # unit step
    
    
    #### E-step: 
    
    # per-group moments, first moments, second moments
    
    # < per-group moments >
    # fix g
    tau.cns.stack = array(NA, c(N, G)) 
    tau.noncns.stack = array(0, c(N, G))
    mu.noncns.stack = array(NA, c(N, D-CnsIndx, G))
    Sigma.noncns.stack = array(NA, c(D-CnsIndx, D-CnsIndx, G))
    for (g in 1:G) {
      
      ## continuous part
      if (CnsIndx > 0) {
        tau.cns <- mvtnorm::dmvnorm(matrix(Y[, 1:CnsIndx], 
          nrow = N), mean = mu[1:CnsIndx, g], sigma = Sigma.stack.old[1:CnsIndx, 
            1:CnsIndx, g], log = TRUE)
      } else {
        tau.cns <- rep(0, N)
      }
      
      ############ non-conti part ###############
      tmp.Sigma = Sigma.stack.old[ , ,g]
      tmp.mu = mu[ ,g]
      if (CnsIndx < D) {
        tmp.interaction = tmp.Sigma[(CnsIndx+1):D, 1:CnsIndx] %*% solve(tmp.Sigma[1:CnsIndx, 1:CnsIndx])
        # mu_noncns : N by (D - CnsIndx) matrix,
        mu.noncns = t(tmp.mu[(CnsIndx+1):D] + tmp.interaction %*% (t(Y[ ,1:CnsIndx]) - tmp.mu[1:CnsIndx]))
        mu.noncns.stack[ , , g] = mu.noncns 
        # Sigma_noncns : (D - CnsIndx) by (D - CnsIndx) matrix
        Sigma.noncns = tmp.Sigma[(CnsIndx+1):D, (CnsIndx+1):D] - 
          tmp.interaction %*% tmp.Sigma[1:CnsIndx, (CnsIndx+1):D]
        Sigma.noncns.stack[ , , g] = Sigma.noncns
        # testing mvNCdf
      }
      
      
      tau.cns.stack[ ,g] = tau.cns      
      if (CnsIndx < D) {
        if (verbose) cat(sprintf("Calculating integral for group %d\n", g)) 
        eigens = eigen(Sigma.noncns)
        # numerical stability
        if (min(eigens$values) < 10^-6) {
          eigens$values = eigens$values + 0.01 
          Sigma.noncns = Sigma.noncns + 0.01 * diag(ncol(Sigma.noncns))
          Sigma.noncns.stack[ , , g] = Sigma.noncns
        }
        if (D - CnsIndx > 1) {
          sqrt.Sigma.noncns = eigens$vectors %*% diag(sqrt(eigens$values)) %*% t(eigens$vectors) 
        } else {
          sqrt.Sigma.noncns = sqrt(Sigma.noncns)
        }
        prob = cal_prob( zlimits_lower = t(zlimits[ , (CnsIndx+1):D, 1]), 
          zlimits_upper = t(zlimits[ , (CnsIndx+1):D, 2]),
          mu_noncns = t(mu.noncns), 
          Sigma_noncns_sqrt = sqrt.Sigma.noncns,
          Nnorms = 1000000, dim = D - CnsIndx, sampleSize = N
        )$prob
        prob[is.na(prob) | is.nan(prob) | is.null(prob)] = 0.00001 ;
        prob[prob < 0.00001] = 0.00001 ; prob[prob > 0.99999] = 0.99999 ;
        tau.noncns = as.numeric(log(prob)) ;
        tau.noncns.stack[ ,g] = tau.noncns
      }
    }
    
    log.tau.unnormalized = log(matrix(rep(pi.vec, rep(N,G)), nrow=N)) + 
      tau.cns.stack + tau.noncns.stack
    tau.unnormalized = exp(log.tau.unnormalized)

    tau = tau.unnormalized / rowSums(tau.unnormalized)

    tau[tau < 1e-6] = 1e-6 ;     
    tau[tau > 1 - 1e-6] = 1 - 1e-6 ;         
    tau.sum = colSums(tau) ; # = N_g

    m.stack = array(NA, c(N, D - CnsIndx, G))
    V.stack = array(NA, c(N, D - CnsIndx, D - CnsIndx, G))
    
    if (CnsIndx < D) {
      for (g in 1:G) {
        for (i in 1:N) {
          if ((i %% 100 == 1) & verbose) cat(sprintf("Calculating moments for the %dth sample of %dth group\n", i, g)) 
          lower = zlimits[i, (CnsIndx+1):D, 1] - mu.noncns.stack[i, ,g]
          upper = zlimits[i, (CnsIndx+1):D, 2] - mu.noncns.stack[i, ,g]

          meanvarTMD = MomTrunc::meanvarTMD(lower = lower, upper = upper, mu=rep(0,length(lower)),
            Sigma = as.matrix(Sigma.noncns.stack[, , g]), dist="normal")
          tmp.mean = meanvarTMD$mean
          tmp.var = meanvarTMD$varcov
          tmp.zlimits = rbind(lower,upper)
          tmp.inds = apply( tmp.zlimits, 2, function(x) which.min(abs(x)) )
          if (sum(is.nan(tmp.mean)) > 0) stop ("error:: sum(is.nan(tmp.mean)) > 0")
          if (sum(is.nan(tmp.var)) > 0) stop ("error:: sum(is.nan(tmp.mean)) > 0") 
          m.stack[i, ,g] = tmp.mean
          V.stack[i, , ,g] = tmp.var + tmp.mean %*% t(tmp.mean)
        }
      }
    }
    
    if (verbose) cat(sprintf("calculating old gram matrix (to calculate the obj ftn of the E-step)...\n")) 
    gram.raw = array(NA, c(N, D, D, G))
    cov.glasso.E = array(NA, c(D, D, G))
    objFtn.stack = array(NA, G)
    loglik.stack = array(NA, G)
    for (g in 1:G) {
      for (i in 1:N) {
        gram.raw[i, 1:CnsIndx, 1:CnsIndx, g] = 
          (Y[i, 1:CnsIndx] - mu[1:CnsIndx, g]) %*%  t(Y[i, 1:CnsIndx] - mu[1:CnsIndx, g])
        if (CnsIndx < D) {
          gram.raw[i, 1:CnsIndx, (CnsIndx+1):D, g] = 
            (Y[i, 1:CnsIndx] - mu[1:CnsIndx, g]) %*%  t(m.stack[i, ,g] - mu[(CnsIndx+1):D, g])
          gram.raw[i, (CnsIndx+1):D, 1:CnsIndx, g] = t( gram.raw[i, 1:CnsIndx, (CnsIndx+1):D, g] ) 
          gram.raw[i, (CnsIndx+1):D, (CnsIndx+1):D, g] =
            V.stack[i, , ,g] - 2 * m.stack[i, ,g] %*% t(mu[(CnsIndx+1):D, g]) +
            mu[(CnsIndx+1):D, g] %*%  t(mu[(CnsIndx+1):D, g])
        }
      }
    }    
    for (g in 1:G) {
      cov.glasso.E[ , ,g] = apply( tau[ ,g] * gram.raw[ , , ,g] , c(2,3), sum ) / sum(tau[ ,g])
    }    
    
    for (g in 1:G) {
      loglik.stack[g] = log( pi.vec[g] ) + log(det(Omega.stack.old[ , ,g]))/2 - 
        sum(diag(cov.glasso.E[ , ,g] %*% Omega.stack.old[ , ,g]))/2
      objFtn.stack[g] = loglik.stack[g] - lambda  / 2 * sum(abs(  Omega.stack.old[ , ,g] )) 
    }
    objFtn.new = sum(tau.sum * objFtn.stack)
    loglik.new = sum(tau.sum * loglik.stack)
    if (verbose) cat(("After E-step\n"))
    if (verbose) print(objFtn.new)
    
    #### M-step
    Omega.stack = array(NA, c(D, D, G))
    Sigma.stack = array(NA, c(D, D, G))
    niter.stack = array(NA, G)
    
    # updated probability (pi_g)
    pi.vec.new = colSums(tau) / sum(tau)
    ind.new = apply(tau, 1, which.max)
    
    
    for (g in 1:G) {
      loglik.stack[g] = log( pi.vec.new[g] ) + log(det(Omega.stack.old[ , ,g]))/2 - 
        sum(diag(cov.glasso.E[ , ,g] %*% Omega.stack.old[ , ,g]))/2
      objFtn.stack[g] = loglik.stack[g] - lambda  / 2 * sum(abs(  Omega.stack.old[ , ,g] )) 
    }
    objFtn.new = sum(tau.sum * objFtn.stack)
    loglik.new = sum(tau.sum * loglik.stack)
    if (verbose) cat(("After pi update\n"))
    if (verbose) print(objFtn.new)
    
       
    # updated mean (mu_g) 
    mu.new = array(NA, c(D, G))
    for (g in 1:G) {
      mu.new[1:CnsIndx, g] = colSums(tau[ , g] * Y[ , 1:CnsIndx]) / sum(tau[ , g])
      if (CnsIndx < D){ mu.new[(CnsIndx+1):D, g] = colSums(as.matrix(tau[ , g] * m.stack[ , ,g])) / sum(tau[ , g]) }
    }
    
    # covariance matrix: for the input of glasso
    gram.raw = array(NA, c(N, D, D, G))
    cov.glasso = array(NA, c(D, D, G))
    
    if (verbose) cat(sprintf("calculating new gram matrix...\n")) 
    for (g in 1:G) {
      for (i in 1:N) {
        gram.raw[i, 1:CnsIndx, 1:CnsIndx, g] = 
          (Y[i, 1:CnsIndx] - mu.new[1:CnsIndx, g]) %*%  t(Y[i, 1:CnsIndx] - mu.new[1:CnsIndx, g])
        if (CnsIndx < D) {
          gram.raw[i, 1:CnsIndx, (CnsIndx+1):D, g] = 
            (Y[i, 1:CnsIndx] - mu.new[1:CnsIndx, g]) %*%  t(m.stack[i, ,g] - mu.new[(CnsIndx+1):D, g])
          gram.raw[i, (CnsIndx+1):D, 1:CnsIndx, g] = t( gram.raw[i, 1:CnsIndx, (CnsIndx+1):D, g] ) 
          gram.raw[i, (CnsIndx+1):D, (CnsIndx+1):D, g] =
            V.stack[i, , ,g] - 2 * m.stack[i, ,g] %*% t(mu.new[(CnsIndx+1):D, g]) +
            mu.new[(CnsIndx+1):D, g] %*%  t(mu.new[(CnsIndx+1):D, g])
        }
      }
    }

    if (verbose) cat(sprintf("Fitting graphical lassos...\n"))
    for (g in 1:G) {
      cov.glasso[ , ,g] = apply( tau[ ,g] * gram.raw[ , , ,g] , c(2,3), sum ) / sum(tau[ ,g])
    }
    
    
    for (g in 1:G) {
      loglik.stack[g] = log( pi.vec.new[g] ) + log(det(Omega.stack.old[ , ,g]))/2 - 
        sum(diag(cov.glasso[ , ,g] %*% Omega.stack.old[ , ,g]))/2
      objFtn.stack[g] = loglik.stack[g] - lambda / 2 * sum(abs(  Omega.stack.old[ , ,g] )) 
    }
    objFtn.new = sum(tau.sum * objFtn.stack)
    loglik.new = sum(tau.sum * loglik.stack)
    if (verbose) cat(("After mu update\n"))
    if (verbose) print(objFtn.new)
    
    for (g in 1:G) {
      quic.obj = QUIC(S = cov.glasso[ , ,g], rho = lambda, maxIter=500, msg=0,
        W.init = Sigma.stack.old[ , ,g], X.init = Omega.stack.old[ , ,g]) # iteration할때 warm으로 바꾸고 w.init이랑 wi.init 명시
      Omega.stack[ , ,g] = quic.obj$X
      if (min(eigen(Omega.stack[ , ,g])$values) < 1e-3) { Omega.stack[ , ,g] = Omega.stack[ , ,g] + diag(1e-3, D)}
      if ( sum(diag(Omega.stack[ , ,g]) < 1e-3) > 0 ) { is.zero.diag = TRUE } 
      Sigma.stack[ , ,g] = quic.obj$W
      if (min(eigen(Sigma.stack[ , ,g])$values) < 1e-3) { Sigma.stack[ , ,g] = Sigma.stack[ , ,g] + diag(1e-3, D)}
      loglik.stack[g] = log( pi.vec.new[g] ) + log(det(Omega.stack[ , ,g]))/2 - 
        sum(diag(cov.glasso[ , ,g] %*% Omega.stack[ , ,g]))/2
      objFtn.stack[g] = loglik.stack[g] - lambda  / 2 * sum(abs(  Omega.stack[ , ,g] )) 
      niter.stack[g] = quic.obj$iter
    }
    
    if (verbose) cat(sprintf("Summarizing results...\n")) 
    objFtn.new = sum(tau.sum * objFtn.stack)
    loglik.new = sum(tau.sum * loglik.stack)
    diff = (objFtn.new - objFtn.old) / objFtn.old
    
    if (verbose) cat(("After sigma update\n"))
    if (verbose) print(objFtn.new)
    
    acc = NA
    if (!is.null(grp)) {
      acc = cal_acc(ind.new, grp)
      acc.store[t] = acc
    }
    if (store.params == TRUE) {
      tau.store[ , ,t] = tau
      ind.store[ ,t] = ind.new
      mu.store[ , ,t] = mu.new
      Sigma.store[ , , ,t] = Sigma.stack
      Omega.store[ , , ,t] = Omega.stack
    }
    objFtn.store[t] = objFtn.new 
    loglik.store[t] = loglik.new
    
    
    #### convergence check
    
    # likelihood value
    if (verbose) cat(sprintf("objFtn diff: %.4f \n", abs(diff)))
    if ( abs(diff) < stop.tol ) break
    
    # parameter value
    diff2 = max( c(max(abs(mu.new - mu)), max(abs(pi.vec.new - pi.vec)), 
      max(abs(Sigma.stack - Sigma.stack.old)), max(abs(Omega.stack - Omega.stack.old)) ) )
    if (verbose) cat(sprintf("parameter max diff: %.4f \n", diff2))
    if (verbose & (!is.null(grp))) cat(sprintf("acc: %.2f \n", acc.store[t]))
    if ( abs(diff2) < stop.tol ) break
    
    if ( is.zero.diag ) break
    
    # else
    ind = ind.new
    pi.vec = pi.vec.new
    tau.old = tau
    Sigma.stack.old = Sigma.stack
    Omega.stack.old = Omega.stack
    objFtn.old = objFtn.new
    
    # unit step ends
  }
  ## end of the loop
  
  probs = tau # tau_ig E(l_{ig})
  Omega = Omega.stack # size D by D by G Omega_g
  Sigma = Sigma.stack # size D by D by G Sigma_g
  ind = ind # n by 1 vector of the label
  pi.vec = pi.vec.new
  mu = mu.new  # D by G
  loglik = loglik.new 
  df = G + G*D + G*D + (sum( abs(Omega) > 10^-4 ) - G*D)/2
  BIC = -2 * loglik + df * log( sum(tau.sum) )
  if (is.zero.diag) BIC = NA
  objFtn = objFtn.new
  if (store.params == TRUE) {
    stored = list(tau=tau.store, mu=mu.store, Sigma=Sigma.store, Omega=Omega.store)
  }
  
  time2 = Sys.time()
  est.time = time2 - time1
  print(est.time)
  return( list( ind=ind, probs=probs, pi.vec = pi.vec,
    mu=mu, Omega=Omega, Sigma=Sigma, loglik=loglik, objFtn=objFtn, 
    lambda = lambda, BIC=BIC, est.time=est.time, iter=t, acc=acc, acc.store=acc.store,
    list_inits = list_inits, objFtn.store=objFtn.store, loglik.store=loglik.store, stored=stored ))
}


regClustMD.BIC = function(X, G, CnsIndx, OrdIndx, lambdaseq = 2^seq(from=1, to=-15, by=-1),
  Nnorms = 1e+05, MaxIter = 100, store.params = TRUE, scale = TRUE,
  start.mode="cold", stop.tol=0.001, verbose.outer=TRUE, verbose.inner=FALSE,
  seed=NULL, ind.init=NULL, Omega.stack.init=NULL, Sigma.stack.init=NULL,
  grp=NULL) {
  
 
  lambda.len = length(lambdaseq)
  if (is.null(seed)) { seedmat = 1:lambda.len } else { seedmat = seq(from=seed, by=1, length=lambda.len) }
  res.obj.seq = list() 
  BIC.seq = c()
  lambdas = c()
  # for each lambda
  for (k in 1:lambda.len) {
    lambda = lambdaseq[k]
    myseed = seedmat[k]
    if (verbose.outer) cat(sprintf("Lambda %d/%d (%.6f)\n", k, lambda.len, lambda))
    # initialize
    Omega.stack.init = Sigma.stack.init = NULL
    res.obj.try = try(regClustMD(X=X, G=G, CnsIndx=CnsIndx, OrdIndx=OrdIndx, lambda=lambda, 
      Nnorms=Nnorms, MaxIter=MaxIter, store.params=store.params, scale=scale,
      stop.tol=stop.tol, seed=myseed, verbose=verbose.inner,
      Omega.stack.init=Omega.stack.init, Sigma.stack.init=Sigma.stack.init, ind.init=ind.init, grp=grp))
    if (class(res.obj.try) != "try-error") {
      res.obj.seq[[k]] = res.obj.try
      BIC.seq[k] = res.obj.try$BIC  
      if (verbose.outer) { cat(sprintf("Iter %d, Elapsed %.2fs, BIC %.4f, acc %.2f\n", 
        res.obj.try$iter, res.obj.try$est.time, res.obj.try$BIC, res.obj.try$acc)) }
      if (is.na(res.obj.try$BIC) & verbose.outer) {
        cat(sprintf("Warning: Small lambda may have caused BIC=NA due to numerical unstability.\n"))
        cat(sprintf("So terminating lambda search and returning the best lambda so far\n"))
        break }
    } else {
      res.obj.seq[[k]] = list( ind=NA, probs=NA, pi.vec =NA,
        mu=NA, Omega=NA, Sigma=NA, loglik=NA, objFtn=NA, 
        lambda = lambda, BIC=NA, est.time=NA, iter=NA, acc=NA, 
        acc.store=NA, objFtn.store=NA, loglik.store=NA, stored=NA )
      BIC.seq[k] = NA
      if (verbose.outer) cat(sprintf("Warning: error occurred, moving to next step...\n"))
    }
    
  }
  
  ind.opt = which.min(BIC.seq)
  res.obj.opt = res.obj.seq[[ind.opt]]
  
  return(c(res.obj.opt, lambdaseq=list(lambdaseq), BICseq=list(BIC.seq), allres=list(res.obj.seq) ))
}

