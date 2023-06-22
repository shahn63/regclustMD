
##############################################################################################################
create.dist_kpa <- function(lambda = lambdaest , data1, data2){
  data <- rbind(data1, data2)
  numvars <- sapply(data, is.numeric)
  catvars <- sapply(data, is.factor)
  
  d1 <- (data[1,numvars, drop = FALSE]-data[2,numvars, drop = FALSE])^2
  d1[is.na(d1)] <- 0
  if(length(lambda) == 1) d1 <- sum(d1)
  if(length(lambda) > 1) d1 <- as.matrix(d1, nrow = 1) %*% matrix(lambda[numvars], ncol = 1)
  
  d2 <- sapply(which(catvars), function(j) return(data[1,j] != data[2,j]))
  d2[is.na(d2)] <- FALSE
  if(length(lambda) == 1) d2 <- lambda * sum(d2)
  if(length(lambda) > 1) d2 <- matrix(d2, nrow = 1) %*% matrix(lambda[catvars], ncol = 1)
  
  return(d1 + d2)
}

##############################################################################################################
create.N_w <- function(size){
  sum((size*(size-1))/2)
}

##############################################################################################################
lambdaest = function(x = NULL) {
  numvars = sapply(x, is.numeric)
  anynum = any(numvars)
  catvars = sapply(x, is.factor)
  anyfact = any(catvars)
  vnum <- mean(sapply(x[,numvars, drop = FALSE], var, na.rm = TRUE))
  vcat <- mean(sapply(x[,catvars, drop = FALSE], function(z) return(1-sum((table(z)/sum(!is.na(z)))^2))))
  lambdaest = vnum/vcat
  return(lambdaest)
}

##############################################################################################################
create.S_w_kpa <- function(data,size){
  numvars <- sapply(data, is.numeric)
  catvars <- sapply(data, is.factor)
  S_w <- 0
  for(k in seq_along(size)){
    x_k <- data[which(cluster==k),]
    if(size[k]>1){

      d1 <- (x_k[rep(1:(size[k]-1),times=(size[k]-1):1, each=1), numvars, drop = FALSE] -
               x_k[unlist(lapply(2:size[k], seq, to=size[k])), numvars, drop = FALSE])^2
      d1[is.na(d1)] <- 0
      if(length(lambda) == 1) d1 <- rowSums(d1)
      if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
      
      d2 <- x_k[rep(1:(size[k]-1),times=(size[k]-1):1, each=1),
                catvars, drop = FALSE] != x_k[unlist(lapply(2:size[k], seq, to=size[k])), catvars, drop = FALSE]
      d2[is.na(d2)] <- FALSE
      if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
      if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
      
      S_w <- S_w + sum(d1 + d2)
    }
  }
  return(S_w)
}

##############################################################################################################
create.S_b_kpa <- function(data,size){
  numvars = sapply(data, is.numeric)
  catvars = sapply(data, is.factor)
  S_b <- 0
  for(k in 1:(length(size)-1)){
    x_k <- data[which(cluster==k),]
    for(l in (k+1):length(size)){
      x_l <- data[which(cluster==l),]
      
      d1 <- (x_k[rep(1:size[k], each=size[l]), numvars, drop = FALSE] -
               x_l[rep(1:size[l], times=size[k]), numvars, drop = FALSE])^2
      d1[is.na(d1)] <- 0
      if(length(lambda) == 1) d1 <- rowSums(d1)
      if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]

      d2 <- x_k[rep(1:size[k], each=size[l]),
                catvars, drop = FALSE] != x_l[rep(1:size[l], times = size[k]), catvars, drop = FALSE]
      d2[is.na(d2)] <- FALSE
      if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
      if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]

      S_b <- S_b + sum(d1 + d2)
    }
  }
  return(S_b)
}

##############################################################################################################
c_index = function(data = NULL, k = NULL, #S_sort = NULL, 
           lambda = NULL){

     n <- length(cluster)    #####
      S_all <- matrix(numeric(), ncol = n, nrow = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          S_all[i,j] <- create.dist_kpa(lambda, data1 = data[i,], data2 = data[j,])
        }
      }
      S_sort <- sort(S_all)
    
    S_w_kpa <- create.S_w_kpa(data,size)
    
    N_w <- create.N_w(size)
    S_min <- sum(head(S_sort,n = N_w))
    S_max <- sum(tail(S_sort,n = N_w))
    
    if(S_min != S_max){
      index <- (S_w_kpa - S_min)/(S_max - S_min)
      }
    return(index)  
    }

##############################################################################################################
dunn_index = function(#object = NULL, 
data = NULL, k = NULL, #kp_obj = "optimal", 
lambda = NULL){
  
    numvars <- sapply(data, is.numeric)
    cond <- all(!as.logical(numvars*lambda))
    
    k = length(size)
    
    #determine d(C_i,C_j)
    min_CiCj <- matrix(numeric(k*k), ncol = k, nrow = k)
    for(i in 1:(k-1)){
      xi <- data[which(cluster == i),]
      for(j in (i+1):k){
        xj <- data[which(cluster == j),]
        min_ij <- create.dist_kpa(lambda, data1 = xi[1,], data2 = xj[1,])
        for(l in 1:size[i]){
          for(m in 1:size[j]){
            min_neu <- create.dist_kpa(lambda, data1 = xi[l,], data2 = xj[m,])
            if(min_neu < min_ij){
              min_ij <- min_neu
            }
          }
        }
        min_CiCj[i,j] <- min_ij
      }
    }
    
    if(length(min_CiCj[min_CiCj > 0]) > 0){
      numerator <- min(min_CiCj[min_CiCj > 0])
    }else{
      return(NA)
    }

    
    #determine diam(C_k)
    max_diam <- numeric(k)
    for(p in 1:k){
      xi <- data[which(cluster == p),]
      if(size[p] > 1){
        max_ij <- create.dist_kpa(lambda, data1 = xi[1,], data2 = xi[2,])
        for(l in 1:(size[p]-1)){
          for(m in (l+1):(size[p])){
            max_neu <- create.dist_kpa(lambda, data1 = xi[l,], data2 = xi[m,])
            if(max_neu > max_ij){
              max_ij <- max_neu
            }
          }
        }
        max_diam[p] <- max_ij
      }
    }
    denominator <- max(max_diam)
    
    if(is.finite(numerator/denominator)){
      return(numerator/denominator)
    }else{
      return(NA)
    }
    
   }

##############################################################################################################

gamma_index = function(#object = NULL, 
data = NULL, k = NULL, dists = NULL, #kp_obj = "optimal", 
lambda = NULL){
  

      n <- nrow(data)
      dists <- matrix(numeric(), nrow = n, ncol = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          dists[i,j] <- create.dist_kpa(lambda, data1 = data[i,], data2 = data[j,])
        }
      }
    
    s_plus = 0; s_minus = 0; dist_within = 0; dist_between = 0
    for(k in as.numeric(names(size))){
      dist_within <- c(dist_within, na.omit(as.vector(dists[as.numeric(which(cluster == k)), as.numeric(which(cluster == k))])))
      dist_between <- c(dist_between, na.omit(as.vector(dists[as.numeric(which(cluster == k)), as.numeric(which(cluster != k))])))
    }
    dist_within <- dist_within[-1]
    for(m in 1:length(dist_within)){
      s_plus <- s_plus + sum(dist_within[m] < dist_between)
      s_minus <- s_minus + sum(dist_within[m] > dist_between)
    }
    index <- (s_plus - s_minus)/(s_plus + s_minus)
    
    return(index)
 # }
}
   
##############################################################################################################

gplus_index = function(#object = NULL, 
data = NULL, k = NULL, dists = NULL, #kp_obj = "optimal", 
lambda = NULL){

    n <- nrow(data)
    if(is.null(dists)){
      dists <- matrix(numeric(), nrow=n, ncol=n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          dists[i,j] <- create.dist_kpa(lambda, data1 = data[i,], data2 = data[j,])
        }
      }
    }
    
    s_minus = 0; dist_within = 0; dist_between = 0
    for(k in as.numeric(names(size))){
      dist_within <- c(dist_within, na.omit(as.vector(dists[as.numeric(which(cluster == k)),as.numeric(which(cluster == k))])))
      dist_between <- c(dist_between, na.omit(as.vector(dists[as.numeric(which(cluster == k)),as.numeric(which(cluster != k))])))
    }
    dist_within <- dist_within[-1]
    for(m in 1:length(dist_within)){
      s_minus <- s_minus + sum(dist_within[m] > dist_between)
    }
    N_t <- n*(n - 1)/2
    
    index <- (2 * s_minus)/(N_t * (N_t - 1))
    
    return(index)
  }

##############################################################################################################

mcclain_index <- function(#object = NULL, 
data = NULL, k = NULL, kp_obj = "optimal", lambda = NULL){

    n <- nrow(data)
    S_w_kpa <- create.S_w_kpa(data,size)
    S_b_kpa <- create.S_b_kpa(data,size)
    N_w <- create.N_w(size)
    N_t <- n*(n - 1)/2
    N_b <- N_t-N_w
    
    index <- (S_w_kpa/N_w)/(S_b_kpa/N_b)
    
    return(index)

}

##############################################################################################################

 
ptbiserial_index <- function(#object = NULL, 
data = NULL, k = NULL, #s_d = NULL, kp_obj = "optimal", 
lambda = NULL){
    
    n = nrow(data)
      S_all <- matrix(numeric(), ncol = n, nrow = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          S_all[i,j] <- create.dist_kpa(lambda = lambda, data1 = data[i,], data2 = data[j,])
        }
      }
      s_d <- sd(S_all, na.rm = TRUE)
    
    S_w_kpa = create.S_w_kpa(data,size)
    S_b_kpa = create.S_b_kpa(data,size)
    N_w = create.N_w(size)
    N_t = n*(n-1)/2
    N_b = N_t-N_w
    
    index = ((S_b_kpa/N_b) - (S_w_kpa/N_w)) * sqrt(N_w * N_b/N_t^2) / s_d
    
    return(index)
}

##############################################################################################################

sil_index = function(#object = NULL, 
data = NULL, k = NULL, #kp_obj = "optimal", 
lambda = NULL){
  
    n = nrow(data); #x = data
    k = length(table(cluster))
    numvars = sapply(data, is.numeric)
    catvars = sapply(data, is.factor)
       
    protos = data; nrows = nrow(data)
    dists = matrix(NA, nrow = nrows, ncol = nrows)
    for(i in 1:nrows){
      #distances of the numeric variables
      d1 <- (data[,numvars, drop = FALSE] - matrix(rep(as.numeric(protos[i, numvars, drop = FALSE]), nrows), nrow = nrows, byrow = TRUE))^2
      d1[is.na(d1)] <- 0
      if(length(lambda) == 1) d1 <- rowSums(d1)
      if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
      
      #distances of the categorical variances
      d2 <- sapply(which(catvars), function(j) return(data[,j] != rep(protos[i,j], nrows)))
      d2[is.na(d2)] <- FALSE
      if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
      if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
      
      dists[,i] <- d1 + d2
    }
    
    cluster_dists <- matrix(numeric(nrows*k), nrow = nrows, ncol = k)
    for(i in 1:k){
      if(!(length(which(cluster == i)) == 1)){
        cluster_dists[,i] <- rowMeans(dists[,which(cluster == i)])
      }else{
        cluster_dists[,i] <- dists[,which(cluster == i)]
      }
    }
    
    #determine ai, bi and si
    a = numeric(nrows); b = numeric(nrows); s = numeric(nrows)
    for(i in 1:nrows){
      if(is.na(cluster[i])){
        # special case: no cluster assignment cluster[i] (usually resulting of all variables NA in x[i])
        s[i] <- NA
      }else{
        a[i] <- cluster_dists[i, cluster[i]]
        b[i] <- min(cluster_dists[i, -cluster[i]])
        
        if(max(a[i], b[i], na.rm = TRUE) == 0){
          # special case: a[i]=0 and b[i]=0 => s = 0, since x[i] lies equally far away (distance = 0) from both the clusters
          s[i] <- 0
        }else{
          s[i] <- (b[i] - a[i])/max(a[i],b[i], na.rm = TRUE)
        }
      }
      }
      
    if(any(table(cluster) == 1)){
      for(i in which(cluster %in% as.integer(which(table(cluster) == 1)))){
        s[i] <- 0
      }
      cat(length(which(cluster %in% as.integer(which(table(cluster) == 1))))," cluster with only one observation\n")
    }
    
    index <- mean(s, na.rm = TRUE)
    
    return(index)
 # }
}

##############################################################################################################


tau_index = function(data = NULL, k = NULL, dists = NULL, #kp_obj = "optimal", 
lambda = NULL){
  
      n = nrow(data)
      dists = matrix(numeric(), nrow = n, ncol = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          dists[i,j] <- create.dist_kpa(lambda = lambda, data1 = data[i,], data2 = data[j,])
        }
      }
    
    s_plus = 0; s_minus = 0; dist_within = 0;dist_between = 0
    for(k in as.numeric(names(size))){
      dist_within <- c(dist_within, na.omit(as.vector(dists[as.numeric(which(cluster == k)),as.numeric(which(cluster == k))])))
      dist_between <- c(dist_between, na.omit(as.vector(dists[as.numeric(which(cluster == k)),as.numeric(which(cluster != k))])))
    }
    dist_within <- dist_within[-1]
    for(m in 1:length(dist_within)){
      s_plus <- s_plus + sum(dist_within[m] < dist_between)
      s_minus <- s_minus + sum(dist_within[m] > dist_between)
    }
    
    N_t = n * (n - 1)/2
    M = combn(1:n, m = 2)
    M = rbind(M, apply(X = M, MARGIN = 2, function(x){
      if(any(is.na(c(cluster[x[1]], cluster[x[2]])))){return(-1)}else{
        if(cluster[x[1]] == cluster[x[2]]){return(1)}else{return(-1)}
      }
    }))
    t = 0
    for(i in 1:(ncol(M)-1)){
      t = t + length(which(M[3,i]*M[3,-(1:i)] > 0))
    }
    
    index <- (s_plus - s_minus)/sqrt((N_t * (N_t - 1) * 0.5 - t) * N_t * (N_t-1) * 0.5)
    
    return(index)
}

##############################################################################################################
##############################################################################################################
