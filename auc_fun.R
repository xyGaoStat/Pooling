## this file is used to source the following functions

library(Matrix)
library(knitr)
library(ggplot2)

## get MLE functions (for different pooling designs)
getMLEinGammaDist.equalPool <- function(x, p) {
  # find MLEs of alpha and beta in Gamma dist.
  # x, observed data
  # p: group size in equal pool design; 
  # **** if p=1: random individual design ****
  # output: a.hat, b.hat: MLEs of alpha and beta
  
  ###############################################
  obj.fun <- function(par){
    n = length(x)
    beta = 1 / (n*par) * sum(x)
    derivative = n*p*log(p) - n*p*digamma(p*par) - n*p*log(beta) +p*sum(log(x))
    return(derivative)
  }
  
  # MLEs 
  # uniroot
  alpha.hat <- uniroot(obj.fun, lower = 0.0001, upper = 40, tol = 1e-9)$root
  beta.hat <- 1 / (length(x)*alpha.hat) * sum(x)
  
  return(list("alpha"=alpha.hat, "beta"=beta.hat))
}


getMLEinGammaDist.onePool <- function(N, p, X_ind, X_pool) {
  ########################################################
  # find MLEs of alpha and beta in Gamma dist (one pool design).
  # N: total sample size
  # p: group size in one pool design
  # X_ind: individual data, N-p from N
  # X_pool: averaged data from the remaining p data
  # output: a.hat, b.hat: MLEs of alpha and beta
  ########################################################
  obj.fun <- function(par) {
    beta = 1 / (N*par) * (sum(X_ind) + p*X_pool)
    deriv = -(N-p)*digamma(par) - p*digamma(p*par) - (N-p)*log(beta) - p*log(beta/p) + sum(log(X_ind)) +p*log(X_pool)
    return(deriv)
  }
  
  # MLEs
  alpha.hat <- uniroot(obj.fun, lower = 0.000001, upper =40, tol = 1e-9)$root
  beta.hat <- 1 / (N*alpha.hat) * (sum(X_ind) + p*X_pool)
  return(list("alpha"=alpha.hat, "beta"=beta.hat))
}


getGroupSize.Montone <- function(N, n, a=0.5, pattern="Linear") {
  ##############################################################
  # N: Total specimens
  # n: Total sample size of pooled data
  # a: q_k = a*q_{k-1}+b (AR) or q_k = a*k + b (Linear)
  # pattern: "Linear" or "AR"
  # output: group sizes
  ##############################################################
  
  if (pattern == "AR") {
    if (a == 1) {
      b = 2*N / (n*(n+1))
    }
    else {
      b = N*(a-1)^2 / (a^(n+1) +n - a*(n+1))
    }
    q0 = 0
    q <- rep(0, n)
    q[1] <- floor(a*q0 + b)
    #print(paste("a=",a, "b=", b))
    if (q[1]==0) {
      print("warning: invalid group size (1st group size=0)")
    }
    for (i in 2:(n-1)) {
      q[i] = floor(a* q[i-1] + b) 
    }
    q[n] = N - sum(q)
    # in case q=1
    if (q[1] == 1) {
      q <- q+1
      q[n] <- q[n] - 1 - (n-1)
    }
  }
  else {
    b = N/n - a/2 * (n+1)
    if (a+b <= 1) {
      a <- 0.2
      b = N/n - a/2 * (n+1)
      if (a+b <= 1) {
        return(getGroupSize.Montone(N, n, a=0.5, pattern="AR"))
      }
    }
    q <- rep(0, n)
    for (i in 1:(n-1)) {
      q[i] = floor(a*i + b)
    }
    q[n] <- N - sum(q)
  }
  
  return(q)
}



getMLEinGammaDist.MonotoneGsize <- function(N, p, Z){
  ##############################################################
  # find MLEs of alpha and beta in Gamma dist (monotonic groupSize design).
  # N: total sample size
  # p: group sizes
  # Z: observed pooled data
  # output: a.hat, b.hat: MLEs of alpha and beta
  ##############################################################
  
  ##############################################################
  #MLEs by uniroot
  obj.fun <- function(par) {
    beta <- 1 / (N*par) * (sum(p*Z))
    deriv <- - sum(p*digamma(par*p)) - sum(p*log(beta/p)) + sum(p*log(Z))
    return(deriv)
  }
  # Results
  alpha.hat <- uniroot(obj.fun, lower = 0.000001, upper = 40, tol = 1e-9)$root
  beta.hat <- 1 / (N*alpha.hat) * (sum(p*Z))
  return(list("alpha"=alpha.hat, "beta"=beta.hat))
}

## get AUC function
getAUC <- function(a.x, b.x, a.y, b.y){
  #####################################################
  # this function returns AUC.
  # a.x: (MLE) shape parameter for X
  # b.x: (MLE) scale parameter for X
  # a.y: (MLE) shape parameter for Y
  # b.y: (MLE) scale parameter for Y
  #####################################################
  if (a.x <=0 | a.y<=0 | b.x<=0 | b.y <=0) {
    print("warning: input params should be non-negative values")
  }
  Q = b.x / (b.x + b.y) # upper bound of intergral
  auc = pbeta(Q, shape1=a.y, shape2=a.x, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  return(auc)
}
# e.g.
# getAUC(a.x=1, b.x=5, a.y=1, b.y=3)


getParaFromAUC <- function(auc) {
  ###########################################
  # auc: AUC = Pr(X>Y)
  # output: params in Gamma dist. for X and Y
  ###########################################
  
  # fix three of four params
  a.x <- 2
  b.x <- 5
  a.y <- 2
  
  obj.fun <- function(par) {
    obj.val <- getAUC(a.x=a.x, b.x=b.x, a.y=a.y, b.y=par) - auc
    return(obj.val)
  }
  b.y <- uniroot(obj.fun, lower=0.0001, upper=100, tol=1e-9)$root
  return(list("a.x"=a.x, "b.x"=b.x, "a.y"=a.y, "b.y"=b.y))
}

getParaFromAUC.fixAlpha.X <- function(auc) {
  ###########################################
  # auc: AUC = Pr(X>Y)
  # fix alpha.X, beta.X, and alpha.Y
  # output: params in Gamma dist. for X and Y
  ###########################################
  
  # fix three of four params
  a.x <- 0.4
  a.y <- 0.4
  b.y <- 5
  
  obj.fun <- function(par) {
    obj.val <- getAUC(a.x=a.x, b.x=par, a.y=a.y, b.y=b.y) - auc
    return(obj.val)
  }
  b.x <- uniroot(obj.fun, lower=1e-6, upper=1e6, tol=1e-9)$root
  return(list("a.x"=a.x, "b.x"=b.x, "a.y"=a.y, "b.y"=b.y))
}


getParaFromAUC.fixBeta.X <- function(auc) {
  ###########################################
  # auc: AUC = Pr(X>Y)
  # fix alpha.X, beta.X, and alpha.Y
  # output: params in Gamma dist. for X and Y
  ###########################################
  
  # fix three of four params
  b.x <- 5
  a.y <- 0.4
  b.y <- 5
  
  obj.fun <- function(par) {
    obj.val <- getAUC(a.x=par, b.x=b.x, a.y=a.y, b.y=b.y) - auc
    return(obj.val)
  }
  a.x <- uniroot(obj.fun, lower=1e-6, upper=1e6, tol=1e-9)$root
  return(list("a.x"=a.x, "b.x"=b.x, "a.y"=a.y, "b.y"=b.y))
}

getParaFromAUC.fixAlpha.Y <- function(auc) {
  ###########################################
  # auc: AUC = Pr(X>Y)
  # fix alpha.X, beta.X, and alpha.Y
  # output: params in Gamma dist. for X and Y
  ###########################################
  
  # fix three of four params
  a.x <- 0.4
  b.x <- 5
  a.y <- 0.4
  
  obj.fun <- function(par) {
    obj.val <- getAUC(a.x=a.x, b.x=b.x, a.y=a.y, b.y=par) - auc
    #obj.val <- getAUC(a.x=par, b.x=b.x, a.y=a.y, b.y=b.y) - auc
    return(obj.val)
  }
  b.y <- uniroot(obj.fun, lower=1e-6, upper=1e6, tol=1e-9)$root
  return(list("a.x"=a.x, "b.x"=b.x, "a.y"=a.y, "b.y"=b.y))
}


getParaFromAUC.fixBeta.Y <- function(auc) {
  ###########################################
  # auc: AUC = Pr(X>Y)
  # fix alpha.X, beta.X, and beta.Y
  # output: params in Gamma dist. for X and Y
  ###########################################
  
  # fix three of four params
  a.x <- 0.4
  b.x <- 5
  b.y <- 5
  
  obj.fun <- function(par) {
    obj.val <- getAUC(a.x=a.x, b.x=b.x, a.y=par, b.y=b.y) - auc
    return(obj.val)
  }
  a.y <- uniroot(obj.fun, lower=1e-6, upper=1e6, tol=1e-9)$root
  return(list("a.x"=a.x, "b.x"=b.x, "a.y"=a.y, "b.y"=b.y))
}


## get Partial Derivatives of AUC w.r.t params 
getPartialDerivs <- function(a.x, b.x, a.y, b.y) {
  #####################################################
  # this function returns the partial derivatives of AUC w.r.t params.
  # a.x: (MLE) shape parameter for X
  # b.x: (MLE) scale parameter for X
  # a.y: (MLE) shape parameter for Y
  # b.y: (MLE) scale parameter for Y
  #####################################################
  Q <- b.x / (b.x + b.y)
  
  ###### 1. partial A / partial alpha.x ######
  p11.denominator <- beta(a.x, a.y)^2
  p11.numerator <- 0
  integrand1 <- function(t) {
    t^(a.y-1) * log(1-t) * (1-t)^(a.x-1)
  }
  integrand2 <- function(t) {
    t^(a.y-1) * (1-t)^(a.x-1)
  }
  res1 = integrate(integrand1, lower = 0, upper = Q)$value
  res2 = integrate(integrand2, lower = 0, upper = Q)$value
  p11.numerator = res1 * beta(a.x, a.y) - res2 * beta(a.x, a.y) * (digamma(a.x) - digamma(a.x+a.y))
  p11 <- p11.numerator / p11.denominator
  
  ###### 2. partial A / partial beta.x ######
  p12 <- 1 / beta(a.x, a.y) * (b.x / (b.x+b.y))^(a.y-1) * (b.y / (b.x+b.y))^a.x * 1 / (b.x+b.y)
  
  ###### 3. partial A / partial alpha.y ######
  p21.denominator <- beta(a.x, a.y)^2
  p21.numerator <- 0
  integrand3 <- function(t) {
    log(t) * t^(a.y-1) * (1-t)^(a.x-1)
  }
  res3 <- integrate(integrand3, lower = 0, upper = Q, subdivisions=2000)$value
  #res3 <- integrate(integrand3, lower = 0, upper = Q)$value
  p21.numerator = res3*beta(a.x, a.y) - res2 * beta(a.x, a.y) * (digamma(a.y) - digamma(a.x+a.y))
  p21 <- p21.numerator / p21.denominator
  
  ###### 4. partial A / partial beta.y ######
  p22 <- (-1) / beta(a.x, a.y) * (b.x / (b.x + b.y))^a.y * (b.y / (b.x+b.y))^(a.x-1) * 1 / (b.x + b.y)
  
  # concatenate all partial derivatives 
  derivs <- matrix(c(p11, p12, p21, p22))
  return(derivs)
  
}
# e.g.
# getPartialDerivs(a.x=1, b.x=3, a.y=1, b.y=2)


## data generating functions (for different pooling designs)
dataGen4Poolings <- function(method, N, n, a.true, b.true) {
  #################################################################################
  # method: different numbers represent pooling methods:
  #         method 1: "random sampling design"; method 2: "equal-pool design"; 
  #         method 3: "one-pool design"; 
  #         method 4: "monontonically increasing (decreasing) design".
  # N: total sample size of original data
  # n: total sample size of pooled data
  # a.true: (true) shape parameter of Gamma Dist.
  # b.true: (true) scale parameter of Gamma Dist.
  #################################################################################
  #additional_args <- list(...)
  
  #################################################################################
  # Data Generating
  all.data <- rgamma(N, shape = a.true, scale = b.true)
  
  # random sampling design
  if (method == 1) {
    #print("random sampling design")
    #n <- additional_args$n # total sample size of random sampling 
    #print(paste("Total sample size=", N, "/", "sample size of random sampling=", n))
    data <- sample(all.data, n)
    return(data)
  }
  
  # equal-pool design
  else if (method == 2) {
    #print("equal-pool design")
    #n <- additional_args$n # total sample size of equal-pool design
    p = floor(N/n) # equal group size
    #print(paste("Total sample size=", N, "/","equal group size=", p, "/","sample size of equal pool=", n))
    data <- rep(0, n)
    for (k in 1: (n-1)) {
      data[k] = mean(all.data[(p*(k-1)+1): (p*k)])
    }
    data[n] <- mean(all.data[(p*(n-1)+1):N])
    return(list("data"=data, "group_size"=p))
  }
  
  # one-pool design
  else if (method == 3) {
    #print("one-pool design")
    #p <- additional_args$p # p: group size in the pool-data
    #print(paste("Total sample size=",N,"/", "individual sample size=", N-p,"/", "pooled sample size=1"))
    p <- N - (n-1) # group size of one-pool
    data.ind = all.data[1:(N-p)] # N-p individual sample
    data.pool = mean(all.data[(N-p+1):N])
    return(list("ind"=data.ind, "pool"=data.pool))
  }
  
  else if (method == 4) {
    #print("monontonically increasing (decreasing) design")
    p <- getGroupSize.Montone(N, n) # group sizes
    #print(p)
    Z <- rep(0, n)
    r_k = 0
    for (k in 1:n) {
      if (p[k] == 1) {
        Z[k] = all.data[r_k+1] }
      else {
        Z[k] <- mean(all.data[(r_k+1):(r_k+p[k])])
      }
      r_k = r_k + p[k]
    }
    return(Z)
  }
  
  else {
    print("warning: method does not exist or have not been developed!")
  }
}


## getMLE function
getMLE.fun <- function(method, N.x, N.y, n.x, n.y,
                       a.X.true, b.X.true, a.Y.true, b.Y.true) {
  #################################################################################
  # This function is to derive MLE of params given simulated data which generated from Gamma dist. 
  
  # Common inputs for different methods:
  # method: different numbers represent pooling methods:
  #         method 1: "random sampling design"; method 2: "equal-pool design"; 
  #         method 3: "one-pool design"; 
  #         method 4: "monontonically increasing (decreasing) design".
  # N.x: total sample size for X (diseased group);
  # N.y: tital sample size for Y (Non-diseased group)
  # n.x: total sample siZe of pooled data for X
  # n.y: total sample size of pooled data for Y 
  # a.X.true: (true) shape parameter of Gamma Dist. for X (Diseased group)
  # b.X.true: (true) scale parameter of Gamma Dist. for X (Diseased group)
  # a.Y.true: (true) shape parameter of Gamma Dist. for Y (Non-diseased group)
  # b.Y.true: (true) scale parameter of Gamma Dist. for Y (Non-diseased group)
  
  #################################################################################
  #other_args <- list(...)
  
  #################################################################################
  a.X.MLE <- NULL
  b.x.MLE <- NULL
  a.Y.MLE <- NULL
  b.Y.MLE <- NULL
  
  # method 1: "random sampling design"
  if (method == 1) {
    #n.x <- other_args$n.x # total sample size for random sampling design for X
    #n.y <- other_args$n.y # total sample size for random sampling design for Y
    # MLEs for X
    X.ind = dataGen4Poolings(method=method, N=N.x, n=n.x, a.true=a.X.true, b.true=b.X.true)
    res.X.mle = getMLEinGammaDist.equalPool(X.ind, 1)
    a.X.MLE = res.X.mle$alpha
    b.X.MLE = res.X.mle$beta
    # MLEs for Y
    Y.ind = dataGen4Poolings(method=method, N=N.y, n=n.y, a.true=a.Y.true, b.true=b.Y.true)
    res.Y.mle = getMLEinGammaDist.equalPool(Y.ind, 1)
    a.Y.MLE = res.Y.mle$alpha
    b.Y.MLE = res.Y.mle$beta
  }
  
  # method 2: "equal-pool design"
  else if (method==2) {
    #n.x <- other_args$n.x # total sample size of equal-pool design
    #n.y <- other_args$n.y
    # MLEs for X
    gen.X.data = dataGen4Poolings(method=method, N=N.x, n=n.x, a.true=a.X.true, b.true=b.X.true)
    ZX = gen.X.data$data
    group.size.x = gen.X.data$group_size
    res.X.mle = getMLEinGammaDist.equalPool(ZX, group.size.x)
    a.X.MLE = res.X.mle$alpha
    b.X.MLE = res.X.mle$beta
    # MLEs for Y
    gen.Y.data = dataGen4Poolings(method=method, N=N.y, n=n.y, a.true=a.Y.true, b.true=b.Y.true)
    ZY = gen.Y.data$data
    group.size.y = gen.Y.data$group_size
    res.Y.mle = getMLEinGammaDist.equalPool(ZY, group.size.y)
    a.Y.MLE = res.Y.mle$alpha
    b.Y.MLE = res.Y.mle$beta
  }
  
  # method 3: "one-pool design"
  else if (method==3) {
    #p.x = other_args$p.x  # the number of sample combined to be a one-pool for X
    #p.y = other_args$p.y  # the number of sample combined to be a one-pool for Y
    # MLEs for X
    gen.X.data <- dataGen4Poolings(method=method, N=N.x, n=n.x, a.true=a.X.true, b.true=b.X.true)
    X.ind = gen.X.data$ind
    X.pool = gen.X.data$pool
    p.x = N.x - (n.x - 1)
    res.X.mle = getMLEinGammaDist.onePool(N=N.x, p=p.x, X_ind=X.ind, X_pool=X.pool)
    a.X.MLE = res.X.mle$alpha
    b.X.MLE = res.X.mle$beta
    # MLEs for Y
    gen.Y.data <- dataGen4Poolings(method=method, N=N.y, n=n.y, a.true=a.Y.true, b.true=b.Y.true)
    Y.ind = gen.Y.data$ind
    Y.pool = gen.Y.data$pool
    p.y = N.y - (n.y - 1)
    res.Y.mle = getMLEinGammaDist.onePool(N=N.y, p=p.y, X_ind=Y.ind, X_pool=Y.pool)
    a.Y.MLE = res.Y.mle$alpha
    b.Y.MLE = res.Y.mle$beta
  }
  
  # method 4: "monontonically increasing (decreasing) design"
  else if (method==4) {
    # MLEs for X
    #print("----X----")
    Z.X = dataGen4Poolings(method=method, N=N.x, n=n.x, a.true=a.X.true, b.true=b.X.true)
    p.x = getGroupSize.Montone(N=N.x, n=n.x)
    res.X.mle = getMLEinGammaDist.MonotoneGsize(N=N.x, p=p.x, Z=Z.X)
    a.X.MLE = res.X.mle$alpha
    b.X.MLE = res.X.mle$beta
    # MLEs for Y
    #print("----Y----")
    Z.Y = dataGen4Poolings(method=method, N=N.y, n=n.y, a.true=a.Y.true, b.true=b.Y.true)
    p.y = getGroupSize.Montone(N=N.y, n=n.y)
    res.Y.mle = getMLEinGammaDist.MonotoneGsize(N=N.y, p=p.y, Z=Z.Y)
    a.Y.MLE = res.Y.mle$alpha
    b.Y.MLE = res.Y.mle$beta
  }
  
  else {
    print("Warning: This method has not been developed!")
  }
  
  return(list("a.X.MLE"=a.X.MLE, "b.X.MLE"=b.X.MLE, "a.Y.MLE"=a.Y.MLE, "b.Y.MLE"=b.Y.MLE))
}


get.THAsymVar.MLEparams <- function(method, N, n, alpha.true, beta.true) {
  #################################################################################
  # This function is to derive the asym. variance of MLE of params in Gamma Dist. 
  # (only for theoretical results)
  # method: different numbers represent pooling methods:
  #         method 1: "random sampling design"; method 2: "equal-pool design"; 
  #         method 3: "one-pool design"; 
  #         method 4: "monontonically increasing (decreasing) design".
  # N: total sample size of original specimens;
  # n: total sample size of pooled data; 
  # alpha.true: (true) shape parameter of Gamma Dist.
  # beta.true: (true) scale parameter of Gamma Dist.
  #################################################################################
  
  # method 1: "random sampling design"
  if (method == 1) {
    # asym. var function
    asymVar.randomSampling <- function(omega, a, b) {
      common.coef <- (1/omega) * b^2 / (a*trigamma(a) - 1)
      a11 <- common.coef * a / b^2
      a12 <-  - common.coef / b
      a22 <- common.coef *trigamma(a)
      mat.asym <- matrix(c(a11, a12, a12, a22), 2, 2)
      return(mat.asym)
    }
    
    omega <- n / N
    asym.var <- asymVar.randomSampling(omega, a=alpha.true, b=beta.true)
    
  }
  
  # method 2: "equal-pool design"
  if (method==2) {
    group.size = floor(N / n)
    # asym. var function
    asymVar.equalPool <- function(p, a, b) {
      common.coef <- b^2 / (p*a*trigamma(p*a) - 1)
      a11 <- common.coef * a / b^2
      a12 <-  - common.coef / b
      a22 <- common.coef * p*trigamma(p*a)
      mat.asym <- matrix(c(a11, a12, a12, a22), 2, 2)
      return(mat.asym)
    }
    
    asym.var <- asymVar.equalPool(p=group.size, a=alpha.true, b=beta.true)
    
  }
  
  # method 3: "one-pool design"
  if (method==3) {
    p = N - (n - 1) # group size of one-pool 
    # asym. var function
    asymVar.onePool <- function(N, p, a, b) {
      w <- p/N # the ratio of the number of sample in the one-pool group to the total sample size
      common.coef <- b^2 / (a*((1-w)*trigamma(a)+w*p*trigamma(p*a))-1)        
      a11 <- common.coef * a / b^2
      a12 <-  - common.coef / b
      a22 <- common.coef * ((1-w)*trigamma(a) +w*p*trigamma(p*a))
      mat.asym <- matrix(c(a11, a12, a12, a22), 2, 2)
      return(mat.asym)
    }
    
    asym.var <- asymVar.onePool(N=N, p=p, a=alpha.true, b=beta.true)
  }
  
  # method 4: "monontonically increasing (decreasing) design"
  if (method==4) {
    p <- getGroupSize.Montone(N=N, n=n)
    
    # asym. var function
    asymVar.MonGsize <- function(N, p, a, b) {
      w <- p/N
      common.coef <- b^2 / (a*sum(w*p*trigamma(p*a)) - 1)
      a11 <- common.coef * a / b^2
      a12 <- - common.coef / b
      a22 <- common.coef * sum(w*p*trigamma(p*a))
      mat.asym <- matrix(c(a11, a12, a12, a22), 2, 2)
      return(mat.asym)
    }
    
    asym.var <- asymVar.MonGsize(N=N, p=p, a=alpha.true, b=beta.true)
    
  }
  
  return(asym.var)
  
}



## get AUC (MLE & true) and its asym. variance function
getAUCandAsymVar <- function(method, N.x, N.y, n.x, n.y,
                             a.X.true, b.X.true, a.Y.true, b.Y.true) {
  #################################################################################
  # This function is to derive MLE of AUC and its asym. variance given simulated data,
  # which generated from Gamma dist. through asymptotic property.
  
  # Common Inputs for different methods:
  # method: different numbers represent pooling methods:
  #         method 1: "random sampling design"; method 2: "equal-pool design"; 
  #         method 3: "one-pool design"; 
  #         method 4: "monontonically increasing (decreasing) design".
  # N.x: total sample size for X (diseased group);
  # N.y: tital sample size for Y (Non-diseased group)
  # n.x: total sample siZe of pooled data for X
  # n.y: total sample size of pooled data for Y 
  # a.X.true: (true) shape parameter of Gamma Dist. for X (Diseased group)
  # b.X.true: (true) scale parameter of Gamma Dist. for X (Diseased group)
  # a.Y.true: (true) shape parameter of Gamma Dist. for Y (Non-diseased group)
  # b.Y.true: (true) scale parameter of Gamma Dist. for Y (Non-diseased group)
  
  #################################################################################
  #other_args <- list(...)
  
  #################################################################################
  ## step 1: get MLEs and their asym. variances
  res.MLE <- getMLE.fun(method=method, N.x=N.x, N.y=N.y, n.x=n.x, n.y=n.y,
                        a.X.true=a.X.true, b.X.true=b.X.true, a.Y.true=a.Y.true, b.Y.true=b.Y.true)
  # MLEs for X
  a.X.MLE = res.MLE$a.X.MLE
  b.X.MLE = res.MLE$b.X.MLE
  # MLEs for Y
  a.Y.MLE = res.MLE$a.Y.MLE
  b.Y.MLE = res.MLE$b.Y.MLE
  
  
  # method 1: "random sampling design"
  if (method == 1) {
    # print("random sampling design")
    # print(paste("Total sample size (for X)=", N.x, "/", "sample size of random sampling (for X)=", n.x, "/n",
    #             "Total sample size (for Y)=", N.y, "/", "sample size of random sampling (for Y)=", n.y))
    
    # asym. var function
    asymVar.randomSampling <- function(omega, a, b) {
      common.coef <- (1/omega) * b^2 / (a*trigamma(a) - 1)
      a11 <- common.coef * a / b^2
      a12 <-  - common.coef / b
      a22 <- common.coef *trigamma(a)
      mat.asym <- matrix(c(a11, a12, a12, a22), 2, 2)
      return(mat.asym)
    }
    
    ### asym. var of MLE for X
    omega.x <- n.x / N.x
    # plug in true parameters
    asym.var.X.true <- asymVar.randomSampling(omega.x, a=a.X.true, b=b.X.true)
    # plug in MLEs of parameters
    asym.var.X.MLE <- asymVar.randomSampling(omega.x, a=a.X.MLE, b=b.X.MLE)
    
    ## asym. var of MLE for Y
    omega.y <- n.y / N.y
    # plug in true parameters
    asym.var.Y.true <- asymVar.randomSampling(omega.y, a=a.Y.true, b=b.Y.true)
    # plug in MLEs of parameters
    asym.var.Y.MLE <- asymVar.randomSampling(omega.y, a=a.Y.MLE, b=b.Y.MLE)
  }
  
  # method 2: "equal-pool design"
  if (method==2) {
    # print("equal-pool design")
    group.size.x = floor(N.x / n.x)
    group.size.y = floor(N.y / n.y)
    # print(paste("Total sample size (for X)=", N.x, "/","equal group size (for X)=", group.size.x, "/","sample size of equal pool (for X)=", n.x, "/n", "Total sample size (for Y)=", N.y, "/","equal group size (for Y)=",group.size.y, "/","sample size of equal pool (for Y)=", n.y))
    
    
    # asym. var function
    asymVar.equalPool <- function(p, a, b) {
      common.coef <- b^2 / (p*a*trigamma(p*a) - 1)
      a11 <- common.coef * a / b^2
      a12 <-  - common.coef / b
      a22 <- common.coef * p*trigamma(p*a)
      mat.asym <- matrix(c(a11, a12, a12, a22), 2, 2)
      return(mat.asym)
    }
    ### asym. var of MLE for X
    # plug in true parameters
    asym.var.X.true <- asymVar.equalPool(p=group.size.x, a=a.X.true, b=b.X.true)
    # plug in MLEs of parameters
    asym.var.X.MLE <- asymVar.equalPool(p=group.size.x, a=a.X.MLE, b=b.X.MLE)
    
    ## asym. var of MLE for Y
    # plug in true parameters
    asym.var.Y.true <- asymVar.equalPool(p=group.size.y, a=a.Y.true, b=b.Y.true)
    # plug in MLEs of parameters
    asym.var.Y.MLE <- asymVar.equalPool(p=group.size.y, a=a.Y.MLE, b=b.Y.MLE)
  }
  
  # method 3: "one-pool design"
  if (method==3) {
    #print("one-pool design")
    p.x = N.x - (n.x - 1) # group size of one-pool for X
    p.y = N.y - (n.y - 1) # group size of one-pool for Y
    #print(paste("Total sample size (for X)=",N.x,"/", "individual sample size (for X)=", N.x-p.x,"/", "pooled sample size=1", "/n", "Total sample size (for Y)=",N.y,"/", "individual sample size (for Y)=", N.y-p.y,"/", "pooled sample size=1"))
    
    # asym. var function
    asymVar.onePool <- function(N, p, a, b) {
      w <- p/N # the ratio of the number of sample in the one-pool group to the total sample size
      common.coef <- b^2 / (a*((1-w)*trigamma(a)+w*p*trigamma(p*a))-1)        
      a11 <- common.coef * a / b^2
      a12 <-  - common.coef / b
      a22 <- common.coef * ((1-w)*trigamma(a) +w*p*trigamma(p*a))
      mat.asym <- matrix(c(a11, a12, a12, a22), 2, 2)
      return(mat.asym)
    }
    ### asym. var of MLE for X
    # plug in true parameters
    asym.var.X.true <- asymVar.onePool(N=N.x, p=p.x, a=a.X.true, b=b.X.true)
    # plug in MLEs of parameters
    asym.var.X.MLE <- asymVar.onePool(N=N.x, p=p.x, a=a.X.MLE, b=b.X.MLE)
    ## asym. var of MLE for Y
    # plug in true parameters
    asym.var.Y.true <- asymVar.onePool(N=N.y, p=p.y, a=a.Y.true, b=b.Y.true)
    asym.var.Y.MLE <- asymVar.onePool(N=N.y, p=p.y, a=a.Y.MLE, b=b.Y.MLE)
  }
  
  # method 4: "monontonically increasing (decreasing) design"
  if (method==4) {
    #print("Monotone GroupSize Design")
    p.x <- getGroupSize.Montone(N=N.x, n=n.x)
    p.y <- getGroupSize.Montone(N=N.y, n=n.y)
    
    # asym. var function
    asymVar.MonGsize <- function(N, p, a, b) {
      w <- p/N
      common.coef <- b^2 / (a*sum(w*p*trigamma(p*a)) - 1)
      a11 <- common.coef * a / b^2
      a12 <- - common.coef / b
      a22 <- common.coef * sum(w*p*trigamma(p*a))
      mat.asym <- matrix(c(a11, a12, a12, a22), 2, 2)
      return(mat.asym)
    }
    ### asym. var of MLE for X
    # plug in true parameters
    asym.var.X.true <- asymVar.MonGsize(N=N.x, p=p.x, a=a.X.true, b=b.X.true)
    # plug in MLEs of parameters
    asym.var.X.MLE <- asymVar.MonGsize(N=N.x, p=p.x, a=a.X.MLE, b=b.X.MLE)
    ## asym. var of MLE for Y
    # plug in true parameters
    asym.var.Y.true <- asymVar.MonGsize(N=N.y, p=p.y, a=a.Y.true, b=b.Y.true)
    asym.var.Y.MLE <- asymVar.MonGsize(N=N.y, p=p.y, a=a.Y.MLE, b=b.Y.MLE)
    
  }
  
  ## Construct joint asym. variance 
  # total sample size for X and Y 
  N.total = N.x + N.y
  prop.p = N.x / N.total # consider the joint dist. for different total sample sizes
  # plug in true params
  joint.asym.var.true = bdiag(asym.var.X.true/prop.p, asym.var.Y.true/(1-prop.p))
  # plug in MLEs
  joint.asym.var.MLE = bdiag(asym.var.X.MLE/prop.p, asym.var.Y.MLE/(1-prop.p))
  #################################################################################
  
  
  #################################################################################
  ## step 2: find MLE of AUC
  # plug in true params
  AUC.true = getAUC(a.x=a.X.true, b.x=b.X.true, a.y=a.Y.true, b.y=b.Y.true)
  # plug in MLEs
  AUC.MLE = getAUC(a.x=a.X.MLE, b.x=b.X.MLE, a.y=a.Y.MLE, b.y=b.Y.MLE)
  
  
  #################################################################################
  # step 3: find asym. var of MLE of AUC
  # in this step, all we need to do is to find the partial derivatives of AUC w.r.t. params.
  
  # Partial derivatives
  # plug in true params
  p.derivs.true = getPartialDerivs(a.x=a.X.true, b.x=b.X.true, a.y=a.Y.true, b.y=b.Y.true)
  # plug in MLEs
  p.derivs.MLE = getPartialDerivs(a.x=a.X.MLE, b.x=b.X.MLE, a.y=a.Y.MLE, b.y=b.Y.MLE)
  
  # asym. var of AUC
  # plug in true params
  asym.var.AUC.true = t(p.derivs.true) %*% joint.asym.var.true %*% p.derivs.true
  # plug in MLES
  asym.var.AUC.MLE = t(p.derivs.MLE) %*% joint.asym.var.MLE %*% p.derivs.MLE
  #################################################################################
  
  
  return(list("a.X.MLE"=a.X.MLE, "b.X.MLE"=b.X.MLE, "a.Y.MLE"=a.Y.MLE, "b.Y.MLE"=b.Y.MLE,
              "True AUC"=AUC.true, "MLE of AUC"=AUC.MLE,
              "true asym. var for X"=asym.var.X.true, "true asym. var for Y"=asym.var.Y.true,
              "MLE asym. var for X"=asym.var.X.MLE, "MLE asym. var for Y"=asym.var.Y.MLE,
              "true joint asym. var"=joint.asym.var.true, "MLE joint asym. var"=joint.asym.var.MLE,
              "true asym var of AUC"=asym.var.AUC.true, "MLE asym var of AUC"=asym.var.AUC.MLE))
  
}


## MC study

#### did not use sapply function
#getAsymVarofAUC.MC <- function(method, MC, N.x, N.y, n.x, n.y,
#                                a.X.true, b.X.true, a.Y.true, b.Y.true) {
#   #################################################################################
#   # This function is to derive the asym. variance of AUC given simulated data,
#   # which generated from Gamma dist. through Monte Carlo method.
#   
#   # Common Inputs for different methods:
#   # method: different numbers represent pooling methods:
#   #         method 1: "random sampling design"; method 2: "equal-pool design"; 
#   #         method 3: "one-pool design"; method 4: "monontonically increasing (decreasing) design";
#   #         method 5: "non-monotonical design".
#   # MC: the number of repeat times;
#   # N.x: total sample size for X (diseased group);
#   # N.y: tital sample size for Y (Non-diseased group)
#   # a.X.true: (true) shape parameter of Gamma Dist. for X (Diseased group)
#   # b.X.true: (true) scale parameter of Gamma Dist. for X (Diseased group)
#   # a.Y.true: (true) shape parameter of Gamma Dist. for Y (Non-diseased group)
#   # b.Y.true: (true) scale parameter of Gamma Dist. for Y (Non-diseased group)
#   # ...: other arguments. Different pooling designs may have additional aruguments
#   
#   #################################################################################
#   #other_args <- list(...)
#   # true AUC value
#   AUC.true = getAUC(a.x=a.X.true, b.x=b.X.true, a.y=a.Y.true, b.y=b.Y.true)
#   A = rep(0, MC)
#   
#   for (i in 1:MC) {
#     # step 1: get MLEs
#     res.MLEs <- getMLE.fun(method, N.x, N.y, n.x, n.y, a.X.true, b.X.true, a.Y.true, b.Y.true)
#     a.X.MLE = res.MLEs$a.X.MLE
#     b.X.MLE = res.MLEs$b.X.MLE
#     a.Y.MLE = res.MLEs$a.Y.MLE
#     b.Y.MLE = res.MLEs$b.Y.MLE
#     
#     # step 2: get MLE of AUC
#     AUC.hat <- getAUC(a.x=a.X.MLE, b.x=b.X.MLE, a.y=a.Y.MLE, b.y=b.Y.MLE)
#     
#     # step 3: record the diff
#     A[i] = AUC.hat - AUC.true
#     
#   }
#   
#   N.total = N.x + N.y
#   res = N.total * mean(A^2)
#   return("Var of AUC (MC)"=res)
#   
# }




### use sapply function
getAsymVarofAUC.MC <- function(method, MC, N.x, N.y, n.x, n.y,
                               a.X.true, b.X.true, a.Y.true, b.Y.true) {
  #################################################################################
  # This function is to derive the asym. variance of AUC given simulated data,
  # which generated from Gamma dist. through Monte Carlo method.
  
  # Common Inputs for different methods:
  # method: different numbers represent pooling methods:
  #         method 1: "random sampling design"; method 2: "equal-pool design"; 
  #         method 3: "one-pool design"; method 4: "monontonically increasing (decreasing) design";
  #         method 5: "non-monotonical design".
  # MC: the number of repeat times;
  # N.x: total sample size for X (diseased group);
  # N.y: tital sample size for Y (Non-diseased group)
  # a.X.true: (true) shape parameter of Gamma Dist. for X (Diseased group)
  # b.X.true: (true) scale parameter of Gamma Dist. for X (Diseased group)
  # a.Y.true: (true) shape parameter of Gamma Dist. for Y (Non-diseased group)
  # b.Y.true: (true) scale parameter of Gamma Dist. for Y (Non-diseased group)
  # ...: other arguments. Different pooling designs may have additional aruguments
  
  #################################################################################
  #other_args <- list(...)
  # true AUC value
  AUC.true = getAUC(a.x=a.X.true, b.x=b.X.true, a.y=a.Y.true, b.y=b.Y.true)
  
  
  singleRun <- function(hh) {
    set.seed(hh)
    # step 1: get MLEs
    res.MLEs <- getMLE.fun(method=method, N.x=N.x, N.y=N.y, n.x=n.x, n.y=n.y, 
                           a.X.true=a.X.true, b.X.true=b.X.true, 
                           a.Y.true=a.Y.true, b.Y.true=b.Y.true)
    a.X.MLE = res.MLEs$a.X.MLE
    b.X.MLE = res.MLEs$b.X.MLE
    a.Y.MLE = res.MLEs$a.Y.MLE
    b.Y.MLE = res.MLEs$b.Y.MLE
    # step 2: get MLE of AUC
    AUC.hat <- getAUC(a.x=a.X.MLE, b.x=b.X.MLE, a.y=a.Y.MLE, b.y=b.Y.MLE)
    return(AUC.hat - AUC.true)
  }
  
  A = sapply(1:MC, singleRun)
  
  N.total = N.x + N.y
  res = N.total * mean(A^2)
  return("Var of AUC (MC)"=res)
}


#------------------------Bootstrap Functionds----------------------------------#

pooled_data.Gen4Bootstrap <- function(method, seed_id, data, N, n) {
  set.seed(seed_id)
  resampled.data <- sample(data, replace = TRUE)
  
  # method 1: rand sampling
  if (method == 1) {
    pooled.data <- sample(resampled.data, n)
    return(pooled.data)
  }
  
  # method 2: Equal-pool
  else if (method == 2) {
    p <- floor(N/n)
    pooled.data <- rep(0, n)
    for (k in 1: (n-1)) {
      pooled.data[k] = mean(resampled.data[(p*(k-1)+1): (p*k)])
    }
    pooled.data[n] <- mean(resampled.data[(p*(n-1)+1):N])
    return(list("data"=pooled.data, "group_size"=p))
  }
  
  # method 3: One-pool
  else if (method == 3) {
    p <- N - (n-1) # group size of one-pool
    data.ind = resampled.data[1:(N-p)] # N-p individual sample
    data.pool = mean(resampled.data[(N-p+1):N])
    return(list("ind"=data.ind, "pool"=data.pool))
  }
  
  # method 4: monotone G size
  else if (method == 4) {
    p <- getGroupSize.Montone(N, n)
    Z <- rep(0, n)
    r_k = 0
    for (k in 1:n) {
      if (p[k] == 1) {
        Z[k] = resampled.data[r_k+1] }
      else {
        Z[k] <- mean(resampled.data[(r_k+1):(r_k+p[k])])
      }
      r_k = r_k + p[k]
    }
    return(Z)
  }
}


getMLE_Bootstrap <- function(method, seed_id, data.X, data.Y, N.x, N.y, n.x, n.y) {
  
  a.X.MLE <- NULL; b.x.MLE <- NULL; a.Y.MLE <- NULL; b.Y.MLE <- NULL
  
  # method 1: "random sampling design"
  if (method == 1) {
    # MLEs for X
    X.ind = pooled_data.Gen4Bootstrap(method=method, seed_id=seed_id, data=data.X, N=N.x, n=n.x)
    res.X.mle = getMLEinGammaDist.equalPool(X.ind, 1)
    a.X.MLE = res.X.mle$alpha
    b.X.MLE = res.X.mle$beta
    # MLEs for Y
    Y.ind = pooled_data.Gen4Bootstrap(method=method, seed_id=seed_id, data=data.Y, N=N.y, n=n.y)
    res.Y.mle = getMLEinGammaDist.equalPool(Y.ind, 1)
    a.Y.MLE = res.Y.mle$alpha
    b.Y.MLE = res.Y.mle$beta
  }
  
  # method 2: "equal-pool design"
  else if (method==2) {
    #n.x <- other_args$n.x # total sample size of equal-pool design
    #n.y <- other_args$n.y
    # MLEs for X
    gen.X.data = pooled_data.Gen4Bootstrap(method=method, seed_id=seed_id, data=data.X, N=N.x, n=n.x)
    ZX = gen.X.data$data
    group.size.x = gen.X.data$group_size
    res.X.mle = getMLEinGammaDist.equalPool(ZX, group.size.x)
    a.X.MLE = res.X.mle$alpha
    b.X.MLE = res.X.mle$beta
    # MLEs for Y
    gen.Y.data = pooled_data.Gen4Bootstrap(method=method, seed_id=seed_id, data=data.Y, N=N.y, n=n.y)
    ZY = gen.Y.data$data
    group.size.y = gen.Y.data$group_size
    res.Y.mle = getMLEinGammaDist.equalPool(ZY, group.size.y)
    a.Y.MLE = res.Y.mle$alpha
    b.Y.MLE = res.Y.mle$beta
  }
  
  # method 3: "one-pool design"
  else if (method==3) {
    #p.x = other_args$p.x  # the number of sample combined to be a one-pool for X
    #p.y = other_args$p.y  # the number of sample combined to be a one-pool for Y
    # MLEs for X
    gen.X.data <- pooled_data.Gen4Bootstrap(method=method, seed_id=seed_id, data=data.X, N=N.x, n=n.x)
    X.ind = gen.X.data$ind
    X.pool = gen.X.data$pool
    p.x = N.x - (n.x - 1)
    res.X.mle = getMLEinGammaDist.onePool(N=N.x, p=p.x, X_ind=X.ind, X_pool=X.pool)
    a.X.MLE = res.X.mle$alpha
    b.X.MLE = res.X.mle$beta
    # MLEs for Y
    gen.Y.data <- pooled_data.Gen4Bootstrap(method=method, seed_id=seed_id, data=data.Y, N=N.y, n=n.y)
    Y.ind = gen.Y.data$ind
    Y.pool = gen.Y.data$pool
    p.y = N.y - (n.y - 1)
    res.Y.mle = getMLEinGammaDist.onePool(N=N.y, p=p.y, X_ind=Y.ind, X_pool=Y.pool)
    a.Y.MLE = res.Y.mle$alpha
    b.Y.MLE = res.Y.mle$beta
  }
  
  else if (method==4) {
    # MLEs for X
    Z.X = pooled_data.Gen4Bootstrap(method=method, seed_id=seed_id, data=data.X, N=N.x, n=n.x)
    p.x = getGroupSize.Montone(N=N.x, n=n.x)
    res.X.mle = getMLEinGammaDist.MonotoneGsize(N=N.x, p=p.x, Z=Z.X)
    a.X.MLE = res.X.mle$alpha
    b.X.MLE = res.X.mle$beta
    # MLEs for Y
    Z.Y = pooled_data.Gen4Bootstrap(method=method, seed_id=seed_id, data=data.Y, N=N.y, n=n.y)
    p.y = getGroupSize.Montone(N=N.y, n=n.y)
    res.Y.mle = getMLEinGammaDist.MonotoneGsize(N=N.y, p=p.y, Z=Z.Y)
    a.Y.MLE = res.Y.mle$alpha
    b.Y.MLE = res.Y.mle$beta
  }
  
  else {
    print("Warning: This method has not been developed!")
  }
  
  return(list("a.X.MLE"=a.X.MLE, "b.X.MLE"=b.X.MLE, "a.Y.MLE"=a.Y.MLE, "b.Y.MLE"=b.Y.MLE))
  
}




