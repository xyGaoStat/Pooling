# This script will be invoked by Bi-Gamma research
library(rootSolve)
library(Matrix)
library(numDeriv)
#library(VGAM)


########################################################################################################################
# data generating functions
get_RV_biGamma <- function(n.size, shape1, shape2, scalePar) {
  #############################################################
  # shape1: a, shape2: b, scalePar: mu
  # returns random vectors which follows McKay's bi-Gamma
  #############################################################
  mdata <- data.frame(w1 = rgamma(nn <- n.size, shape=shape1, scale=scalePar))
  mdata <- transform(mdata, zedd = rgamma(nn, shape=shape2, scale=scalePar))
  mdata <- transform(mdata, w2 = w1 + zedd) # Z is defined as W2- w1 | W1=w
  #print(paste0("simulation ", dim(mdata)[1], " bi-Gamma samples"))
  return(mdata[, c("w1", "w2")])
}

data.generating.fullDat <- function(N, a, b, mu) {
  ##########################################################
  # This is a data generating function for Full data
  ##########################################################
  dat <- get_RV_biGamma(n.size=N, shape1=a, shape2=b, scalePar=mu)
  return(dat)
}


data.generating.randSamp <- function(N, n, a, b, mu) {
  ##########################################################
  # This is a data generating function for random sampling
  ##########################################################
  
  # Raw data
  raw.dat <- get_RV_biGamma(n.size=N, shape1=a, shape2=b, scalePar=mu)
  id.rand <- sample(1:N, n)
  dat <- raw.dat[id.rand, ]
  return(dat)
}

data.generating.equal_Pool <- function(N, n, a, b, mu) {
  ##########################################################
  # This is a data generating function for equal pooling design
  ##########################################################
  
  # Raw data
  raw.dat <- get_RV_biGamma(n.size=N, shape1=a, shape2=b, scalePar=mu)
  p = floor(N/n) # group size in equal-pooling design
  dat <- data.frame(w1=rep(NA, n), w2=rep(NA, n))
  # fill in obs
  for (k in 1:(n-1)) {
    dat[k, ] <- colMeans(raw.dat[(p*(k-1)+1):(p*k), ]) 
  }
  dat[n, ] <- colMeans(raw.dat[(p*(n-1)+1):N, ]) 
  return(dat)
}

data.generating.one_Pool <- function(N, n, a, b, mu) {
  ##########################################################
  # This is a data generating function for one pool design
  ##########################################################
  
  # Raw data
  raw.dat <- get_RV_biGamma(n.size=N, shape1=a, shape2=b, scalePar=mu)
  p <- N - n + 1  # group size in one-pool design
  dat.ind <- raw.dat[1:(N-p), ]
  dat.pool <- colMeans(raw.dat[(N-p+1):N, ])
  return(list("ind"=dat.ind, "pool"=dat.pool))
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

data.generating.MonotoneGSize <- function(N, n, a, b, mu) {
  ##########################################################
  # This is a data generating function for one pool design
  ##########################################################
  
  # Raw data
  raw.dat <- get_RV_biGamma(n.size=N, shape1=a, shape2=b, scalePar=mu)
  p <- getGroupSize.Montone(N=N,n=n) # group size in Monotone Group Size Design
  dat <- data.frame(w1=rep(NA, n), w2=rep(NA, n))
  r_k = 0
  for (k in 1:n) {
    if (p[k] == 1) {
      dat[k, ] <- raw.dat[r_k+1, ]
    }
    else {
      dat[k, ] <- colMeans(raw.dat[(r_k+1):(r_k+p[k]), ])
    }
    r_k = r_k + p[k]
  }
  return(dat)
}

data.generating.4Poolings <- function(method, N, n, a, b, mu) {
  #################################################################################
  # method: different numbers represent pooling methods:
  #         method 0: "Full data"
  #         method 1: "random sampling design"; method 2: "equal-pool design"; 
  #         method 3: "one-pool design"; method 4: "monontonically increasing design".
  # N: total sample size of original data
  # n: total sample size of pooled data
  #################################################################################
  # Full data
  if (method == 0) { data.info <- data.generating.fullDat(N=N, a=a, b=b, mu=mu)}
  # random sampling
  else if (method == 1) { data.info <- data.generating.randSamp(N = N, n = n, a=a, b=b, mu=mu) }
  # equal-pool
  else if (method == 2) {data.info <- data.generating.equal_Pool(N = N, n = n, a=a, b=b, mu=mu)  }
  # one-pool
  else if (method == 3) {data.info <- data.generating.one_Pool(N = N, n = n, a=a, b=b, mu=mu)  }
  # monontone
  else if (method == 4) {data.info <- data.generating.MonotoneGSize(N = N, n = n, a=a, b=b, mu=mu)}
  return(data.info)
}

# e.g. d <- data.generating.4Poolings(method=4, N=100, n=20, a=1,b=2,mu=1.5 )

########################################################################################################################
# find MLEs based on pooling samples
get_MLE_biGamma.AllDat.vglm <- function(dat) {
  w1 <- unlist(dat[, 1]); w2 <- unlist(dat[, 2])
  fit <- vglm(cbind(w1, w2) ~ 1, bigamma.mckay)
  res <- Coef(fit)
  return(list("a.hat"=as.numeric(res[2]), "b.hat"=as.numeric(res[3]), "mu.hat"=as.numeric(res[1])))
}

get_MLE_biGamma.AllDat <- function(data, start.point=c(1,1)) {
  # get MLE of biGamma based on ALL data / random samples using 
  obj.fun <- function(par) {
    w1 <- data[,1]; w2 <- data[, 2]
    n <- length(w1)
    a <- par[1]; b <- par[2]
    mu <- 1 / (n*(a+b)) * sum(w2)
    deriv1 <- -n*log(mu) - n * digamma(a) + sum(log(w1));
    deriv2 <- -n*log(mu) - n * digamma(b) + sum(log(w2-w1));
    return(c(deriv1, deriv2))
  }
  res.multiroot <- multiroot(obj.fun, start = start.point)
  w2 <- data[, 2]; n <- length(w2)
  a.hat <- res.multiroot$root[1]
  b.hat <- res.multiroot$root[2]
  mu.hat <- 1 / (n*(a.hat+b.hat)) * sum(w2)
  return(list("a.hat"=a.hat, "b.hat"=b.hat, "mu.hat"=mu.hat))
}
# e.g.
# MLE.res <- get_MLE_biGamma.AllDat(data = dat, start.point = c(a,b)) # true value for a and b

get_MLE_biGamma.equalPool <- function(N, n, data, start.point=c(1,1)) {
  # get MLE of biGamma based on equal_Pool sample
  p = floor(N/n)
  obj.fun <- function(par) {
    w1 <- data[,1]; w2 <- data[, 2]
    a <- par[1]; b <- par[2]
    mu <- 1 / (n*(a+b)) * sum(w2)
    #print(paste0("a=",a, ", b=", b, ", mu=", mu))
    deriv1 <- n*p*log(p) - n*p*log(mu) - n*p * digamma(p*a) + p*sum(log(w1));
    deriv2 <- n*p*log(p) - n*p*log(mu) - n*p * digamma(p*b) + p*sum(log(w2-w1));
    return(c(deriv1, deriv2))
  }
  res.multiroot <- multiroot(obj.fun, start = start.point)
  w2 <- data[, 2]; n <- length(w2)
  a.hat <- res.multiroot$root[1]
  b.hat <- res.multiroot$root[2]
  mu.hat <- 1 / (n*(a.hat+b.hat)) * sum(w2)
  return(list("a.hat"=a.hat, "b.hat"=b.hat, "mu.hat"=mu.hat))
}
# e.g.
# MLE.res <- get_MLE_biGamma.equalPool(N = N, n = n, data = dat,start.point = c(a,b)) # true value for a and b

# set.seed(279)
# data.test <- data.generating.4Poolings(method=2, N=N.y, n=n.y, a=a.y, b=b.y, mu=mu.y)
# start.point=c(a.y+rnorm(1,0,0.1),b.y+rnorm(1,0,0.1))
# get_MLE_biGamma.equalPool(N=N.y, n=n.y, data=data.test, start.point=start.point)
#   
  

get_MLE_biGamma.onePool <- function(N, n, data, start.point=c(1,1)) {
  # get MLE of biGamma based on equal_Pool sample
  p = N - n + 1
  data.ind <- data$ind; data.pool <- data$pool
  w1.ind <- data.ind[,1]; w2.ind <- data.ind[, 2]
  w1.pool <- data.pool[1]; w2.pool <- data.pool[2]
  
  obj.fun <- function(par) {
    a <- par[1]; b <- par[2]
    mu <- 1 / (N*(a+b)) * (sum(w2.ind) + p*w2.pool)
    deriv1 <- p*log(p) - p*log(mu) - p*digamma(p*a) + p*log(w1.pool) - 
      (N-p)*log(mu) - (N-p)*digamma(a) + sum(log(w1.ind))   
    deriv2 <- p*log(p) - p*log(mu) - p*digamma(p*b) + p*log(w2.pool-w1.pool) - 
      (N-p)*log(mu) - (N-p)*digamma(b) + sum(log(w2.ind - w1.ind))
    return(c(deriv1, deriv2))
  }
  res.multiroot <- multiroot(obj.fun, start = start.point)
  a.hat <- res.multiroot$root[1]
  b.hat <- res.multiroot$root[2]
  mu.hat <- 1 / (N*(a.hat+b.hat)) * (sum(w2.ind) + p*w2.pool)
  mu.hat <- as.numeric(mu.hat)
  return(list("a.hat"=a.hat, "b.hat"=b.hat, "mu.hat"=mu.hat))
}
# e.g.
# MLE.res <- get_MLE_biGamma.onePool(N = N, n = n, data = dat, start.point = c(a,b))

get_MLE_biGamma.MonotoneGSize <- function(N, n, data, start.point=c(1,1)) {
  # get MLE of biGamma based on Monotone Group size
  p <- getGroupSize.Montone(N=N,n=n)
  w1 <- data[,1]; w2 <- data[, 2]
  obj.fun <- function(par) {
    a <- par[1]; b <- par[2]
    mu <- 1 / (N*(a+b)) * sum(p*w2)
    deriv1 <- sum(p*log(p)) - N*log(mu) - sum(p*digamma(p*a)) + sum(p*log(w1))
    deriv2 <- sum(p*log(p)) - N*log(mu) - sum(p*digamma(p*b)) + sum(p*log(w2 - w1))
    return(c(deriv1, deriv2))
  }
  res.multiroot <- multiroot(obj.fun, start = start.point)
  a.hat <- res.multiroot$root[1]
  b.hat <- res.multiroot$root[2]
  mu.hat <- 1 / (N*(a.hat+b.hat)) * sum(p*w2)
  return(list("a.hat"=a.hat, "b.hat"=b.hat, "mu.hat"=mu.hat))
}
# e.g.
# MLE.res <- get_MLE_biGamma.MonotoneGSize(N = N, n = n, data = dat, start.point = c(a,b))

get_MLE_biGamma.4Poolings <- function(method, N, n, a, b, mu, start.point=c(1,1)) {
  ##################################################################################
  # method: different numbers represent pooling methods:
  #         method 0: "Full data"
  #         method 1: "random sampling design"; method 2: "equal-pool design"; 
  #         method 3: "one-pool design"; method 4: "monontonically increasing design".
  # a, b, mu: true values of Parameters
  ###################################################################################
  
  # step 1: generate pooling samples based on true parameters
  data <- data.generating.4Poolings(method=method, N=N, n=n, a=a, b=b, mu=mu)
  
  if (method == 0) {res.MLEs <- get_MLE_biGamma.AllDat(data = data, start.point = start.point) }
  #if (method == 0) {res.MLEs <-   get_MLE_biGamma.AllDat.vglm(data) }
  # random sampling shares the same getMLE function
  else if (method == 1) {res.MLEs <- get_MLE_biGamma.AllDat(data = data, start.point = start.point) } 
  #else if (method == 1) {res.MLEs <- get_MLE_biGamma.AllDat.vglm(data) } 
  else if (method == 2) {res.MLEs <- get_MLE_biGamma.equalPool(N = N, n = n, data = data, start.point = start.point)}
  else if (method == 3) {res.MLEs <- get_MLE_biGamma.onePool(N = N, n = n, data = data, start.point = start.point)}
  else if (method == 4) {res.MLEs <- get_MLE_biGamma.MonotoneGSize(N = N, n = n, data = data, start.point = start.point)}
  return(res.MLEs)
}
# e.g.
# method=3
# N=1000;n=100;a=2;b=3;mu=2.5;
# res.MLEs <- get_MLE_biGamma.4Poolings(method, N, n, a, b, mu, start.point=c(a,b))

########################################################################################################################
# asym. Variance of MLEs

asymVar.randSamp <- function(N, n, a, b, mu) {
  ####################################################################
  # asym. Var of MLEs using random sampling design
  # if N = n, it is equivalent to full data
  ####################################################################
  omega <- n/N
  common.coef <- 1/omega * mu^2 / ((a+b)*trigamma(a)*trigamma(b) - trigamma(a) - trigamma(b))
  a11 <- common.coef * ((a+b)*trigamma(b) - 1) / mu^2
  a12 <- common.coef * 1 / mu^2 ; a21 <- a12
  a13 <- common.coef * (-1) * trigamma(b) / mu; a31 <- a13
  a22 <- common.coef * ((a+b)*trigamma(a) - 1) / mu^2
  a23 <- common.coef * (-1) * trigamma(a) / mu; a32 <- a23
  a33 <- common.coef * trigamma(a)*trigamma(b)
  mat <- matrix(c(a11, a12, a13, a21, a22, a23, a31, a32,a33), 3, 3)
  return(mat)
}

asymVar.fullData <- function(N, a, b, mu) {
  mat <- asymVar.randSamp(N=N, n=N, a=a, b=b, mu=mu)
  return(mat)
}

asymVar.equalPool <- function(N, n, a, b, mu) {
  ####################################################################
  # asym. Var of MLEs using equal-pooling design
  ####################################################################
  p <- floor(N/n)
  common.coef.denom <- p^3 * (p*(a+b)*trigamma(p*a)*trigamma(p*b) - trigamma(p*a) - trigamma(p*b))
  common.coef.numerator <- mu^2
  common.coef <- common.coef.numerator / common.coef.denom
  #det.val <- 1 / common.coef
  a11 <- common.coef * (p*(a+b)*trigamma(p*b) - 1) * p^2 / mu^2
  a12 <- common.coef * p^2 / mu^2 ; a21 <- a12
  a13 <- common.coef * (-1) * p^3 * trigamma(p*b) / mu; a31 <- a13
  a22 <- common.coef * (p*(a+b)*trigamma(p*a) - 1) * p^2 / mu^2
  a23 <- common.coef * (-1) * p^3 * trigamma(p*a) / mu; a32 <- a23
  a33 <- common.coef * p^4 * trigamma(p*a)*trigamma(p*b)
  mat <- matrix(c(a11, a12, a13, a21, a22, a23, a31, a32,a33), 3, 3)
  return(mat)
}


asymVar.onePool <- function(N, n, a, b, mu) {
  ####################################################################
  # asym. Var of MLEs using one-pool design
  # N: total sample size; n: size of one-pooled data
  ####################################################################
  # group size p
  p = N - n + 1; omega = p / N
  # entries of matrix
  fun.pattern.a <- omega*p*trigamma(p*a) + (1-omega)*trigamma(a)
  fun.pattern.b <- omega*p*trigamma(p*b) + (1-omega)*trigamma(b)
  common.coef.denom <- (a+b)*fun.pattern.a*fun.pattern.b - fun.pattern.a - fun.pattern.b
  common.coef.numerator <- mu^2
  common.coef <- common.coef.numerator / common.coef.denom
  #det.val <- 1 / common.coef
  a11 <- common.coef * ((a+b)*fun.pattern.b - 1) / mu^2
  a12 <- common.coef * 1 / mu^2 ; a21 <- a12
  a13 <- common.coef * (-1) * fun.pattern.b / mu; a31 <- a13
  a22 <- common.coef * ((a+b)*fun.pattern.a - 1) / mu^2
  a23 <- common.coef * (-1) * fun.pattern.a / mu; a32 <- a23
  a33 <- common.coef * fun.pattern.a * fun.pattern.b
  mat <- matrix(c(a11, a12, a13, a21, a22, a23, a31, a32,a33), 3, 3)
  return(mat)
}


asymVar.Monotone <- function(N, n, a, b, mu) {
  ####################################################################
  # asym. Var of MLEs using Monotone Group Size design
  ####################################################################
  # group size p
  p.vec <- getGroupSize.Montone(N = N, n=n)
  omega.vec <- p.vec / N
  fun.pattern.a <- sum(omega.vec * p.vec * trigamma(p.vec*a))
  fun.pattern.b <- sum(omega.vec * p.vec * trigamma(p.vec*b))
  # entries of matrix
  common.coef.denom <- (a+b)*fun.pattern.a*fun.pattern.b - fun.pattern.a - fun.pattern.b
  common.coef.numerator <- mu^2
  common.coef <- common.coef.numerator / common.coef.denom
  #det.val <- 1 / common.coef
  a11 <- common.coef * ((a+b)*fun.pattern.b - 1) / mu^2
  a12 <- common.coef * 1 / mu^2 ; a21 <- a12
  a13 <- common.coef * (-1) * fun.pattern.b / mu; a31 <- a13
  a22 <- common.coef * ((a+b)*fun.pattern.a - 1) / mu^2
  a23 <- common.coef * (-1) * fun.pattern.a / mu; a32 <- a23
  a33 <- common.coef * fun.pattern.a * fun.pattern.b
  mat <- matrix(c(a11, a12, a13, a21, a22, a23, a31, a32,a33), 3, 3)
  return(mat)
}

asymVar.4Poolings <- function(method, N, n, a, b, mu) {
  ##################################################################################
  # method: different numbers represent pooling methods:
  #         method 0: "Full data"
  #         method 1: "random sampling design"; method 2: "equal-pool design"; 
  #         method 3: "one-pool design"; method 4: "monontonically increasing design".
  ###################################################################################
  if (method == 0) { mat <- asymVar.fullData(N=N, a=a, b=b, mu=mu)}
  else if (method == 1) {mat <- asymVar.randSamp(N=N, n=n, a=a, b=b, mu=mu)}
  else if (method == 2) {mat <- asymVar.equalPool(N=N, n=n, a=a, b=b, mu=mu)}
  else if (method == 3) {mat <- asymVar.onePool(N=N, n=n, a=a, b=b, mu=mu)}
  else if (method == 4) {mat <- asymVar.Monotone(N=N, n=n, a=a, b=b, mu=mu)}
  return(mat)
}



######################################################################################################
# closed-form AUC

getAUC.biGamma.muxGTmuy <- function(a.x, b.x, mu.x, a.y, b.y, mu.y) {
  ############################################################################
  # This function is used only when a.x - a.y + b.x - b.y + 1 > 0 and mu.x > mu.y
  ############################################################################
  
  get_Sol <- function(t, z, a.x, b.x, mu.x, a.y, b.y, mu.y) {
    # get (only one) solution of the function g(v;theta)=c(t,z,theta)
    # where g(x)=x^(a.x-a.y+b.x-b.y) * exp((mu.x - mu.y)/(mu.x*mu.y)*x)
    
    # c function
    denom <- t; nom <- z^(a.x - a.y) * (1-z)^(b.x - b.y) 
    c.fun <- denom / nom
    
    # obj.fun: g(x) - c.fun
    obj.fun <- function(x) {
      obj.val <- x^(a.x-a.y+b.x-b.y)*exp((mu.x - mu.y)/(mu.x*mu.y)*x) - c.fun
      return(obj.val)
    }
    # find root via uniroot function
    root.val <- uniroot(obj.fun, lower = 0, upper = 1e3, extendInt = "yes")$root 
    return(root.val)
  }
  
  
  get_Prob_hGTt.given_t <- function(t, a.x, b.x, mu.x, a.y, b.y, mu.y) {
    ####################################################################################
    # Pr{h() >= t} for a given t, which is equivalent to
    # find the integral of [1 - F_Gamma(s(z1, t))] * f_Beta(z1) where z1 from 0 to 1;
    # F_Gamma is CDF of Gamma distribution with shape a.x+b.x and scale mu.x;
    # f_Beta is PDF of Beta distribution with shape1=a.x and shape2=b.x.
    ####################################################################################
    
    get_1_Minus_FGamma <- function(t, z, a.x, b.x, mu.x, a.y, b.y, mu.y) {
      # 1 - F_Gamma(s)
      res.root <- get_Sol(t=t, z=z, a.x=a.x, b.x=b.x, mu.x=mu.x, 
                          a.y=a.y, b.y=b.y, mu.y=mu.y)
      out <- 1 - pgamma(res.root, shape=a.x+b.x, scale = mu.x)
      return(out)
    }
    
    integrand <- function(z) {
      val.1_Minus_FGamma <- get_1_Minus_FGamma(t = t, z = z,
                                               a.x = a.x, b.x = b.x, mu.x = mu.x,
                                               a.y = a.y, b.y = b.y, mu.y = mu.y)
      out <- val.1_Minus_FGamma * dbeta(z, shape1 = a.x, shape2 = b.x)
      return(out)
    }
    integrand1.vec <- Vectorize(integrand, vectorize.args = "z")
    int.val <- integrate(integrand1.vec, lower = 0, upper = 1, subdivisions = 2000)$value
    return(int.val)
  }
  
  
  get_density_TIMES_partialDeriv <- function(t, z, a.x, b.x, mu.x, a.y, b.y, mu.y) {
    #############################################################################################
    # this function calculates f_Gamma(s) * deriv of S
    #############################################################################################
    get_partialDeriv.S.wrt.t <- function(root, z, a.x, b.x, mu.x, a.y, b.y, mu.y) {
      # this inner function "get_partialDeriv.S.wrt.t" is to calculate the partial 
      # derivative of S w.r.t t, where S is the solution of the equation g() = c.fun
      deriv.g <- function(s) {
        # this inner function "deriv.g" calculates the derivative of g(x), i.e. g'(x)
        val <- (a.x-a.y+b.x-b.y) * s^(a.x-a.y+b.x-b.y-1) * exp((mu.x-mu.y)/(mu.x*mu.y)*s) +
          (mu.x-mu.y)/(mu.x*mu.y) * s^(a.x-a.y+b.x-b.y) * exp((mu.x-mu.y)/(mu.x*mu.y)*s)
        return(val)
      }
      # "deriv.C_fun": deriv of C_fun w.r.t. t
      deriv.C_fun <- 1 / (z^(a.x - a.y) * (1-z)^(b.x - b.y)) 
      out <- deriv.C_fun / deriv.g(root)
      return(out)
    }
    
    root.val <- get_Sol(t=t, z=z, a.x=a.x, b.x=b.x, mu.x=mu.x, 
                        a.y=a.y, b.y=b.y, mu.y=mu.y)
    density_TIMES_partialDeriv <- dgamma(root.val, shape = a.y+b.y, scale = mu.y) * 
      get_partialDeriv.S.wrt.t(root=root.val, z=z, a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y)
    return(density_TIMES_partialDeriv)
  }
  
  get_deriv_Prob_hLEt.given_t <- function(t, a.x, b.x, mu.x, a.y, b.y, mu.y) {
    #################################################################################
    # dPr{h(Y) <= t}
    # find the integral of [f_Gamma(s)*deriv_of_s * f_Beta(z) where z from 0 to 1
    #################################################################################
    integrand <- function(z) {
      densityTIMESpartialDeriv <- get_density_TIMES_partialDeriv(t = t, z = z,
                                                                 a.x = a.x, b.x = b.x, mu.x = mu.x,
                                                                 a.y = a.y, b.y = b.y, mu.y = mu.y)
      out <- densityTIMESpartialDeriv * dbeta(z, shape1 = a.y, shape2 = b.y)
      #print(paste0("out: ", out, " ------"))
      return(out)
    }
    integrand.vec <- Vectorize(integrand, vectorize.args = "z")
    int.val <- try(integrate(integrand.vec, lower = 0, upper = 1, subdivisions = 2000)$value)
    return(int.val)
  }
  
  AUC.get_Prob_LX.GT.LY <- function(a.x, b.x, mu.x, a.y, b.y, mu.y) {
    #################################################################################
    # AUC: Pr{L(X1,X2; theta) >= L(Y1,Y2, theta)}
    #################################################################################
    integrand <- function(t) {
      #print(paste0("t:", t))
      # out1: Pr{h(X) > t}
      out1 <- get_Prob_hGTt.given_t(t = t, a.x = a.x, b.x = b.x, mu.x = mu.x, a.y = a.y, b.y = b.y, mu.y = mu.y)
      # out2: dPr{h(Y) < t}
      out2 <- get_deriv_Prob_hLEt.given_t(t = t, a.x = a.x, b.x = b.x, mu.x = mu.x, a.y = a.y, b.y = b.y, mu.y = mu.y)
      out <- out1 * out2
      #print(paste0("t=",  t," ||out1: ", out1, " ||out2: ",out2, " ||out: ", out))
      #print("-----------")
      return(out)
    }
    
    integrand.vec <- Vectorize(integrand, vectorize.args = "t")
    int.val <- integrate(integrand.vec, lower = 0, upper = Inf)$value
    return(int.val)
    
  }
  
  auc <- AUC.get_Prob_LX.GT.LY(a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y)
  return(auc)
}

# e.g.
# a.x <- 3.8; b.x <- 4.3; mu.x <- 3.8;
# a.y <- 3.2; b.y <- 3.5; mu.y <- 3.5
# auc <- getAUC.biGamma.muxGTmuy(a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y)
# f <- function(x) {
#   a.x <- x[1]; b.x <- x[2]; mu.x <- x[3]
#   a.y <- x[4]; b.y <- x[5]; mu.y <- x[6]
#   return(getAUC.biGamma.muxGTmuy(a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y))
# }
# start_time <- Sys.time()
# grad(f,x=c(a.x, b.x, mu.x, a.y, b.y, mu.y), method = "simple")
# end_time <- Sys.time()
# end_time - start_time


######################################################################################################
# AUC based on non-parametric method

getAUC.biGamma.nonParametric <- function(a.x, b.x, mu.x, a.y, b.y, mu.y) {
  ##########################################################################
  # get AUC based on non-parametric method
  ##########################################################################
  MC <- 300; NN.x <- 1000; NN.y <- 1000
  
  get_RV_biGamma <- function(n, shape1, shape2, scalePar) {
    #############################################################
    # shape1: a, shape2: b, scalePar: mu
    # returns random vectors which follows McKay's bi-Gamma
    #############################################################
    mdata <- data.frame(w1 = rgamma(nn <- n, shape=shape1, scale=scalePar))
    mdata <- transform(mdata, zedd = rgamma(nn, shape=shape2, scale=scalePar))
    mdata <- transform(mdata, w2 = w1 + zedd) # Z is defined as W2- w1 | W1=w
    #print(paste0("simulation ", dim(mdata)[1], " bi-Gamma samples"))
    return(mdata[, c("w1", "w2")])
  }
  
  best.combination.fun <- function(w1, w2, a1, b1, mu1, a2, b2, mu2) {
    #####################################################################
    # L(w1,w2; a1,b1,mu1,a2,b2,mu2)
    #####################################################################
    denom <- mu2^(a2+b2)*gamma(a2)*gamma(b2)
    nom <- mu1^(a1+b1)*gamma(a1)*gamma(b1)
    L.fun <- denom / nom * w1^(a1-a2) * (w2-w1)^(b1-b2) * exp((mu1-mu2)/(mu1*mu2)*w2)
    return(L.fun)
  }
  
  best.combination.fun.vector <- function(data, a1, b1, mu1, a2, b2, mu2) {
    inner.fun <- function(i) {best.combination.fun(w1 = data[i,1], w2 = data[i,2],
                                                   a1=a1,b1=b1,mu1=mu1,a2=a2,b2=b2,mu2=mu2)}
    n.sample <- dim(data)[1]
    return(sapply(1:n.sample,inner.fun))
  }
  
  getAUC.mannWhitney <- function(NN.x, NN.y, a.x, b.x, mu.x, a.y, b.y, mu.y, seed.num=1) {
    # get AUC using mann-whitney stat
    set.seed(seed.num)
    dat.x <- get_RV_biGamma(NN.x, shape1=a.x, shape2=b.x, scalePar=mu.x)
    set.seed(seed.num)
    dat.y <- get_RV_biGamma(NN.y, shape1=a.y, shape2=b.y, scalePar=mu.y)
    # L(X)
    L_x <- best.combination.fun.vector(data = dat.x, a1=a.x,b1=b.x,mu1=mu.x,a2=a.y, b2=b.y, mu2=mu.y)
    # L(Y)
    L_y <- best.combination.fun.vector(data = dat.y, a1=a.x,b1=b.x,mu1=mu.x,a2=a.y,b2=b.y, mu2=mu.y)
    test.res <- wilcox.test(L_x, L_y)
    return( test.res$statistic / (NN.x * NN.y))
  }
  
  
  res.AUC.MC <- rep(NA, MC) 
  for (jj in 1:MC){
    # if (jj %% floor(MC/20) == 0) {
    #   print(paste0("Job i: ", jj))
    # }
    res.AUC.MC[jj] <- getAUC.mannWhitney(NN.x, NN.y, a.x, b.x, mu.x, a.y, b.y, mu.y, seed.num = jj)
  }
  
  auc.nonParam <- mean(res.AUC.MC)
  return(auc.nonParam)
  
}

# # e.g.
# N.x <- 1000; N.y <- 1000; MC <- 500
# a.x <- 3.2; b.x <- 3.5; mu.x <- 3.5;
# a.y <- 3.2; b.y <- 3.5; mu.y <- 3.5
# getAUC.biGamma.nonParametric(a.x, b.x, mu.x, a.y, b.y, mu.y)




######################################################################################################
# asymVar of AUC (true value)
get.AsymVar_AUC.TrueParam <- function(method.pooling, N.x, N.y, n.x, n.y,
                                      a.x, b.x, mu.x, a.y, b.y, mu.y) {
  #########################################################################################
  # this function is used to calculate the asym Variance plugging in true parameters
  # method.pooling: different numbers represent pooling methods:
  #         method 0: "Full data"
  #         method 1: "random sampling design"; method 2: "equal-pool design"; 
  #         method 3: "one-pool design"; method 4: "monontonically increasing design".
  #########################################################################################
  
  # part 1: AUC
  AUC.true <- try(getAUC.biGamma.muxGTmuy(a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y))
  if  (inherits(AUC.true, "try-error")) {
    print(paste("method=", method, "with params=", a.x, b.x, mu.x, a.y, b.y, mu.y, "closed-form AUC failed, using non-param AUC instead"))
    AUC.true <- getAUC.biGamma.nonParametric(a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y)
  }
  
  # part 2: partial deriv of AUC with respect to all parameters
  f <- function(x) {
    a.x <- x[1]; b.x <- x[2]; mu.x <- x[3]
    a.y <- x[4]; b.y <- x[5]; mu.y <- x[6]
    auc.with_param <- try(getAUC.biGamma.muxGTmuy(a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y))
    if  (inherits(auc.with_param, "try-error"))  {
      auc.with_param <- getAUC.biGamma.nonParametric(a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y)
    }
    return(auc.with_param)
  }
  partial.deriv <- grad(f,x=c(a.x, b.x, mu.x, a.y, b.y, mu.y), method = "simple")
  #partial.deriv <- grad(f,x=c(a.x, b.x, mu.x, a.y, b.y, mu.y))
  
  # part 3: joint asym var of MLEs of 6 parameters
  # asym Var For Diseased group X
  asym.Var.X.true <- asymVar.4Poolings(method=method.pooling, N=N.x, n=n.x, a=a.x, b=b.x, mu=mu.x)
  # asym Var For non-Diseased group Y
  asym.Var.Y.true <- asymVar.4Poolings(method=method.pooling, N=N.y, n=n.y, a=a.y, b=b.y, mu=mu.y)
  # stack them
  N.total <- N.x + N.y; prop.p <- N.x / N.total
  joint.asym.var.true <- bdiag(asym.Var.X.true / prop.p, asym.Var.Y.true / (1-prop.p) )
  
  # part 4: asym.Var of MLEs of AUC (plug in true values)
  asym.Var.AUC <- t(partial.deriv) %*% joint.asym.var.true %*% partial.deriv
  asym.Var.AUC <- as.numeric(asym.Var.AUC)
  return(asym.Var.AUC)
}




# e.g.
# method.pooling <- 0
# N.x <- 1000; N.y <- 1000; n.x <- 100; n.y <- 100
# a.x <- 4.3; b.x <- 4.6; mu.x <- 5.2;
# a.y <- 3.2; b.y <- 3.5; mu.y <- 3.5
# start_time <- Sys.time()
# get.AsymVar_AUC.TrueParam(method.pooling, N.x, N.y, n.x, n.y,
#                                       a.x, b.x, mu.x, a.y, b.y, mu.y)
# end_time <- Sys.time()
# end_time - start_time



# # asymVar of AUC (MLE)
# get.AsymVar_AUC.plugInMLE <- function(method.pooling, N.x, N.y, n.x, n.y,
#                                       a.X.MLE, b.X.MLE, mu.X.MLE,  
#                                       a.Y.MLE, b.Y.MLE, mu.Y.MLE) {
#   #########################################################################################
#   # this function is used to calculate the asym Variance plugging in true parameters
#   # method.pooling: different numbers represent pooling methods:
#   #         method 0: "Full data"
#   #         method 1: "random sampling design"; method 2: "equal-pool design"; 
#   #         method 3: "one-pool design"; method 4: "monontonically increasing design".
#   #########################################################################################
#   
#   # a.x - a.y + b.x - b.y + 1 > 0 and mu.x > mu.y
#   if ((a.X.MLE + b.X.MLE - a.Y.MLE - b.Y.MLE + 1) <= 0 | mu.x <= mu.y) {
#     print("Does not satisfy a.x - a.y + b.x - b.y + 1 > 0 and mu.x > mu.y!")
#     return(return(list("AUC.hat" = NA, "asymVar.hat"=NA)))
#   }
#   
#   # part 1: AUC
#   AUC.MLE <- try(getAUC.biGamma.muxGTmuy(a.x=a.X.MLE, b.x=b.X.MLE, mu.x=mu.X.MLE, 
#                                           a.y=a.Y.MLE, b.y=b.Y.MLE, mu.y=mu.Y.MLE))
#   if  (inherits(AUC.MLE, "try-error")) {
#     print(paste0("method=", method, "with params=", a.X.MLE, b.X.MLE, mu.X.MLE,
#                 a.Y.MLE, b.Y.MLE, mu.Y.MLE, "closed-form AUC failed, using non-param AUC instead"))
#     AUC.MLE <- getAUC.biGamma.nonParametric(a.x=a.X.MLE, b.x=b.X.MLE, mu.x=mu.X.MLE, 
#                                              a.y=a.Y.MLE, b.y=b.Y.MLE, mu.y=mu.Y.MLE)
#   }
#   # part 2: partial deriv of AUC with respect to all parameters
#   f <- function(x) {
#     a.x <- x[1]; b.x <- x[2]; mu.x <- x[3]
#     a.y <- x[4]; b.y <- x[5]; mu.y <- x[6]
#     auc.with_param <- try(getAUC.biGamma.muxGTmuy(a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y))
#     if  (inherits(auc.with_param, "try-error"))  {
#       auc.with_param <- getAUC.biGamma.nonParametric(a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y)
#     }
#     return(auc.with_param)
#   }
#   print("---start--pd")
#   partial.deriv <- grad(f,x=c(a.x=a.X.MLE, b.x=b.X.MLE, mu.x=mu.X.MLE, 
#                               a.y=a.Y.MLE, b.y=b.Y.MLE, mu.y=mu.Y.MLE), method = "simple")
#   print("---end--pd")
#   # print(partial.deriv)
#   #partial.deriv <- grad(f,x=c(a.x, b.x, mu.x, a.y, b.y, mu.y))
#   
#   # part 3: joint asym var of MLEs of 6 parameters
#   # asym Var For Diseased group X
#   asym.Var.X.MLE <- asymVar.4Poolings(method=method.pooling, N=N.x, n=n.x, a=a.X.MLE, b=b.X.MLE, mu=mu.X.MLE)
#   # asym Var For non-Diseased group Y
#   asym.Var.Y.MLE <- asymVar.4Poolings(method=method.pooling, N=N.y, n=n.y, a=a.Y.MLE, b=b.Y.MLE, mu=mu.Y.MLE)
#   # stack them
#   N.total <- N.x + N.y; prop.p <- N.x / N.total
#   joint.asym.var.MLE <- bdiag(asym.Var.X.MLE / prop.p, asym.Var.Y.MLE / (1-prop.p) )
#   # print(joint.asym.var.MLE)
#   # part 4: asym.Var of MLEs of AUC (plug in true values)
#   asym.Var.AUC <- t(partial.deriv) %*% joint.asym.var.MLE %*% partial.deriv
#   asym.Var.AUC <- as.numeric(asym.Var.AUC)
#   return(list("AUC.hat" = AUC.MLE, "asymVar.hat"=asym.Var.AUC))
# }


get.AsymVar_AUC.plugInMLE <- function(method.pooling, N.x, N.y, n.x, n.y,
                                      a.X.MLE, b.X.MLE, mu.X.MLE,  
                                      a.Y.MLE, b.Y.MLE, mu.Y.MLE) {
  #########################################################################################
  # this function is used to calculate the asym Variance plugging in true parameters
  # method.pooling: different numbers represent pooling methods:
  #         method 0: "Full data"
  #         method 1: "random sampling design"; method 2: "equal-pool design"; 
  #         method 3: "one-pool design"; method 4: "monontonically increasing design".
  #########################################################################################

  # a.x - a.y + b.x - b.y + 1 > 0 and mu.x > mu.y
  if ((a.X.MLE + b.X.MLE - a.Y.MLE - b.Y.MLE + 1) <= 0 | mu.X.MLE <= mu.Y.MLE) {
    print("Does not satisfy a.x - a.y + b.x - b.y + 1 > 0 and mu.x > mu.y!, we will use non-param AUC instead")
    AUC.MLE <-  getAUC.biGamma.nonParametric(a.x=a.X.MLE, b.x=b.X.MLE, mu.x=mu.X.MLE, 
                                 a.y=a.Y.MLE, b.y=b.Y.MLE, mu.y=mu.Y.MLE)
    return(return(list("AUC.hat" = AUC.MLE, "asymVar.hat"=NA)))
  }
  
  # part 1: AUC
  AUC.MLE <- try(getAUC.biGamma.muxGTmuy(a.x=a.X.MLE, b.x=b.X.MLE, mu.x=mu.X.MLE, 
                                         a.y=a.Y.MLE, b.y=b.Y.MLE, mu.y=mu.Y.MLE))
  if  (inherits(AUC.MLE, "try-error")) {
    print(paste("method=", method.pooling, "with params=", a.X.MLE, b.X.MLE, mu.X.MLE,
                 a.Y.MLE, b.Y.MLE, mu.Y.MLE, "closed-form AUC failed, using non-param AUC instead"))
    AUC.MLE <- getAUC.biGamma.nonParametric(a.x=a.X.MLE, b.x=b.X.MLE, mu.x=mu.X.MLE, 
                                            a.y=a.Y.MLE, b.y=b.Y.MLE, mu.y=mu.Y.MLE)
    return(return(list("AUC.hat" = AUC.MLE, "asymVar.hat"=NA)))
  }
  # part 2: partial deriv of AUC with respect to all parameters
  f <- function(x) {
    a.x <- x[1]; b.x <- x[2]; mu.x <- x[3]
    a.y <- x[4]; b.y <- x[5]; mu.y <- x[6]
    auc.with_param <- try(getAUC.biGamma.muxGTmuy(a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y))
    # if  (inherits(auc.with_param, "try-error"))  {
    #   auc.with_param <- getAUC.biGamma.nonParametric(a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y)
    # }
    return(auc.with_param)
  }
  print("---start--pd")
  partial.deriv <- try(grad(f,x=c(a.x=a.X.MLE, b.x=b.X.MLE, mu.x=mu.X.MLE, 
                              a.y=a.Y.MLE, b.y=b.Y.MLE, mu.y=mu.Y.MLE), method = "simple"))
  if (inherits(partial.deriv, "try-error")) {
    print("partial deriv failed, will use true asym Var instead")
    return(return(list("AUC.hat" = AUC.MLE, "asymVar.hat"=NA)))
  }
  print("---end--pd")
  # print(partial.deriv)
  #partial.deriv <- grad(f,x=c(a.x, b.x, mu.x, a.y, b.y, mu.y))
  
  # part 3: joint asym var of MLEs of 6 parameters
  # asym Var For Diseased group X
  asym.Var.X.MLE <- asymVar.4Poolings(method=method.pooling, N=N.x, n=n.x, a=a.X.MLE, b=b.X.MLE, mu=mu.X.MLE)
  # asym Var For non-Diseased group Y
  asym.Var.Y.MLE <- asymVar.4Poolings(method=method.pooling, N=N.y, n=n.y, a=a.Y.MLE, b=b.Y.MLE, mu=mu.Y.MLE)
  # stack them
  N.total <- N.x + N.y; prop.p <- N.x / N.total
  joint.asym.var.MLE <- bdiag(asym.Var.X.MLE / prop.p, asym.Var.Y.MLE / (1-prop.p) )
  # print(joint.asym.var.MLE)
  # part 4: asym.Var of MLEs of AUC (plug in true values)
  asym.Var.AUC <- t(partial.deriv) %*% joint.asym.var.MLE %*% partial.deriv
  asym.Var.AUC <- as.numeric(asym.Var.AUC)
  return(list("AUC.hat" = AUC.MLE, "asymVar.hat"=asym.Var.AUC))
}




######################################################################################################
# variance of the estimated AUC via MC evalutions
get.Var_AUC.MC <- function(method.pooling, N.x, N.y, n.x, n.y,
                           a.x, b.x, mu.x, a.y, b.y, mu.y, start.MC, end.MC) {
  #########################################################################################
  # this function is used to approximate the Variance of the estimated AUC via
  # extensive MC evalutions based on different pooling designs.
  # method.pooling: different numbers represent pooling methods:
  #         method 0: "Full data"
  #         method 1: "random sampling design"; method 2: "equal-pool design"; 
  #         method 3: "one-pool design"; method 4: "monontonically increasing design".
  # start.MC: start_index of MC; end.MC: end_index of MC
  #########################################################################################
  
  AUC.true <- try(getAUC.biGamma.muxGTmuy(a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y))
  if  (inherits(AUC.true, "try-error")) {
    print(paste("method=", method, "with params=", a.x, b.x, mu.x, a.y, b.y, mu.y, "closed-form AUC failed, using non-param AUC instead"))
    # here might be problematic
    AUC.true <- getAUC.biGamma.nonParametric(a.x=a.x, b.x=b.x, mu.x=mu.x, a.y=a.y, b.y=b.y, mu.y=mu.y)
  }
  
  N.total = N.x + N.y
   
  AUC.hat.array <- c()
  for (seed.num in start.MC:end.MC) {
    set.seed(seed.num)
    # res.X.MLEs <- get_MLE_biGamma.4Poolings(method=method.pooling, N=N.x, n=n.x, a=a.x, b=b.x, mu=mu.x, start.point=c(a.x+rnorm(1,0,0.2),b.x+rnorm(1,0,0.2)))
    res.X.MLEs <- get_MLE_biGamma.4Poolings(method=method.pooling, N=N.x, n=n.x, a=a.x, b=b.x, mu=mu.x, start.point=c(1,1))
    a.X.MLE = res.X.MLEs$a.hat; b.X.MLE = res.X.MLEs$b.hat; mu.X.MLE = res.X.MLEs$mu.hat
    # for Y
    # set.seed(seed.num)
    # res.Y.MLEs <- get_MLE_biGamma.4Poolings(method=method.pooling, N=N.y, n=n.y, a=a.y, b=b.y, mu=mu.y, start.point=c(a.y+rnorm(1,0,0.2),b.y+rnorm(1,0,0.2)))
    res.Y.MLEs <- get_MLE_biGamma.4Poolings(method=method.pooling, N=N.y, n=n.y, a=a.y, b=b.y, mu=mu.y, start.point=c(1,1))
    a.Y.MLE = res.Y.MLEs$a.hat; b.Y.MLE = res.Y.MLEs$b.hat; mu.Y.MLE = res.Y.MLEs$mu.hat
    
    
    AUC.hat <- try(getAUC.biGamma.muxGTmuy(a.x=a.X.MLE, b.x=b.X.MLE, mu.x=mu.X.MLE, 
                                           a.y=a.Y.MLE, b.y=b.Y.MLE, mu.y=mu.Y.MLE))
    if  (inherits(AUC.hat, "try-error")) {
      AUC.hat <- getAUC.biGamma.nonParametric(a.x=a.X.MLE, b.x=b.X.MLE, mu.x=mu.X.MLE, 
                                              a.y=a.Y.MLE, b.y=b.Y.MLE, mu.y=mu.Y.MLE)
    }
    AUC.hat.array <- c(AUC.hat.array, AUC.hat)
    print(c(seed.num,  N.total * mean((AUC.hat.array - AUC.true)^2)))
  }
  
  # singleRun <- function(seed.num) {
  #   # for X
  #   print(paste("seed num=", seed.num))
  #   set.seed(seed.num)
  #   # res.X.MLEs <- get_MLE_biGamma.4Poolings(method=method.pooling, N=N.x, n=n.x, a=a.x, b=b.x, mu=mu.x, start.point=c(a.x+rnorm(1,0,0.2),b.x+rnorm(1,0,0.2)))
  #   res.X.MLEs <- get_MLE_biGamma.4Poolings(method=method.pooling, N=N.x, n=n.x, a=a.x, b=b.x, mu=mu.x, start.point=c(1,1))
  #   a.X.MLE = res.X.MLEs$a.hat; b.X.MLE = res.X.MLEs$b.hat; mu.X.MLE = res.X.MLEs$mu.hat
  #   # for Y
  #   # set.seed(seed.num)
  #   # res.Y.MLEs <- get_MLE_biGamma.4Poolings(method=method.pooling, N=N.y, n=n.y, a=a.y, b=b.y, mu=mu.y, start.point=c(a.y+rnorm(1,0,0.2),b.y+rnorm(1,0,0.2)))
  #   res.Y.MLEs <- get_MLE_biGamma.4Poolings(method=method.pooling, N=N.y, n=n.y, a=a.y, b=b.y, mu=mu.y, start.point=c(1,1))
  #   a.Y.MLE = res.Y.MLEs$a.hat; b.Y.MLE = res.Y.MLEs$b.hat; mu.Y.MLE = res.Y.MLEs$mu.hat
  # 
  #   
  #   AUC.hat <- try(getAUC.biGamma.muxGTmuy(a.x=a.X.MLE, b.x=b.X.MLE, mu.x=mu.X.MLE, 
  #                                          a.y=a.Y.MLE, b.y=b.Y.MLE, mu.y=mu.Y.MLE))
  #   if  (inherits(AUC.hat, "try-error")) {
  #     AUC.hat <- getAUC.biGamma.nonParametric(a.x=a.X.MLE, b.x=b.X.MLE, mu.x=mu.X.MLE, 
  #                                             a.y=a.Y.MLE, b.y=b.Y.MLE, mu.y=mu.Y.MLE)
  #   }
  #   # print(c(AUC.hat, a.X.MLE, b.X.MLE, mu.X.MLE, a.Y.MLE, b.Y.MLE, mu.Y.MLE))
  #   return(AUC.hat)
  # }
  # 
  # records <- sapply(start.MC:end.MC, singleRun)
  # res = N.total * mean((records - AUC.true)^2)
  # return(list("var"=res, "records"=records))
  
  res = N.total * mean((AUC.hat.array - AUC.true)^2)
  return(list("var"=res, "AUC.mle.array"=AUC.hat.array))
}

# method.pooling = 2; start.MC=20; end.MC = 80
# N.x <- 1000; N.y <- 1000; n.x <- 100; n.y <- 100
# a.x <- 4.3; b.x <- 4.6; mu.x <- 5.2;
# a.y <- 3.2; b.y <- 3.5; mu.y <- 3.5
# r <- get.Var_AUC.MC(method.pooling, N.x, N.y, n.x, n.y,
#               a.x, b.x, mu.x, a.y, b.y, mu.y, start.MC, end.MC)
# 
# 
# start_time <- Sys.time()
# get.AsymVar_AUC.TrueParam(method.pooling, N.x, N.y, n.x, n.y,
#                           a.x, b.x, mu.x, a.y, b.y, mu.y)
# end_time <- Sys.time()
# end_time - start_time



