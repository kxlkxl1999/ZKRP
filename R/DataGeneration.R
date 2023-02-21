data_generation <- function(n, p, a, b, c, d, e, f, g, h, i, j, seed){
  set.seed(seed)
  x.C <- matrix(NA, nrow=n, ncol=p)   # center point of x
  x.R <- matrix(NA, nrow=n, ncol=p)
  beta0 = c()
  beta1 = c()
  betastar = c()
  for(s in 1:p)
  {
    x.C[,s] = runif(n, a[s], b[s])
    beta1 = c(beta1, runif(1, c[s], d[s]))
    betastar = c(betastar, runif(1, g[s], h[s]))
  }
  beta0 = runif(n, c[1], d[1])
  epsilon = runif(n, e, f)
  epsilonr = runif(n, i, j)
  y.C = beta0 + x.C %*% diag(beta1) + epsilon
  x.R = x.C %*% diag(betastar) + epsilonr
  y.R = y.C %*% diag(betastar) + epsilonr
  return(list(
    xc = x.C,
    xr = x.R,
    yc = y.C,
    yr = y.R,
    xl = x.C - x.R,
    xu = x.C + x.R,
    yl = y.C - y.R,
    yu = y.C + y.R
  ))
}

data_generation_outlier1 <- function(n, a, b, c, d, e, f, g, h, i, j, seed, outlierType="central", alpha=0.1){
  set.seed(seed)
  x.C  = c()   # center point of x
  x.R  = c()
  beta0 = 0
  beta1 = 0
  betastar = 0
  x.C = runif(n, a, b)
  beta1 = runif(1, c, d)
  betastar = runif(1, g, h)
  beta0 = runif(n, c, d)
  epsilon = runif(n, e, f)
  epsilonr = runif(n, i, j)
  y.C = beta0 + x.C * diag(beta1) + epsilon
  x.R = x.C * betastar + epsilonr
  y.R = y.C * betastar + epsilonr
  n1 = floor(n*alpha)
  
  if(outlierType == "central")
  {
    
  }
  
  return(list(
    xc = x.C,
    xr = x.R,
    yc = y.C,
    yr = y.R,
    xl = x.C - x.R,
    xu = x.C + x.R,
    yl = y.C - y.R,
    yu = y.C + y.R
  ))
}