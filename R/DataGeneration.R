data_generation <- function(n, a, b, c, d, e, f, g, h, i, j, seed){
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
    y.C = beta0 + x.C * beta1 + epsilon
    x.R = x.C * betastar + epsilonr
    y.R = y.C * betastar + epsilonr
    return(list(
        yl = y.C - y.R,
        yu = y.C + y.R,
        xl = x.C - x.R,
        xu = x.C + x.R
    ))
}

data_generation0 <- function(n, p, a, b, c, d, e, f, g, h, i, j, seed){
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
  y.C = beta0 + x.C * beta1 + epsilon
  x.R = x.C * betastar + epsilonr
  y.R = y.C * betastar + epsilonr
  n1 = floor(n*alpha)
  sampled.index = sample(1:n,n1)
  
  if(outlierType == "central")
    y.C[sampled.index] = 5
  else if(outlierType == "range")
    y.R[sampled.index] = 5
  else if(outlierType == "central and range")
  {
    y.C[sampled.index] = 5
    y.R[sampled.index] = 5
  }
  
  return(list(
    xc = x.C,
    xr = x.R,
    yc = y.C,
    yr = y.R,
    yl = y.C - y.R,
    yu = y.C + y.R,
    xl = x.C - x.R,
    xu = x.C + x.R
  ))
}

data_generation_outlier2 <- function(n, a, b, c, d, e, f, g, h, i, j, seed, outlierType="central", alpha=0.1){
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
    y.C = beta0 + x.C * beta1 + epsilon
    x.R = x.C * betastar + epsilonr
    y.R = y.C * betastar + epsilonr
    n1 = floor(n*alpha)
    sampled.index = sample(1:n,n1)
    
    if(outlierType == "central")
        y.C[sampled.index] = 1000
    else if(outlierType == "range")
        y.R[sampled.index] = 1000
    else if(outlierType == "central and range")
    {
        y.C[sampled.index] = 1000
        y.R[sampled.index] = 1000
    }
    
    return(list(
        xc = x.C,
        xr = x.R,
        yc = y.C,
        yr = y.R,
        yl = y.C - y.R,
        yu = y.C + y.R,
        xl = x.C - x.R,
        xu = x.C + x.R
    ))
}

data_generation_outlier3 <- function(n, a, b, c, d, e, f, g, h, i, j, seed, outlierType="central", alpha=0.1){
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
    y.C = beta0 + x.C * beta1 + epsilon
    x.R = x.C * betastar + epsilonr
    y.R = y.C * betastar + epsilonr
    n1 = floor(n*alpha)
    sampled.index = sample(1:n,n1)
    
    if(outlierType == "central")
        y.C[sampled.index] = runif(n1, 6, 12)
    else if(outlierType == "range")
        y.R[sampled.index] = runif(n1, 2, 4)
    else if(outlierType == "central and range")
    {
        y.C[sampled.index] = runif(n1, 6, 12)
        y.R[sampled.index] = runif(n1, 2, 4)
    }
    
    return(list(
        xc = x.C,
        xr = x.R,
        yc = y.C,
        yr = y.R,
        yl = y.C - y.R,
        yu = y.C + y.R,
        xl = x.C - x.R,
        xu = x.C + x.R
    ))
}

data_generation_outlier4 <- function(n, a, b, c, d, e, f, g, h, i, j, seed, outlierType="central", alpha=0.1){
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
  y.C = beta0 + x.C * beta1 + epsilon
  x.R = x.C * betastar + epsilonr
  y.R = y.C * betastar + epsilonr
  n1 = floor(n*alpha)
  sampled.index = sample(1:n,n1)
  
  if(outlierType == "central")
    y.C[sampled.index] = y.C[sampled.index] + 5*rt(n1,3)
  else if(outlierType == "range")
    y.R[sampled.index] = y.R[sampled.index] + 5*abs(rt(n1,3))
  else if(outlierType == "central and range")
  {
    y.C[sampled.index] = y.C[sampled.index] + 5*rt(n1,3)
    y.R[sampled.index] = y.R[sampled.index] + 5*abs(rt(n1,3))
  }
  
  return(list(
    xc = x.C,
    xr = x.R,
    yc = y.C,
    yr = y.R,
    yl = y.C - y.R,
    yu = y.C + y.R,
    xl = x.C - x.R,
    xu = x.C + x.R
  ))
}

data_generation_outlier_realdata <- function(data, seed, outlierType=1, alpha=0.1){
    n = nrow(data)
    n1 = floor(n * alpha)
    sampled.index = sample(1:n,n1)
    
    if(outlierType==1)
    {
        data[sampled.index,1] = 0.1
        data[sampled.index,2] = 0.2
    }
    else if(outlierType==2)
    {
        data[sampled.index,1] = 1000
        data[sampled.index,2] = 2000
    }
    else if(outlierType==3)
    {
        yc = runif(n1,0,1)
        yr = runif(n1,0,1)
        data[sampled.index,1] = yc-yr
        data[sampled.index,2] = yc+yr
    }
    else
    {
        yc = data[sampled.index,2]/2 - data[sampled.index,1]/2
        yr = data[sampled.index,2]/2 + data[sampled.index,1]/2
        yc = yc + 2*rt(n1,3)
        yr = yr + 2*abs(rt(n1,3))
        data[sampled.index,1] = yc-yr
        data[sampled.index,2] = yc+yr
    }
    
    return(data)
}