##### Hyper-parameter estimators functions of iETKRR
# S1
sigma2est = function(y1, y2, frac = .5) {
  n = length(y1)
  m = floor(n*frac)
  idx1 = sample(1:n, m, replace = T)
  idx2 = sample(1:n, m, replace = T)
  tmp = (y1[idx1] - y2[idx2])^2
  mean(quantile(tmp[tmp != 0], probs = c(.9, .1)))
}

# S2
sigma2est2.new = function(y1, y2) {
  D = outer(y1,y2,"-")
  D = D^2
  D.no.zero = D[ which(!D == 0)]
  median(D.no.zero)
}

##### Gaussian Kernel
gauss.kern = function(a, b, s)
{
  as.vector(exp(-(1/s)*(a-b)^2))
}


# ETKRR_S1
kernel.reg = function(x, y, tol = 1e-10, maxit = 100)
{
  n = nrow(x)
  # Initialization
  x = cbind(1, x)
  betahat = solve(t(x)%*%x)%*%t(x)%*%y
  yhat = x%*%betahat
  s2 = sigma2est(y, yhat)
  K = gauss.kern(y, yhat, s2)
  S = sum(2-2*K)
  it = 1
  # Model Step
  repeat {
    it = it+1
    betahat = solve(t(x)%*%diag(K)%*%x)%*%t(x)%*%diag(K)%*%y
    yhat = x%*%betahat
    K = gauss.kern(y, yhat, s2)
    S = c(S, sum(2-2*K))
    if (abs(S[it]-S[(it-1)]) <= tol || it >= maxit) break
  }
  (result = list(coef = as.vector(betahat), fitted = as.vector(yhat), criterion = S, weigth = K, iter = it))
}


# ETKRR_S2
kernel.reg2.new = function(x, y, tol = 1e-10, maxit = 100)
{
  x = as.matrix(x)
  n = nrow(x)
  # Initialization
  x = cbind(1, x)
  betahat = solve(t(x)%*%x)%*%t(x)%*%y
  yhat = x%*%betahat
  s2 = sigma2est2.new(y, yhat)
  K = gauss.kern(y, yhat, s2)
  S = sum(2-2*K)
  it = 1
  # Model Step
  repeat {
    it = it+1
    betahat = solve(t(x)%*%diag(K)%*%x)%*%t(x)%*%diag(K)%*%y
    yhat = x%*%betahat
    K = gauss.kern(y, yhat, s2)
    S = c(S, sum(2-2*K))
    if (abs(S[it]-S[(it-1)]) <= tol || it >= maxit) break
  }
  (result = list(coef = as.vector(betahat), fitted = as.vector(yhat), criterion = S, weigth = K, iter = it))
}

# ETKRR_S3
kernel.reg3 = function(x, y, tol = 1e-10, maxit = 100)
{
  x = as.matrix(x)
  n = nrow(x)
  # Initialization
  x = cbind(1, x)
  betahat = solve(t(x)%*%x)%*%t(x)%*%y
  yhat = x%*%betahat
  s2 = sum((y-yhat)^2)/(n-ncol(x))
  K = gauss.kern(y, yhat, s2)
  S = sum(2-2*K)
  it = 1
  # Model Step
  repeat {
    it = it+1
    betahat = solve(t(x)%*%diag(K)%*%x)%*%t(x)%*%diag(K)%*%y
    yhat = x%*%betahat
    K = gauss.kern(y, yhat, s2)
    S = c(S, sum(2-2*K))
    if (abs(S[it]-S[(it-1)]) <= tol || it >= maxit) break
  }
  (result = list(coef = as.vector(betahat), fitted = as.vector(yhat), criterion = S, weigth = K, iter = it))
}

# ETKRR_S4
kernel.reg4 = function(x, y, tol = 1e-10, maxit = 100)
{
  x = as.matrix(x)
  n = nrow(x)
  # Initialization
  x = cbind(1, x)
  betahat = solve(t(x)%*%x)%*%t(x)%*%y
  yhat = as.numeric(x%*%betahat)
  s2 = h.select(y, yhat, method="aicc")
  K = gauss.kern(y, yhat, s2)
  S = sum(2-2*K)
  it = 1
  # Model Step
  repeat {
    it = it+1
    betahat = solve(t(x)%*%diag(K)%*%x)%*%t(x)%*%diag(K)%*%y
    yhat = x%*%betahat
    K = gauss.kern(y, yhat, s2)
    S = c(S, sum(2-2*K))
    if (abs(S[it]-S[(it-1)]) <= tol || it >= maxit) break
  }
  (result = list(coef = as.vector(betahat), fitted = as.vector(yhat), criterion = S, weigth = K, iter = it))
}