HF1 <- function(formula = formula, data = data)
{
  # model = rq(y~x, method = "pfn")
  # return(as.numeric(model$coefficients))
  
  # input formula for regression
  formula1 = as.formula(formula)
  model_interval = model.frame(formula = formula1, data = data)
  
  # response variable interval value
  yinterval = model.response(model_interval)
  y = as.matrix(yinterval)
  y.C = (y[, 1] + y[, 2]) / 2   # center point of y
  y.R = y[, 2] - y[, 1]         # range point of y
  
  # predictor variable interval matrix (each left: Lower, each right: Upper)
  xinterval = model.matrix(attr(model_interval, "terms"), data = data)
  x = as.matrix(xinterval)
  
  p = ncol(x) - 1 # num. of predictor variable except intercept
  n = nrow(x)
  var.names = colnames(x)[(1:{p/2})*2]
  
  x.C <- matrix(NA, nrow=n, ncol=p/2)   # center point of x
  colnames(x.C) <- var.names
  x.R <- matrix(NA, nrow=n, ncol=p/2)   # range point of x
  colnames(x.R) <- var.names
  for (j in 1:{p/2}) {
    for (i in 1:n) {
      x.C[i, j] = (x[i, 2*j] + x[i, 2*j+1]) / 2
      x.R[i, j] = (x[i, 2*j+1] - x[i, 2*j])
    }
  }
  
  model.C = rq(y.C~x.C, method = "pfn")
  model.R = rq(y.R~x.R, method = "pfn")
  coef.C = as.numeric(model.C$coefficients)
  print(x.C)
  print(coef.C)
  print(y.C)
  coef.R = as.numeric(model.R$coefficients)
  
  # fitted values (Lower & Upper)
  fitted_L <- (x.C %*% coef.C[-1]) - ((x.R %*% coef.R[-1]) / 2)
  fitted_U <- (x.C %*% coef.C[-1]) + ((x.R %*% coef.R[-1]) / 2)
  fitted_values <- cbind(fitted_L, fitted_U)
  colnames(fitted_values) <- c("fitted.Lower", "fitted.Upper")
  
  # residuals (Lower & Upper)
  residual_L <- (y.C - y.R / 2) - fitted_L
  residual_U <- (y.C + y.R / 2) - fitted_U
  residuals <- cbind(residual_L, residual_U)
  colnames(residuals) <- c("resid.Lower", "resid.Upper")
  
  result <- list(call = match.call(),
                 response = y,
                 predictor = x,
                 coefficients.Center = coef.C,
                 coefficeints.Range = coef.R,
                 fitted.values = fitted_values,
                 residuals = residuals)
  return(result)
}