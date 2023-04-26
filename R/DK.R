dk<- function(beta, xc, xr, yc, yr){
    n = nrow(xc)
    p = length(beta)
    beta_c = beta[1:(p/2)]
    beta_r = beta[(p/2+1):p]
    yc_hat = cbind(1,xc) %*% beta_c
    yr_hat = cbind(1,xr) %*% beta_r
    yl = yc-yr
    yu = yc+yr
    yl_hat = yc_hat - yr_hat
    yu_hat = yc_hat + yr_hat
    
    return(mean((yu_hat-yu)^2/2+(yl_hat-yl)^2/2+2*(yc_hat-yc)^2))
}

Dk <- function(formula = formula, data = data){
    # model = rq(y~x, method = "pfn")
    # return(as.numeric(model$coefficients))
    
    # input formula for regression
    formula1 = as.formula(formula)
    model_interval = model.frame(formula = formula1, data = data)
    
    # response variable interval value
    yinterval = model.response(model_interval)
    y = as.matrix(yinterval)
    y.C = (y[, 1] + y[, 2]) / 2   # center point of y
    y.R = (y[, 2] - y[, 1]) / 2   # range point of y
    
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
            x.R[i, j] = (x[i, 2*j+1] - x[i, 2*j]) / 2
        }
    }
    hf1 = HF1(formula = formula, data = data)
    beta0 = c(as.vector(hf1$coefficients.Center), as.vector(hf1$coefficeints.Range))
    coef.C=c()
    coef.R=c()
    
    beta_dk = optim(par=beta0, fn=dk, xc=x.C, xr=x.R, yc=y.C, yr=y.R)
    coef.C = beta_dk$par[1:(p/2+1)]
    coef.R = beta_dk$par[(p/2+2):(p+2)]
    
    result <- list(call = match.call(),
                   response = y,
                   predictor = x,
                   coefficients.Center = coef.C,
                   coefficeints.Range = coef.R)
    return(result)
}