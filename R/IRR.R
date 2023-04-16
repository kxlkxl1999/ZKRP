IRR <- function(formula = formula, data = data){
    
    # input formula for regression
    formula1 = as.formula(formula)
    model_interval = model.frame(formula = formula1, data = data)
    
    # response variable interval value
    yinterval = model.response(model_interval)
    y = as.matrix(yinterval)
    Y_c = (y[, 1] + y[, 2]) / 2   # center point of y
    Y_r = (y[, 2] - y[, 1]) / 2   # range point of y
    
    # predictor variable interval matrix (each left: Lower, each right: Upper)
    xinterval = model.matrix(attr(model_interval, "terms"), data = data)
    x = as.matrix(xinterval)
    
    p = (ncol(x) - 1)/2 # num. of predictor variable except intercept
    n = nrow(x)
    var.names = colnames(x)[(1:p)*2]
    
    X_c <- matrix(NA, n, p)   # center point of x
    colnames(X_c) <- var.names
    X_r <- matrix(NA, n, p)   # range point of x
    colnames(X_r) <- var.names
    for (j in 1:p) {
        for (i in 1:n) {
            X_c[i, j] = (x[i, 2*j] + x[i, 2*j+1]) / 2
            X_r[i, j] = (x[i, 2*j+1] - x[i, 2*j]) / 2
        }
    }
    
    X_2r = 2*X_r
    Y_2r = 2*Y_r
    lm.Bisquare.c <- rlm(Y_c ~ X_c,psi = psi.bisquare,maxit = 5000,test.vec = "coef") 
    lm.Bisquare.r <- rlm(Y_2r ~ X_2r,psi = psi.bisquare,maxit = 5000,test.vec = "coef") 
    
    summary_lm.Bisquare.c <- summary(lm.Bisquare.c)
    summary_lm.Bisquare.r <- summary(lm.Bisquare.r)
    coefficients.C <- summary_lm.Bisquare.c$coefficients[1:(p+1)]
    coefficients.R <- summary_lm.Bisquare.r$coefficients[1:(p+1)]
    coef.C = as.matrix(coefficients.C)
    coef.R = as.matrix(coefficients.R)
    
    result <- list(call = match.call(),
                   response = y,
                   predictor = x,
                   coefficients.Center = coef.C,
                   coefficeints.Range = coef.R
    )
    
    return(result)
}