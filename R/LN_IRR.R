LN_IRR <- function(formula = formula, data = data){
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
    
    p = ncol(x) - 1 # num. of predictor variable except intercept
    n = nrow(x)
    var.names = colnames(x)[(1:{p/2})*2]
    
    X_c <- matrix(NA, n, {p/2})   # center point of x
    colnames(X_c) <- var.names
    X_r <- matrix(NA, n, {p/2})   # range point of x
    colnames(X_r) <- var.names
    for (j in 1:{p/2}) {
        for (i in 1:n) {
            X_c[i, j] = (x[i, 2*j] + x[i, 2*j+1]) / 2
            X_r[i, j] = (x[i, 2*j+1] - x[i, 2*j]) / 2
        }
    }
    Y_r_ln <- log(Y_r)
    p = p / 2
    
    data_c <- data.frame(Y_c,X_c,X_r)
    data_r <- data.frame(Y_r_ln,X_c,X_r)
    
    lm.huber.c <- rlm(Y_c ~ .,data = data_c, psi = psi.huber,maxit=5000,test.vec = "coef")
    lm.huber.r <- rlm(Y_r_ln ~ .,data = data_r,psi = psi.huber,maxit=5000,test.vec = "coef") 
    
    summary_lm.huber.c <- summary(lm.huber.c)
    summary_lm.huber.r <- summary(lm.huber.r)
    coefficients.C <- summary_lm.huber.c$coefficients[1:(2*p+1)]
    coefficients.R <- summary_lm.huber.r$coefficients[1:(2*p+1)]
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