Fuzzy <- function(formula = formula, data = data){
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
    p = p / 2
    Y_r_ln <- log(Y_r)
    
    X_all <- cbind(X_c, X_r)
    Y <- c(Y_c, Y_r_ln)
    block <- list()
    block[[1]] <- cbind(X_all,1)
    block[[2]] <- cbind(X_all,1)
    X <- as.matrix(bdiag(block))
    
    fuzzy_result <- rlm(Y ~ X-1, method =  "MM", psi = psi.bsiquare, maxit=5000) 
    summary_fuzzy_result <- summary(fuzzy_result)
    coefficients <- summary_fuzzy_result$coefficients
    
    coefficients.C <- coefficients[1:(2*p+1)]
    coefficients.R <- coefficients[(2*p+2):(4*p+2)]
    
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