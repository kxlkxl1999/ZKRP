SSLR <- function(formula = formula, data = data){
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
    X_2r = 2 * X_r
    Y_2r = 2 * Y_r
    Fisher <- function(data0){
        beta <- list()
        phi <- NULL
        n <- nrow(data0)
        result_SSLR <- lm(Y_c~., data = data0)
        beta[[1]] <- as.matrix(result_SSLR$coefficients,nrow=2)
        phi[1] <- sum(result_SSLR$residuals^2)/n
        u <- list()
        v <- list()
        D <- list()
        v0=2
        Y <- data0[,1]
        X <- as.matrix(cbind(rep(1,n),data0[,-1]))
        
        for (m in 1:1000) {
            r1 <- (Y - X%*%beta[[m]])^2
            u[[m]] <- r1/phi[m]
            v[[m]] <- (v0+1)/(v0+u[[m]])
            D[[m]] <- diag(as.numeric(v[[m]]))
            beta[[m+1]] <- solve(t(X)%*%D[[m]]%*%X)%*%t(X)%*%D[[m]]%*%Y 
            phi[m+1] <- t(r1)%*%v[[m]]/n
            if (dist(t(cbind(beta[[m+1]],beta[[m]]))) < 1e-4 && abs(phi[m+1]-phi[m] < 1e-4)){
                beta_final <- beta[[m+1]]
                phi_final <- phi[m+1]
                m_final <- m
                break}
        }
        
        result_final <- list("beta"=beta_final,"phi"=phi_final)
        
        return(result_final)
        
    }
    
    ## Gaussian distribution is applied for midpoint model
    data_c <- data.frame(Y_c,X_c)
    result_c <- lm(Y_c~., data = data_c)
    coefficients.C <- as.numeric(result_c$coefficients)
    
    #### Student-t distribution is applied for midpoint model 
    # data_c <- data.frame(Y_c,X_c)
    # result_c <- Fisher(data_c)
    # coefficients.C <- result_c$beta[1:(p+1)]
    
    ## Gaussian distribution is applied for range model
    data_r <- data.frame(Y_2r,X_2r)
    result_r <- lm(Y_2r~., data = data_r)
    coefficients.R <- as.numeric(result_r$coefficients)
    
    #### Student-t distribution is applied for range model 
    # data_r <- data.frame("Y_c"=Y_2r,X_2r)
    # result_r <- Fisher(data_r)
    # coefficients.R <- result_r$beta[1:(p+1)]
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