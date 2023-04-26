hf <- function(beta, xc, xr, yc, yr){
    p = length(beta)
    beta_c = beta[1:(p/2)]
    beta_r = beta[(p/2+1):p]
    yc_hat = cbind(1,xc) %*% beta_c
    yr_hat = cbind(1,xr) %*% beta_r
    return(mean(abs(yc-yc_hat) + abs(yr-yr_hat)))
}

hf_normalization1 <- function(beta, xc, xr, yc, yr){
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
    data = cbind(yl,yu,yl_hat,yu_hat)
    interval_union <- function(x) {
        a=x[1]
        b=x[2]
        c=x[3]
        d=x[4]
        return(min(max(b,d)-min(a,c), b-a+d-c))
        }
    u = apply(data, 1, interval_union)
    return(mean((abs(yc-yc_hat) + abs(yr-yr_hat))/u))
}

hf_normalization2 <- function(beta, xc, xr, yc, yr){
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
    return(mean((abs(yc-yc_hat) + abs(yr-yr_hat))/(yl-yu+yl_hat-yu_hat)))
}

hf_sum <- function(beta, xc, xr, yc, yr){
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
    return(mean(abs(yl-yl_hat)+abs(yu-yu_hat)))
    # return(mean(abs(yc-yc_hat) + abs(yr-yr_hat)))
}

hf_sum1 <- function(beta, xc, xr, yc, yr){
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
    
    data = cbind(yl,yu,yl_hat,yu_hat)
    interval_union <- function(x) {
        a=x[1]
        b=x[2]
        c=x[3]
        d=x[4]
        return(min(max(b,d)-min(a,c), b-a+d-c))
    }
    u = apply(data, 1, interval_union)
    return(mean((abs(yl-yl_hat)+abs(yu-yu_hat))/u))
}

hf_sum2 <- function(beta, xc, xr, yc, yr){
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
    return(mean((abs(yl-yl_hat)+abs(yu-yu_hat))/(yl-yu+yl_hat-yu_hat)))
}

hf_modified<- function(beta, xc, xr, yc, yr){
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
    data = cbind(yl,yu,yl_hat,yu_hat)
    hf_m <- function(x) {
        x1=x[1]
        x2=x[2]
        x3=x[3]
        x4=x[4]
        l_union = min(max(x2, x4)-min(x1, x3), x2-x1+x4-x3)
        l1 = x2-x1
        l2 = x4-x3
        if(l_union == l1 + l2)
            {return((abs(x1-x3)+abs(x3-x2))/2 + (abs(x2-x4)+abs(x3-x2))/2)}
        else if(l_union>l1 & l_union>l2)
            {return(max((x1-x3)^2,(x2-x4)^2)/max(2*l1, 2*l2) + min((x1-x3)^2,(x2-x4)^2)/min(2*l1, 2*l2))}
        else
            {return(((x1-x3)^2+(x2-x4)^2)/(2*max(l1, l2)))}
    }
    u = apply(data, 1, hf_m)
    return(mean(u))
}

HF_family <- function(formula = formula, data = data, method='normalization1'){
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
    if(method=='normalization1')
    {
        hf = optim(par=beta0, fn=hf_normalization1, xc=x.C, xr=x.R, yc=y.C, yr=y.R)
        coef.C = hf$par[1:(p/2+1)]
        coef.R = hf$par[(p/2+2):(p+2)]
    }
    else if(method=='normalization2')
    {
        hf = optim(par=beta0, fn=hf_normalization2, xc=x.C, xr=x.R, yc=y.C, yr=y.R)
        coef.C = hf$par[1:(p/2+1)]
        coef.R = hf$par[(p/2+2):(p+2)]
    }
    else if(method=='sum')
    {
        hf = optim(par=beta0, fn=hf_sum, xc=x.C, xr=x.R, yc=y.C, yr=y.R)
        coef.C = hf$par[1:(p/2+1)]
        coef.R = hf$par[(p/2+2):(p+2)]
    }
    else if(method=='sum1')
    {
        hf = optim(par=beta0, fn=hf_sum1, xc=x.C, xr=x.R, yc=y.C, yr=y.R)
        coef.C = hf$par[1:(p/2+1)]
        coef.R = hf$par[(p/2+2):(p+2)]
    }
    else if(method=='sum2')
    {
        hf = optim(par=beta0, fn=hf_sum2, xc=x.C, xr=x.R, yc=y.C, yr=y.R)
        coef.C = hf$par[1:(p/2+1)]
        coef.R = hf$par[(p/2+2):(p+2)]
    }
    else if(method=='modified')
    {
        hf = optim(par=beta0, fn=hf_modified, xc=x.C, xr=x.R, yc=y.C, yr=y.R)
        coef.C = hf$par[1:(p/2+1)]
        coef.R = hf$par[(p/2+2):(p+2)]
    }
    else if(method=='hf')
    {
        hf = optim(par=beta0, fn=hf, xc=x.C, xr=x.R, yc=y.C, yr=y.R)
        coef.C = hf$par[1:(p/2+1)]
        coef.R = hf$par[(p/2+2):(p+2)]
    }
    
    result <- list(call = match.call(),
                   response = y,
                   predictor = x,
                   coefficients.Center = coef.C,
                   coefficeints.Range = coef.R)
    return(result)
}