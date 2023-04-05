source('CM.R')
source('CRM.R')
source('CCRM.R')
source('HF1.R')
source('IRR.R')
source('SSLR.R')
source('Fuzzy.R')
source('LN_IRR.R')
source('HF_family.R')

rmse <- function(yc, yr, yhatc, yhatr){
    yl = yc - yr
    yu = yc + yr
    yhatl = yhatc - yhatr
    yhatu = yhatc + yhatr
    return(c(
        sqrt(mean((yc-yhatc)^2)),
        sqrt(mean((yr-yhatr)^2)),
        sqrt(mean((yl-yhatl)^2)),
        sqrt(mean((yu-yhatu)^2))
    ))
}

eval <- function(model, testdata){
    coef.c = as.vector(model$coefficients.Center)
    coef.r = as.vector(model$coefficeints.Range)
    testdata = as.matrix(testdata)
    m = ncol(testdata)
    nx = (m-2)/2
    yl = testdata[,1]
    yu = testdata[,2]
    xl = testdata[,3:(2+nx)]
    xu = testdata[,(3+nx):m]
    xc = (xu + xl) / 2
    xr = (xu - xl) / 2
    yc = (yu + yl) / 2
    yr = (yu - yl) / 2
    yhatc = cbind(1,xc) %*% coef.c
    yhatr = cbind(1,xr) %*% coef.r
    return(rmse(yc, yr, yhatc, yhatr))
}

eval_fuzzy <- function(model, testdata){
    
    coef.c = as.vector(model$coefficients.Center)
    coef.r = as.vector(model$coefficeints.Range)
    testdata = as.matrix(testdata)
    m = ncol(testdata)
    nx = (m-2)/2
    yl = testdata[,1]
    yu = testdata[,2]
    xl = testdata[,3:(2+nx)]
    xu = testdata[,(3+nx):m]
    xc = (xu + xl) / 2
    xr = (xu - xl) / 2
    yc = (yu + yl) / 2
    yr = (yu - yl) / 2
    X_test <-cbind(xc,xr,1)
    
    yhatc <- as.vector(X_test%*%coef.c)
    Y_r_ln_fit <- as.vector(X_test%*%coef.r)
    yhatr <- exp(Y_r_ln_fit)
    return(rmse(yc, yr, yhatc, yhatr))
}

eval_IRR <- function(model, testdata){
    coef.c = as.vector(model$coefficients.Center)
    coef.r = as.vector(model$coefficeints.Range)
    testdata = as.matrix(testdata)
    m = ncol(testdata)
    nx = (m-2)/2
    yl = testdata[,1]
    yu = testdata[,2]
    xl = testdata[,3:(2+nx)]
    xu = testdata[,(3+nx):m]
    xc = (xu + xl) / 2
    xr = (xu - xl) / 2
    yc = (yu + yl) / 2
    yr = (yu - yl) / 2
    x2r = 2*xr
    # print("111")
    # print(ncol(cbind(1,xc)))
    # print("222")
    # print(length(coef.c))
    
    yhatc <- as.vector(cbind(1,xc)%*%coef.c)
    Y_2r_fit <- as.vector(cbind(1,x2r)%*%coef.r)
    yhatr <- Y_2r_fit/2
    return(rmse(yc, yr, yhatc, yhatr))
}

eval_LN_IRR <- function(model, testdata){
    coef.c = as.vector(model$coefficients.Center)
    coef.r = as.vector(model$coefficeints.Range)
    testdata = as.matrix(testdata)
    m = ncol(testdata)
    nx = (m-2)/2
    yl = testdata[,1]
    yu = testdata[,2]
    xl = testdata[,3:(2+nx)]
    xu = testdata[,(3+nx):m]
    xc = (xu + xl) / 2
    xr = (xu - xl) / 2
    yc = (yu + yl) / 2
    yr = (yu - yl) / 2
    X_test <-cbind(xc,xr,1)
    
    yhatc <- as.vector(X_test%*%coef.c)
    Y_r_ln_fit <- as.vector(X_test%*%coef.r)
    yhatr <- exp(Y_r_ln_fit)
    return(rmse(yc, yr, yhatc, yhatr))
}

eval_SSLR <- function(model, testdata){
    coef.c = as.vector(model$coefficients.Center)
    coef.r = as.vector(model$coefficeints.Range)
    testdata = as.matrix(testdata)
    m = ncol(testdata)
    nx = (m-2)/2
    yl = testdata[,1]
    yu = testdata[,2]
    xl = testdata[,3:(2+nx)]
    xu = testdata[,(3+nx):m]
    xc = (xu + xl) / 2
    xr = (xu - xl) / 2
    yc = (yu + yl) / 2
    yr = (yu - yl) / 2
    x2r = 2 * xr
    yhatc <- as.vector(cbind(1,xc)%*%coef.c)
    Y_2r_fit <- as.vector((cbind(1,x2r)%*%coef.r))
    yhatr <- Y_2r_fit/2
    return(rmse(yc, yr, yhatc, yhatr))
}

comparasion <- function(formula, dataset)
{
    k=10
    set.seed(0)
    folds = createFolds(dataset[,1],k=k)
    rmse_cm = matrix(0,nrow = 4)
    rmse_crm = matrix(0,nrow = 4)
    rmse_ccrm = matrix(0,nrow = 4)
    rmse_irr = matrix(0,nrow = 4)
    rmse_fuzzy = matrix(0,nrow = 4)
    rmse_sslr = matrix(0,nrow = 4)
    rmse_ln_irr = matrix(0,nrow = 4)
    rmse_hf1 = matrix(0,nrow = 4)
    for(i in 1:k)
    {
        train = dataset[-folds[[i]], ]
        test = dataset[folds[[i]], ]
        
        # print("train")
        # print(ncol(train))
        # print("test")
        # print(ncol(test))
        
        model_cm = CM(formula, train)
        model_crm = CRM(formula, train)
        model_ccrm = CCRM(formula, train)
        model_hf1 = HF1(formula, train)
        model_irr = IRR(formula, train)
        model_fuzzy = Fuzzy(formula, train)
        model_sslr = SSLR(formula, train)
        model_ln_irr = LN_IRR(formula, train)
        rmse_cm = cbind(rmse_cm, eval(model_cm, test))
        rmse_crm = cbind(rmse_crm, eval(model_crm, test))
        rmse_ccrm = cbind(rmse_ccrm, eval(model_ccrm, test))
        rmse_hf1 = cbind(rmse_hf1, eval(model_hf1, test))
        rmse_irr = cbind(rmse_cm, eval_IRR(model_irr, test))
        rmse_fuzzy = cbind(rmse_crm, eval_fuzzy(model_fuzzy, test))
        rmse_sslr = cbind(rmse_ccrm, eval_SSLR(model_sslr, test))
        rmse_ln_irr = cbind(rmse_hf1, eval_LN_IRR(model_ln_irr, test))
    }
    
    return(cbind(
        apply(rmse_cm[,-1], 1, mean),
        apply(rmse_crm[,-1], 1, mean),
        apply(rmse_ccrm[,-1], 1, mean),
        apply(rmse_hf1[,-1], 1, mean),
        apply(rmse_irr[,-1], 1, mean),
        apply(rmse_fuzzy[,-1], 1, mean),
        apply(rmse_sslr[,-1], 1, mean),
        apply(rmse_ln_irr[,-1], 1, mean)
        )
    )
}

outlier_comparasion <- function(formula, dataset, seed, outlierType=1, alpha=0.1,k = 10)
{
    set.seed(0)
    folds = createFolds(dataset[,1],k=k)
    rmse_cm = matrix(0,nrow = 4)
    rmse_crm = matrix(0,nrow = 4)
    rmse_ccrm = matrix(0,nrow = 4)
    rmse_irr = matrix(0,nrow = 4)
    rmse_fuzzy = matrix(0,nrow = 4)
    rmse_sslr = matrix(0,nrow = 4)
    rmse_ln_irr = matrix(0,nrow = 4)
    rmse_hf1 = matrix(0,nrow = 4)
    for(i in 1:k)
    {
        train = dataset[-folds[[i]], ]
        test = dataset[folds[[i]], ]
        train = data_generation_outlier_realdata(train, 0, outlierType = outlierType, alpha = alpha)
        
        model_cm = CM(formula, train)
        model_crm = CRM(formula, train)
        model_ccrm = CCRM(formula, train)
        model_hf1 = HF1(formula, train)
        model_irr = IRR(formula, train)
        model_fuzzy = Fuzzy(formula, train)
        model_sslr = SSLR(formula, train)
        model_ln_irr = LN_IRR(formula, train)
        rmse_cm = cbind(rmse_cm, eval(model_cm, test))
        rmse_crm = cbind(rmse_crm, eval(model_crm, test))
        rmse_ccrm = cbind(rmse_ccrm, eval(model_ccrm, test))
        rmse_hf1 = cbind(rmse_hf1, eval(model_hf1, test))
        rmse_irr = cbind(rmse_cm, eval_IRR(model_irr, test))
        rmse_fuzzy = cbind(rmse_crm, eval_fuzzy(model_fuzzy, test))
        rmse_sslr = cbind(rmse_ccrm, eval_SSLR(model_sslr, test))
        rmse_ln_irr = cbind(rmse_hf1, eval_LN_IRR(model_ln_irr, test))
    }
    
    result1 = cbind(
            apply(rmse_cm[,-1], 1, mean),
            apply(rmse_crm[,-1], 1, mean),
            apply(rmse_ccrm[,-1], 1, mean),
            apply(rmse_hf1[,-1], 1, mean),
            apply(rmse_irr[,-1], 1, mean),
            apply(rmse_fuzzy[,-1], 1, mean),
            apply(rmse_sslr[,-1], 1, mean),
            apply(rmse_ln_irr[,-1], 1, mean))
    result2 = cbind(
            apply(rmse_cm[,-1], 1, sd),
            apply(rmse_crm[,-1], 1, sd),
            apply(rmse_ccrm[,-1], 1, sd),
            apply(rmse_hf1[,-1], 1, sd),
            apply(rmse_irr[,-1], 1, sd),
            apply(rmse_fuzzy[,-1], 1, sd),
            apply(rmse_sslr[,-1], 1, sd),
            apply(rmse_ln_irr[,-1], 1, sd))
    colnames(result1) = c("cm","crm",'ccrm','hf1','irr','fuzzy','sslr','ln_irr')
    rownames(result1) = c('rmsec','rmser','rmsel','rmseu')
    colnames(result2) = c("cm","crm",'ccrm','hf1','irr','fuzzy','sslr','ln_irr')
    rownames(result2) = c('rmsec','rmser','rmsel','rmseu')
    
    return(list(result1,result2)
    )
}
