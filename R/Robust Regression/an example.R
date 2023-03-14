
library(iRegression)
library(MASS)
library(quantreg)
library(Rdonlp2)
library(sm)
library(Matrix)

source("synthetic data generation.R")
source("Hyper-parameter estimators functions of iETKRR.R")
source("estimation of different methods.R")
source("prediction of different methods.R")

# n means sample size
# p means the number of independent variables
# out.p means the percentage of outliers
# err means the error variability with Unif(0,err)
# out.index means the existence of Y^c outliers when out.index=1 while Y^r outliers when out.index=2

n=90  # n=90,300,900 
p=1   # p=1,3   
out.p=0.05  #out.p=0.05,0.1,0.15
err=20  #err=20,40
out.index=1  # out.index=1,2 

iteration=100 ## number of repeated Monte Carlo simulations
accuracy <- list()  

## choose appropriate distribution assumption of midpoint and range model for SSLR method
SSLR <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  n <- nrow(X_lower)
  p <- ncol(X_lower)
  X_c <- (X_lower + X_upper)/2
  X_2r <- (X_upper - X_lower)
  Y_c <- (Y_lower + Y_upper)/2
  Y_2r <- (Y_upper - Y_lower)
  
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
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
  
}
## choose right functions of CRM and CCRM corresponding to p
CRM_pre <- CRM_pre_1
CCRM_pre <- CCRM_pre_1
## choose best hyper-parameter estimator for iETKRR method
iETKRR_pre <- iETKRR11_pre

source("method compare.R")

for (l in 1:iteration) {
  dataset <- data_generate_II(n,p,out.p,err,out.index)
  accuracy[[l]] <- method_compare(dataset)
}

## calculate the empirical average and standard deviation of each measurement of each method
mes <- NULL
accuracy_mean <- matrix(0,10,4)
accuracy_sd <- matrix(0,10,4)
for (i in 1:10) {
  for (j in 1:4) {
    for (l in 1:iteration) {
      mes[l] <- accuracy[[l]][i,j]
    }
    accuracy_mean[i,j] <- mean(mes)
    accuracy_sd[i,j] <- sd(mes)
  }
}
rownames(accuracy_mean) <- c("CRM","CCRM","CCRJM","ID","SSLR","IRR","IQR","iETKRR","Fuzzy-RR","LN-IRR")
colnames(accuracy_mean) <- c("RMSE_L","RMSE_U","RMSE_H","AR")
rownames(accuracy_sd) <- c("CRM","CCRM","CCRJM","ID","SSLR","IRR","IQR","iETKRR","Fuzzy-RR","LN-IRR")
colnames(accuracy_sd) <- c("RMSE_L","RMSE_U","RMSE_H","AR")
## the empirical average of each measurement of each method
accuracy_mean
## the empirical standard deviation of each measurement of each method
accuracy_sd


### paired t-test between each competitor and proposed method
index1 <-NULL
index2 <-NULL
result_CI <- data.frame()
for (k in 1:4) {
  for (j in 1:9) {
    for (l in 1:iteration) {
      index1[l] <- accuracy[[l]][j,k]
      index2[l] <- accuracy[[l]][10,k]
    }
    result.t <- t.test(index1,index2,alternative = "two.sided",paired = T)
    conf <- as.vector(result.t$conf.int)
    conf <- round(conf,4)
    result_CI[k,j] <- paste("[",conf[1],",",conf[2],"]")
  }
}
colnames(result_CI) <- c("CRM","CCRM","CCRJM","ID","SSLR","IRR","IQR","iETKRR","Fuzzy")
rownames(result_CI) <- c("RMSE_L","RMSE_U","RMSE_H","AR")
### 95% confidence intervals of paired t-test between each competitor and proposed method
result_CI
