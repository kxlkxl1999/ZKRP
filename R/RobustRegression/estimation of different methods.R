
######  input: data, a list of length 4
# data[[1]] means the matrix of lower bounds of independent variables 
# data[[2]] means the matrix of upper bounds of independent variables 
# data[[3]] means the vector of lower bounds of dependent variable 
# data[[4]] means the vector of upper bounds of dependent variable 

######  output: parameters estimated 
 
############## p=1 for CRM
CRM_1 <- function(data){
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  X_c <- (X_lower + X_upper)/2
  X_2r <- X_upper - X_lower
  Y_c <- (Y_lower + Y_upper)/2
  Y_2r <- Y_upper - Y_lower
  
  data_crm <- data.frame(Y_c,X_c,Y_2r,X_2r)
  colnames(data_crm) <- c("Y_c","X1_c","Y_2r","X1_2r")
  
  res.crm <- crm("Y_c ~ X1_c", "Y_2r ~ X1_2r", data = data_crm)
  res.crm.summary <- summary(res.crm)
  coefficients.C <- res.crm.summary$coefficients.C[1:(p+1)]
  coefficients.R <- res.crm.summary$coefficients.R[1:(p+1)]
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}
############## p=3 for CRM
CRM_3 <- function(data){
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  X_c <- (X_lower + X_upper)/2
  X_2r <- X_upper - X_lower
  Y_c <- (Y_lower + Y_upper)/2
  Y_2r <- Y_upper - Y_lower
  
  data_crm <- data.frame(Y_c,X_c,Y_2r,X_2r)
  colnames(data_crm) <- c("Y_c","X1_c","X2_c","X3_c","Y_2r","X1_2r","X2_2r","X3_2r")
  
  res.crm <- crm("Y_c ~ X1_c + X2_c + X3_c", "Y_2r ~ X1_2r + X2_2r + X3_2r", data = data_crm)
  res.crm.summary <- summary(res.crm)
  coefficients.C <- res.crm.summary$coefficients.C[1:(p+1)]
  coefficients.R <- res.crm.summary$coefficients.R[1:(p+1)]
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}
############## p=6 for CRM
CRM_6 <- function(data){
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  X_c <- (X_lower + X_upper)/2
  X_2r <- X_upper - X_lower
  Y_c <- (Y_lower + Y_upper)/2
  Y_2r <- Y_upper - Y_lower
  
  data_crm <- data.frame(Y_c,X_c,Y_2r,X_2r)
  colnames(data_crm) <- c("Y_c","X1_c","X2_c","X3_c","X4_c","X5_c","X6_c","Y_2r","X1_2r","X2_2r","X3_2r","X4_2r","X5_2r","X6_2r")
  
  res.crm <- crm("Y_c ~ X1_c + X2_c + X3_c + X4_c + X5_c + X6_c", "Y_2r ~ X1_2r + X2_2r + X3_2r + X4_2r + X5_2r + X6_2r", data = data_crm)
  res.crm.summary <- summary(res.crm)
  coefficients.C <- res.crm.summary$coefficients.C[1:(p+1)]
  coefficients.R <- res.crm.summary$coefficients.R[1:(p+1)]
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}


############## p=1 for CCRM
CCRM_1 <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  X_c <- (X_lower + X_upper)/2
  X_2r <- X_upper - X_lower
  Y_c <- (Y_lower + Y_upper)/2
  Y_2r <- Y_upper - Y_lower
  
  data_ccrm <- data.frame(Y_c,X_c,Y_2r,X_2r)
  colnames(data_ccrm) <- c("Y_c","X1_c","Y_2r","X1_2r")
  
  res.ccrm <- ccrm("Y_c ~ X1_c", "Y_2r ~ X1_2r", data = data_ccrm)
  res.ccrm.summary <- summary(res.ccrm)
  coefficients.C <- res.ccrm.summary$coefficients.C[1:(p+1)]
  coefficients.R <- res.ccrm.summary$coefficients.R[1:(p+1)]
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}
############## p=3 for CCRM
CCRM_3 <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  X_c <- (X_lower + X_upper)/2
  X_2r <- X_upper - X_lower
  Y_c <- (Y_lower + Y_upper)/2
  Y_2r <- Y_upper - Y_lower
  
  data_ccrm <- data.frame(Y_c,X_c,Y_2r,X_2r)
  colnames(data_ccrm) <- c("Y_c","X1_c","X2_c","X3_c","Y_2r","X1_2r","X2_2r","X3_2r")
  
  res.ccrm <- ccrm("Y_c ~ X1_c + X2_c + X3_c", "Y_2r ~ X1_2r + X2_2r + X3_2r", data = data_ccrm)
  res.ccrm.summary <- summary(res.ccrm)
  coefficients.C <- res.ccrm.summary$coefficients.C[1:(p+1)]
  coefficients.R <- res.ccrm.summary$coefficients.R[1:(p+1)]
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}
############## p=6 for CCRM
CCRM_6 <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  X_c <- (X_lower + X_upper)/2
  X_2r <- X_upper - X_lower
  Y_c <- (Y_lower + Y_upper)/2
  Y_2r <- Y_upper - Y_lower
  
  data_ccrm <- data.frame(Y_c,X_c,Y_2r,X_2r)
  colnames(data_ccrm) <- c("Y_c","X1_c","X2_c","X3_c","X4_c","X5_c","X6_c","Y_2r","X1_2r","X2_2r","X3_2r","X4_2r","X5_2r","X6_2r")
  
  res.ccrm <- ccrm("Y_c ~ X1_c + X2_c + X3_c + X4_c + X5_c + X6_c", "Y_2r ~ X1_2r + X2_2r + X3_2r + X4_2r + X5_2r + X6_2r", data = data_ccrm)
  res.ccrm.summary <- summary(res.ccrm)
  coefficients.C <- res.ccrm.summary$coefficients.C[1:(p+1)]
  coefficients.R <- res.ccrm.summary$coefficients.R[1:(p+1)]
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}


############## CCRJM
CCRJM <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  n <- nrow(X_lower)
  p <- ncol(X_lower)
  X_c <- (X_lower + X_upper) / 2
  X_r <- (X_upper - X_lower) / 2
  Y_c <- (Y_lower + Y_upper) / 2
  Y_r <- (Y_upper - Y_lower) / 2
  
  
  Y <- c(Y_c,Y_r)
  X <- cbind(1,X_c,X_r)
  Z <- rbind(cbind(X,matrix(0,n,ncol(X))),cbind(matrix(0,n,ncol(X)),X))
  G <- rbind(cbind(-X,-X),cbind(X,-X),cbind(matrix(0,n,ncol(X)),-X))
  h <- c(Y_r - Y_c , Y_r + Y_c, rep(0,n))
  beta <- rep(1,2*(2*p+1))
  
  # objective function
  fn = function(beta){
    sum((Y-Z%*%beta)^2)
  }
  
  # solution
  ret = donlp2(par=beta, fn,
               par.upper = rep(+Inf, length(beta)),
               par.lower = rep(-Inf, length(beta)),
               A=G,
               lin.lower = rep(-Inf,3*(n)),
               lin.upper = h)
  
  # output
  beta_fit_CCRJM <- ret$par
  coefficients.C <- beta_fit_CCRJM[1:(2*p+1)]
  coefficients.R <- beta_fit_CCRJM[(2*p+2):(2*(2*p+1))]
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}


############## ID
ID <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  n <- nrow(X_lower)
  p <- ncol(X_lower)
  X_c <- (X_lower + X_upper) / 2
  X_r <- (X_upper - X_lower) / 2
  Y_c <- (Y_lower + Y_upper) / 2
  Y_r <- (Y_upper - Y_lower) / 2
  
  
  Y <- c(Y_c,Y_r/(sqrt(3)))
  X.c <- NULL
  X.r <- NULL
  for (j in 1:(2*p)) {
    X.c <- cbind(X.c,-X_c[,ceiling(j/2)]*((-1)^(j)))
    X.r <- cbind(X.r,X_r[,ceiling(j/2)])
  }
  
  X.c <- cbind(X.c,1)
  X.r <- cbind(X.r,0)
  X <- rbind(X.c,X.r/(sqrt(3)))
  G <- diag( c(rep(1,2*p),0), (2*p+1), (2*p+1))
  
  beta <- rep(0.5,(2*p+1))
  
  fn = function(beta){
    sum((Y-X%*%beta)^2)
  }
  
  ret = donlp2(par=beta, fn, par.upper = rep(+Inf, (2*p+1)), par.lower = c(rep(0, 2*p),-Inf))
  
  beta_fit_IR <- ret$par
  a <- beta_fit_IR[seq(1,2*p,by=2)]
  b <- beta_fit_IR[seq(2,2*p,by=2)]
  v <- beta_fit_IR[2*p+1]
  coefficients.C <- c(v,a-b)
  coefficients.R <- a+b
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}


############## SSLR
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


############## IRR
IRR <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  n <- nrow(X_lower)
  p <- ncol(X_lower)
  X_c <- (X_lower + X_upper) / 2
  X_2r <- (X_upper - X_lower) 
  Y_c <- (Y_lower + Y_upper) / 2
  Y_2r <- (Y_upper - Y_lower) 
  
  lm.Bisquare.c <- rlm(Y_c ~ X_c,psi = psi.bisquare,maxit = 5000,test.vec = "coef") 
  lm.Bisquare.r <- rlm(Y_2r ~ X_2r,psi = psi.bisquare,maxit = 5000,test.vec = "coef") 
  
  summary_lm.Bisquare.c <- summary(lm.Bisquare.c)
  summary_lm.Bisquare.r <- summary(lm.Bisquare.r)
  coefficients.C <- summary_lm.Bisquare.c$coefficients[1:(p+1)]
  coefficients.R <- summary_lm.Bisquare.r$coefficients[1:(p+1)]
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}


############## IQR
IQR <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  n <- nrow(X_lower)
  p <- ncol(X_lower)
  X_c <- (X_lower + X_upper) / 2
  X_r <- (X_upper - X_lower) / 2
  Y_c <- (Y_lower + Y_upper) / 2
  Y_r <- (Y_upper - Y_lower) / 2
  
  data_c <- data.frame(Y_c,X_c)
  result_c <- rq(Y_c~.,tau=.5,data = data_c)
  coefficients.C <-  as.numeric(result_c$coefficients)
  
  data_r <- data.frame(Y_r,X_r)
  result_r <- rq(Y_r~.,tau=.5, data = data_r)
  coefficients.R <- as.numeric(result_r$coefficients)
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}


############## different hyper-parameter estimator for iETKRR such as

## hyper-parameter estimators (S1,S1) for iETKRR method in midpoints and in ranges, respectively
iETKRR11 <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  n <- nrow(X_lower)
  p <- ncol(X_lower)
  X_c <- (X_lower + X_upper) / 2
  X_r <- (X_upper - X_lower) / 2
  Y_c <- (Y_lower + Y_upper) / 2
  Y_r <- (Y_upper - Y_lower) / 2
  
  mod.k.pm = kernel.reg(X_c, Y_c, 1e-10, 100)
  mod.k.h = kernel.reg(X_r, Y_r, 1e-10, 100)
  
  coefficients.C <- mod.k.pm$coef
  coefficients.R <- mod.k.h$coef
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}

## hyper-parameter estimators (S3,S1) for iETKRR method in midpoints and in ranges, respectively
iETKRR31 <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  n <- nrow(X_lower)
  p <- ncol(X_lower)
  X_c <- (X_lower + X_upper) / 2
  X_r <- (X_upper - X_lower) / 2
  Y_c <- (Y_lower + Y_upper) / 2
  Y_r <- (Y_upper - Y_lower) / 2
  
  mod.k3.pm = kernel.reg3(X_c, Y_c, 1e-10, 100)
  mod.k.h = kernel.reg(X_r, Y_r, 1e-10, 100)
  
  coefficients.C <- mod.k3.pm$coef
  coefficients.R <- mod.k.h$coef
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}

## hyper-parameter estimators (S1,S3) for iETKRR method in midpoints and in ranges, respectively
iETKRR13 <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  n <- nrow(X_lower)
  p <- ncol(X_lower)
  X_c <- (X_lower + X_upper) / 2
  X_r <- (X_upper - X_lower) / 2
  Y_c <- (Y_lower + Y_upper) / 2
  Y_r <- (Y_upper - Y_lower) / 2
  
  mod.k.pm = kernel.reg(X_c, Y_c, 1e-10, 100)
  mod.k3.h = kernel.reg3(X_r, Y_r, 1e-10, 100)
  
  coefficients.C <- mod.k.pm$coef
  coefficients.R <- mod.k3.h$coef
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}

## hyper-parameter estimators (S1,S4) for iETKRR method in midpoints and in ranges, respectively
iETKRR14 <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  n <- nrow(X_lower)
  p <- ncol(X_lower)
  X_c <- (X_lower + X_upper) / 2
  X_r <- (X_upper - X_lower) / 2
  Y_c <- (Y_lower + Y_upper) / 2
  Y_r <- (Y_upper - Y_lower) / 2
  
  mod.k.pm = kernel.reg(X_c, Y_c, 1e-10, 100)
  mod.k4.h = kernel.reg4(X_r, Y_r, 1e-10, 100)
  
  coefficients.C <- mod.k.pm$coef
  coefficients.R <- mod.k4.h$coef
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}

## hyper-parameter estimators (S3,S4) for iETKRR method in midpoints and in ranges, respectively
iETKRR34 <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  n <- nrow(X_lower)
  p <- ncol(X_lower)
  X_c <- (X_lower + X_upper) / 2
  X_r <- (X_upper - X_lower) / 2
  Y_c <- (Y_lower + Y_upper) / 2
  Y_r <- (Y_upper - Y_lower) / 2
  
  mod.k3.pm = kernel.reg3(X_c, Y_c, 1e-10, 100)
  mod.k4.h = kernel.reg4(X_r, Y_r, 1e-10, 100)
  
  coefficients.C <- mod.k3.pm$coef
  coefficients.R <- mod.k4.h$coef
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}


############## Fuzzy-IR
Fuzzy <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  n <- nrow(X_lower)
  p <- ncol(X_lower)
  print(p)
  X_c <- (X_lower + X_upper) / 2
  X_r <- (X_upper - X_lower) / 2
  Y_c <- (Y_lower + Y_upper) / 2
  Y_r <- (Y_upper - Y_lower) / 2
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
  print(length(coefficients.C))
  coefficients.R <- coefficients[(2*p+2):(4*p+2)]
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}


############## LN-IRR
LN_IRR <- function(data){
  
  X_lower <- as.matrix(data[[1]])
  X_upper <- as.matrix(data[[2]])
  Y_lower <- data[[3]]
  Y_upper <- data[[4]]
  n <- nrow(X_lower)
  p <- ncol(X_lower)
  X_c <- (X_lower + X_upper) / 2
  X_r <- (X_upper - X_lower) / 2
  Y_c <- (Y_lower + Y_upper) / 2
  Y_r <- (Y_upper - Y_lower) / 2
  Y_r_ln <- log(Y_r)
  
  data_c <- data.frame(Y_c,X_c,X_r)
  data_r <- data.frame(Y_r_ln,X_c,X_r)
  
  lm.huber.c <- rlm(Y_c ~ .,data = data_c, psi = psi.huber,maxit=5000,test.vec = "coef")
  lm.huber.r <- rlm(Y_r_ln ~ .,data = data_r,psi = psi.huber,maxit=5000,test.vec = "coef") 
  
  summary_lm.huber.c <- summary(lm.huber.c)
  summary_lm.huber.r <- summary(lm.huber.r)
  coefficients.C <- summary_lm.huber.c$coefficients[1:(2*p+1)]
  coefficients.R <- summary_lm.huber.r$coefficients[1:(2*p+1)]
  
  result <- list("coefficients.C"=coefficients.C,"coefficients.R"=coefficients.R)
  return(result)
}






