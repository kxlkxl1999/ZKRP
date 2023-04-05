
# n means sample size
# p means the number of independent variables
# out.p means the percentage of outliers
# err means the error variability with Unif(0,err)
# out.index means the existence of Y^c outliers when out.index=1 while Y^r outliers when out.index=2

## data Configuration I where Y^c is merely related to X^c, and Y^r is merely related to X^r.
data_generate_I <- function(n,p,out.p,err,out.index){
  
  X_lower <- matrix(0, n, p)
  Y_lower <- NULL
  X_upper <- matrix(0, n, p)
  Y_upper <- NULL
  X_c <- matrix(0, n, p)
  X_r <- matrix(0, n, p)
  e_c <- NULL
  e_r <- NULL
  Y_c <- NULL
  Y_r <- NULL
  
  for (k in 1:n) {
    
    for (j in 1:p) {
      X_c[k, ] <- runif(p, 0, 80)
      X_r[k, ] <- runif(p, 5, 10)
    }
  
    ### parameters 
    beta_c <- runif(p+1, 2, 2.5)
    beta_r <- runif(p+1, 3, 3.5)
    
    ### errors
    e_c[k] <- runif(1,-err,err)
    e_r[k] <- runif(1,0, err)
    
    Y_c[k] <- c(1, X_c[k, ]) %*% beta_c + e_c[k]
    Y_r[k] <- c(1, X_r[k, ]) %*% beta_r + e_r[k]
  }
  
  X_lower <- X_c - X_r
  X_upper <- X_c + X_r
  Y_lower <- Y_c - Y_r
  Y_upper <- Y_c + Y_r
  
  ##### the synthetic data set is divided into a training set and a test set randomly with two thirds 
  ##### and one third of the samples, respectively.
  train_ind <- sample(seq_len(n), size = n*2/3 )
  
  data_test <- list(X_lower[-train_ind,],X_upper[-train_ind,],Y_lower[-train_ind],Y_upper[-train_ind])
  
  S_Y <- sd(Y_c[train_ind])
  data <- cbind(Y_c,Y_r,X_c,X_r)[train_ind,]
  data <- data[order(data[,1],decreasing=F),]

  ##### generate outliers
  if (out.p>0) {
    n0 <- n *2/3 * out.p;
    for (i in  1:n0) {
      data[i,out.index] <- data[i,out.index] + 9*S_Y
    }
  }
  
  X_lower_tr <- data[,3:(3+p-1)] - data[,(3+p):(3+2*p-1)]
  X_upper_tr <- data[,3:(3+p-1)] + data[,(3+p):(3+2*p-1)]
  Y_lower_tr <- data[,1] - data[,2]
  Y_upper_tr <- data[,1] + data[,2]
  data_train <- list(X_lower_tr,X_upper_tr,Y_lower_tr,Y_upper_tr)
  
  data <- list(data_train,data_test)
  return(data)
}


## data Configuration II where both X^c and X^r are related to both Y^c and Y^r.
data_generate_II <- function(n,p,out.p,err,out.index){
  
  X_lower <- matrix(0, n, p)
  Y_lower <- NULL
  X_upper <- matrix(0, n, p)
  Y_upper <- NULL
  X_c <- matrix(0, n, p)
  X_r <- matrix(0, n, p)
  e_c <- NULL
  e_r <- NULL
  Y_c <- NULL
  Y_r <- NULL
  
  for (k in 1:n) {
    
    for (j in 1:p) {
      X_c[k, ] <- runif(p, 0, 80)
      X_r[k, ] <- runif(p, 5, 10)
    }
    
    ### parameters 
    beta_c <- runif((2*p+1), 2, 2.5)
    beta_r <- runif((2*p+1), 3, 3.5)
    
    ### errors
    e_c[k] <- runif(1,-err,err)
    e_r[k] <- runif(1,0, err)
    
    X_cr <- cbind(X_c,X_r)
    
    Y_c[k] <- c(1, X_cr[k, ]) %*% beta_c + e_c[k]
    Y_r[k] <- c(1, X_cr[k, ]) %*% beta_r + e_r[k]
  }
  
  X_lower <- X_c - X_r
  X_upper <- X_c + X_r
  Y_lower <- Y_c - Y_r
  Y_upper <- Y_c + Y_r
  
  train_ind <- sample(seq_len(n), size = n*2/3 )
  
  data_test <- list(X_lower[-train_ind,],X_upper[-train_ind,],Y_lower[-train_ind],Y_upper[-train_ind])
  
  S_Y <- sd(Y_c[train_ind])
  data <- cbind(Y_c,Y_r,X_c,X_r)[train_ind,]
  data <- data[order(data[,1],decreasing=F),]
  
  if (out.p>0) {
    n0 <- n *2/3 * out.p;
    for (i in  1:n0) {
      data[i,out.index] <- data[i,out.index] + 9*S_Y
    }
  }
  
  X_lower_tr <- data[,3:(3+p-1)] - data[,(3+p):(3+2*p-1)]
  X_upper_tr <- data[,3:(3+p-1)] + data[,(3+p):(3+2*p-1)]
  Y_lower_tr <- data[,1] - data[,2]
  Y_upper_tr <- data[,1] + data[,2]
  data_train <- list(X_lower_tr,X_upper_tr,Y_lower_tr,Y_upper_tr)
  
  data <- list(data_train,data_test)
  return(data)
}


#Configuration I1
#data_generate_I(n,p,out.p=0,err=20,out.index=1)

#Configuration I2
#data_generate_I(n,p,out.p=0,err=40,out.index=1)

#Configuration I3
#data_generate_I(n,p,out.p,err=20,out.index=1)

#Configuration I4
#data_generate_I(n,p,out.p,err=20,out.index=2)

#Configuration II1
#data_generate_II(n,p,out.p=0,err=20,out.index=1)

#Configuration II2
#data_generate_II(n,p,out.p=0,err=40,out.index=1)

#Configuration II3
#data_generate_II(n,p,out.p,err=20,out.index=1)

#Configuration II4
#data_generate_II(n,p,out.p,err=20,out.index=2)
