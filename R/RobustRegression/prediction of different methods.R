
######  input: data, a list of length 2
### data[[1]]: the training set, a list of length 4
### data[[2]]: the test set, a list of length 4
# data[[1]][[1]] means the matrix of lower bounds of independent variables of the training set
# data[[1]][[2]] means the matrix of upper bounds of independent variables of the training set
# data[[1]][[3]] means the vector of lower bounds of dependent variable of the training set
# data[[1]][[4]] means the vector of upper bounds of dependent variable of the training set

# data[[2]][[1]] means the matrix of lower bounds of independent variables of the test set
# data[[2]][[2]] means the matrix of upper bounds of independent variables of the test set
# data[[2]][[3]] means the vector of lower bounds of dependent variable of the test set
# data[[2]][[4]] means the vector of upper bounds of dependent variable of the test set

######  output: 1) parameters estimated obtained from the training set
#############   2) RMSE_L, RMSE_U, RMSE_H, and AR calculated from the test set


############## p=1 for CRM
CRM_pre_1 <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- CRM_1(traindata)
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  X_2r_test <- X_upper_test - X_lower_test
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  Y_2r_test <- Y_upper_test - Y_lower_test
  
  
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_r_fit <- as.vector((cbind(1,X_2r_test)%*%coefficients.R)/2)
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  return(result)
}
############## p=3 for CRM
CRM_pre_3 <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- CRM_3(traindata)
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  X_2r_test <- X_upper_test - X_lower_test
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  Y_2r_test <- Y_upper_test - Y_lower_test
  
  
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_r_fit <- as.vector((cbind(1,X_2r_test)%*%coefficients.R)/2)
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  
  return(result)
}
############## p=6 for CRM
CRM_pre_6 <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- CRM_6(traindata)
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  X_2r_test <- X_upper_test - X_lower_test
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  Y_2r_test <- Y_upper_test - Y_lower_test
  
  
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_r_fit <- as.vector((cbind(1,X_2r_test)%*%coefficients.R)/2)
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  
  return(result)
}


############## p=1 for CCRM
CCRM_pre_1 <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- CCRM_1(traindata)
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  X_2r_test <- X_upper_test - X_lower_test
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  Y_2r_test <- Y_upper_test - Y_lower_test
  
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_r_fit <- as.vector((cbind(1,X_2r_test)%*%coefficients.R)/2)
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  
  return(result)
}
############## p=3 for CCRM
CCRM_pre_3 <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- CCRM_3(traindata)
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  X_2r_test <- X_upper_test - X_lower_test
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  Y_2r_test <- Y_upper_test - Y_lower_test
  
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_r_fit <- as.vector((cbind(1,X_2r_test)%*%coefficients.R)/2)
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  
  return(result)
}
############## p=6 for CCRM
CCRM_pre_6 <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- CCRM_6(traindata)
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  X_2r_test <- X_upper_test - X_lower_test
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  Y_2r_test <- Y_upper_test - Y_lower_test
  
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_r_fit <- as.vector((cbind(1,X_2r_test)%*%coefficients.R)/2)
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  
  return(result)
}


############## CCRJM
CCRJM_pre <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- CCRJM(traindata)
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  X_test <-cbind(1,X_c_test,X_r_test)
  
  Y_c_fit <- as.vector(X_test%*%coefficients.C)
  Y_r_fit <- as.vector(X_test%*%coefficients.R)
  
  Y_l_fit <- (Y_c_fit - Y_r_fit)
  Y_u_fit <- (Y_c_fit + Y_r_fit)
  
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  
  
  return(result)
}


############## ID
ID_pre <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- ID(traindata)
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  
  Y_c_fit <- as.vector(X_c_test%*%coefficients.C[-1])+coefficients.C[1]
  Y_r_fit <- as.vector(X_r_test%*%coefficients.R)
  
  Y_l_fit <- (Y_c_fit - Y_r_fit)
  Y_u_fit <- (Y_c_fit + Y_r_fit)
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i]))^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  
  return(result)
}


############## SSLR
SSLR_pre <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <-SSLR(traindata)
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_2r_test <- X_upper_test - X_lower_test
  X_r_test <- X_2r_test/2
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_2r_test <- Y_upper_test - Y_lower_test
  Y_r_test <- Y_2r_test/2
  
  
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_2r_fit <- as.vector((cbind(1,X_2r_test)%*%coefficients.R))
  Y_r_fit <- Y_2r_fit/2
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  
  return(result)
}


############## IRR
IRR_pre <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- IRR(traindata)
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_2r_test <- X_upper_test - X_lower_test
  X_r_test <- X_2r_test/2
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_2r_test <- Y_upper_test - Y_lower_test
  Y_r_test <- Y_2r_test/2
  
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_2r_fit <- as.vector(cbind(1,X_2r_test)%*%coefficients.R)
  Y_r_fit <- Y_2r_fit/2
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  
  return(result)
}


############## IQR
IQR_pre <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <-IQR(traindata)
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_r_fit <- as.vector((cbind(1,X_r_test)%*%coefficients.R))
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  
  return(result)
  
}


################ different hyper-parameter estimator for iETKRR such as

## hyper-parameter estimators (S1,S1) for iETKRR method in midpoints and in ranges, respectively
iETKRR11_pre <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- iETKRR11(traindata)
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  
  
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_r_fit <- as.vector(cbind(1,X_r_test)%*%coefficients.R)
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  return(result)
}

## hyper-parameter estimators (S3,S1) for iETKRR method in midpoints and in ranges, respectively
iETKRR31_pre <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- iETKRR31(traindata)
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  
  
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_r_fit <- as.vector(cbind(1,X_r_test)%*%coefficients.R)
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  return(result)
}

## hyper-parameter estimators (S1,S3) for iETKRR method in midpoints and in ranges, respectively
iETKRR13_pre <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- iETKRR13(traindata)
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  
  
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_r_fit <- as.vector(cbind(1,X_r_test)%*%coefficients.R)
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  return(result)
}

## hyper-parameter estimators (S1,S4) for iETKRR method in midpoints and in ranges, respectively
iETKRR14_pre <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- iETKRR14(traindata)
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  
  
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_r_fit <- as.vector(cbind(1,X_r_test)%*%coefficients.R)
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  return(result)
}

## hyper-parameter estimators (S3,S4) for iETKRR method in midpoints and in ranges, respectively
iETKRR34_pre <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- iETKRR34(traindata)
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  
  
  Y_c_fit <- as.vector(cbind(1,X_c_test)%*%coefficients.C)
  Y_r_fit <- as.vector(cbind(1,X_r_test)%*%coefficients.R)
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  return(result)
}


############## Fuzzy-IR
Fuzzy_pre <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- Fuzzy(traindata)
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  X_test <-cbind(X_c_test,X_r_test,1)
  
  Y_c_fit <- as.vector(X_test%*%coefficients.C)
  Y_r_ln_fit <- as.vector(X_test%*%coefficients.R)
  Y_r_fit <- exp(Y_r_ln_fit)
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  return(result)
}


############## LN-IRR
LN_IRR_pre <- function(dataset){
  
  traindata <- dataset[[1]]
  testdata <- dataset[[2]]
  trainresult <- LN_IRR(traindata)
  coefficients.C <- trainresult$coefficients.C
  coefficients.R <- trainresult$coefficients.R
  
  X_lower_test <- as.matrix(testdata[[1]])
  X_upper_test <- as.matrix(testdata[[2]])
  Y_lower_test <- testdata[[3]]
  Y_upper_test <- testdata[[4]]
  n_test <- nrow(X_lower_test)
  p <- ncol(X_lower_test)
  X_c_test <- (X_lower_test + X_upper_test)/2
  X_r_test <- (X_upper_test - X_lower_test)/2
  Y_c_test <- (Y_lower_test + Y_upper_test)/2
  Y_r_test <- (Y_upper_test - Y_lower_test)/2
  X_test <-cbind(1,X_c_test,X_r_test)
  
  Y_c_fit <- as.vector(X_test%*%coefficients.C)
  Y_r_ln_fit <- as.vector(X_test%*%coefficients.R)
  Y_r_fit <- exp(Y_r_ln_fit)
  
  Y_l_fit <- Y_c_fit - Y_r_fit
  Y_u_fit <- Y_c_fit + Y_r_fit
  
  
  RMSE_l <- dist(rbind(Y_l_fit, Y_lower_test)) / sqrt(n_test)
  RMSE_u <- dist(rbind(Y_u_fit, Y_upper_test)) / sqrt(n_test)
  
  #### RMSE_h
  sum = 0
  for (i in 1:n_test) {
    sum = sum + (abs(Y_c_fit[i] - Y_c_test[i]) + abs(Y_r_fit[i] - Y_r_test[i])) ^2
  }
  RMSE_h <- sqrt(sum) / sqrt(n_test)
  
  #### AR
  sum1 = 0
  s <- NULL
  b <- NULL
  for (i in 1:n_test) {
    s[i] = max(min(Y_u_fit[i], Y_upper_test[i]) - max(Y_l_fit[i], Y_lower_test[i]),0)
    b[i] = max(Y_u_fit[i], Y_upper_test[i]) - min(Y_l_fit[i], Y_lower_test[i]) - max(max(Y_l_fit[i], Y_lower_test[i]) - min(Y_u_fit[i], Y_upper_test[i]),0)
    sum1 = sum1 + s[i] / b[i]
  }
  AR <- sum1 / n_test
  
  #### MAE
  MAE_L <- median(abs(Y_l_fit - Y_lower_test))
  MAE_U <- median(abs(Y_u_fit - Y_upper_test))
  MAE_C <- median(abs(Y_c_fit - Y_c_test))
  MAE_R <- median(abs(Y_r_fit - Y_r_test))
  
  
  result <- list("RMSE_l"=RMSE_l,"RMSE_u"=RMSE_u,"RMSE_h"=RMSE_h,"AR"=AR,
                 "MAE_L"=MAE_L,"MAE_U"=MAE_U,"MAE_C"=MAE_C,"MAE_R"=MAE_R,
                 "coe.C"=coefficients.C,"coe.R"=coefficients.R)
  return(result)
}
