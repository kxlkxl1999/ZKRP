
method_compare <- function(dataset){
  
  result <- list()
  result[[1]] <- CRM_pre(dataset)
  result[[2]] <- CCRM_pre(dataset)
  result[[3]] <- CCRJM_pre(dataset)
  result[[4]] <- ID_pre(dataset)
  result[[5]] <- SSLR_pre(dataset)
  result[[6]] <- IRR_pre(dataset)
  result[[7]] <- IQR_pre(dataset)
  result[[8]] <- iETKRR_pre(dataset)
  result[[9]] <- Fuzzy_pre(dataset)
  result[[10]] <- LN_IRR_pre(dataset)
  
  result_RMSE_l <- c(result[[1]]$RMSE_l,result[[2]]$RMSE_l,result[[3]]$RMSE_l,result[[4]]$RMSE_l,
                     result[[5]]$RMSE_l,result[[6]]$RMSE_l,result[[7]]$RMSE_l,result[[8]]$RMSE_l,
                     result[[9]]$RMSE_l,result[[10]]$RMSE_l)
  
  result_RMSE_u <- c(result[[1]]$RMSE_u,result[[2]]$RMSE_u,result[[3]]$RMSE_u,result[[4]]$RMSE_u,
                     result[[5]]$RMSE_u,result[[6]]$RMSE_u,result[[7]]$RMSE_u,result[[8]]$RMSE_u,
                     result[[9]]$RMSE_u,result[[10]]$RMSE_u)
  
  result_RMSE_h <- c(result[[1]]$RMSE_h,result[[2]]$RMSE_h,result[[3]]$RMSE_h,result[[4]]$RMSE_h,
                     result[[5]]$RMSE_h,result[[6]]$RMSE_h,result[[7]]$RMSE_h,result[[8]]$RMSE_h,
                     result[[9]]$RMSE_h,result[[10]]$RMSE_h)
  
  result_AR <- c(result[[1]]$AR,result[[2]]$AR,result[[3]]$AR,result[[4]]$AR,result[[5]]$AR,
                 result[[6]]$AR,result[[7]]$AR,result[[8]]$AR,result[[9]]$AR,result[[10]]$AR)
  
  
  result_MAE_l <- c(result[[1]]$MAE_L,result[[2]]$MAE_L,result[[3]]$MAE_L,result[[4]]$MAE_L,
                    result[[5]]$MAE_L,result[[6]]$MAE_L,result[[7]]$MAE_L,result[[8]]$MAE_L,
                    result[[9]]$MAE_L,result[[10]]$MAE_L)
  
  result_MAE_u <- c(result[[1]]$MAE_U,result[[2]]$MAE_U,result[[3]]$MAE_U,result[[4]]$MAE_U,
                    result[[5]]$MAE_U,result[[6]]$MAE_U,result[[7]]$MAE_U,result[[8]]$MAE_U,
                    result[[9]]$MAE_U,result[[10]]$MAE_U)
  
  result_MAE_c <- c(result[[1]]$MAE_C,result[[2]]$MAE_C,result[[3]]$MAE_C,result[[4]]$MAE_C,
                    result[[5]]$MAE_C,result[[6]]$MAE_C,result[[7]]$MAE_C,result[[8]]$MAE_C,
                    result[[9]]$MAE_C,result[[10]]$MAE_C)
  
  result_MAE_r <- c(result[[1]]$MAE_R,result[[2]]$MAE_R,result[[3]]$MAE_R,result[[4]]$MAE_R,
                    result[[5]]$MAE_R,result[[6]]$MAE_R,result[[7]]$MAE_R,result[[8]]$MAE_R,
                    result[[9]]$MAE_R,result[[10]]$MAE_R)
  
  
  
  accuracy <- cbind(result_RMSE_l,result_RMSE_u,result_RMSE_h,result_AR,
                    result_MAE_l,result_MAE_u,result_MAE_c,result_MAE_r)
  rownames(accuracy) <- c("CRM","CCRM","CCRJM","ID","SSLR","IRR","IQR","iETKRR","Fuzzy","LN-IRR")
  colnames(accuracy) <- c("RMSE_l","RMSE_u","RMSE_h","AR",
                          "MAE_L","MAE_U","MAE_C","MAE_R")

  return(accuracy)
}


