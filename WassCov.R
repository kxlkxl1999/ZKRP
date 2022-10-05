## This script takes in a list of distributional observations for a number
## of subjects and produces an overall Wasserstein covariance matrix, representing
## concordance across different densities for the entire cohort, as well as an
## individual Wasserstein covariance matrix for each subject.  These objects are
## those defined in Petersen and Mueller (2019), Wasserstein Covariance for Multiple 
## Random Densities, Biometrika

## As an alternative option, one can create only individual Wasserstein covariance
## matrices that replace the cohort-based Wasserstein means (one per density, computed
## across subjects) used in the paper with a single individual-based Wasserstein
## mean computed across densities.  This is an experimental quantity, but will still
## yield a valid covariance matrix.

## Required packages: frechet, fdadensity

## Inputs:
##
## y      - list of length n (number of subjects), y[[i]] is a list of length
##          p, with y[[i]][[j]] containing the observed sample for density j
##          of subject i.
## f      - list of length n (number of subjects), f[[i]] is a list of length 
##          p, with f[[i]][[j]] a vector of probability density function values 
##          for density j of subject i.
## fin    - list of length n (number of subjects), fin[[i]] is a list of length p,
##          with fin[[i]][[j]] the support grid for the density values in f[[i]][[j]].
## Q      - list of length n (number of subjects), Q[[i]] is a matrix of dimensions
##          p-by-M, with Q[[i]][j, ] containing the quantile function values corresponding
##          to density j of subject i.  All quantile functions for all densities and 
##          all subjects are assumed to be measured on a common grid Qin.
## Qin    - vector of length M containing ordered values spanning the interval [0, 1]
##          representing the support grid of all quantile functions.  If y, f, and Qin
##          are null, but Q is provided, then Qin defaults to an equally spaced grid
##          of length M, where M is inferred from the dimensions of Q.
## UseExp - T/F: Should the experimental individual Wasserstein covariance quantity be
##          used as opposed to the usual one?  Default is FALSE.

## Outputs:
##
## WC_ind - list of length n, giving the individual p-by-p Wasserstein covariance  
##          matrices for all subjects
## WC     - p-by-p covariance matrix for entire cohort, will be NULL if UseExp = TRUE

## Details
##
## At least one of y, f, or Q must be provided; if more than one is provided, only
## one will be used, with y taking precedence over f, and f over Q.
##
## If y is used, then a preliminary density estimation step will be performed on
## the samples using the CreateDensity function from the frechet package, which 
## executes histogram smoothing with a data-driven bandwidth choice.  This can 
## be time-consuming if there are many density samples. Once densities are available, 
## these are converted to quantile functions using the dens2quantile function in the
## fdadensity package. 


WassCov <- function(y = NULL, f = NULL, fin = NULL, Q = NULL, Qin = NULL, UseExp = FALSE){
  
  # perform input checks  
  
  if(is.null(y) && is.null(f) && is.null(Q)) {
    stop("At least one of y, f, or Q must be provided")
  }
  
  if(!is.null(y)) {
    
    inputType <- "samples"
    
    if(!is.null(f) || !is.null(Q)) {
      warning("As y is provided, other inputs f and/or Q will be ignored")
    }
    
    if(!is.list(y)) {
      stop("Input y must be a list, one element per subject")
    }
    
    if(!all(sapply(y, class) == "list")) {
      stop("All elements of y must be lists, one element/sample per density")
    }
    
    if(length(unique(lapply(y, length))) != 1) {
      stop("Must provide the same number of samples for each subject")
    }
    
    if(any(sapply(y, function(yy) sapply(yy, length)) < 10)) {
      warning("Some density samples have very few (< 10) observations")
    }
    
  } else if(!is.null(f)) {
    
    inputType <- "densities"
    
    if(!is.null(Q)) {
      warning("As f is provided, input Q will be ignored")
    }
    
    if(is.null(fin) || length(fin) != length(f)) {
      stop("fin is either not provided or has an incorrect format")
    }
    
    if(!is.list(f)) {
      stop("Input f must be a list, one element per subject")
    }
    
    if(!all(sapply(f, class) == "list")) {
      stop("All elements of f must be lists, one element/sample per density")
    }
    
    if(length(unique(lapply(f, length))) != 1) {
      stop("Must provide the same number of samples for each subject")
    }
    
  } else {
    
    inputType <- "quantiles"
    
    if(!is.list(Q)) {
      stop("Input Q must be a list")
    }
    
    if(any(lapply(Q, function(QQ) class(QQ)[1]) != "matrix")) {
      stop("Each element of Q must be a matrix")
    } else if(length(unique(sapply(Q, nrow))) != 1 || length(unique(sapply(Q, ncol))) != 1) {
      stop("Each element of Q must be a matrix of the same dimensions")
    }
    
    if(is.null(Qin)) {
      Qin <- seq(0, 1, length.out = ncol(Q[[1]]))
    } else if(ncol(Q[[1]]) != length(Qin)) {
      stop("Length of Qin must match number of columns in the matrices Q[[i]]")
    }
    
  }
  
  # Obtain Quantile Functions, if necessary
  
  if(inputType == "samples") {
    
    n <- length(y); p <- length(y[[1]])
    Qin <- seq(0, 1, length.out = 101)
    Q <- lapply(1:n, function(i){
      
      t(sapply(y[[i]], function(yy){
        
        tmp <- CreateDensity(y = yy)
        val <- dens2quantile(dens = tmp$y, dSup = tmp$x, qSup = Qin)
        val
        
      }))
      
    })
    
  } else if(inputType == "densities") {
    
    n <- length(f); p <- length(f[[1]]) 
    Qin <- seq(0, 1, length.out = 101)
    Q <- lapply(1:n, function(i){
      
      t(sapply(1:p, function(j){
        yy <- f[[i]][[j]]; xx <- fin[[i]][[j]]
        val <- dens2quantile(dens = yy, dSup = xx, qSup = Qin)
        val
      }))
      
    })
    
  }
  
# Center quantile functions
  
  M <- length(Qin); n <- length(Q); p <- nrow(Q[[1]])
  
  if(UseExp == FALSE) { # get the usual Wasserstein Covariance
    
    QMn <- sapply(1:p, function(j){
      Qj <- sapply(Q, function(q) q[j, ])
      rowMeans(Qj)
    })
    QCen <- lapply(Q, function(q) q - t(QMn))
    
  } else {
    
    QCen <- lapply(Q, function(q) scale(q, scale = FALSE))
    
  }

# Construct Individual Wasserstein Covariance Matrices
  
  WC_ind <- vector(mode = "list", length = n)
  for(i in 1:n){
    WC_ind[[i]] <- matrix(NA, nrow = p, ncol = p)
    for(j in 1:p) {
      for(k in j:p) {
        WC_ind[[i]][j, k] <- WC_ind[[i]][k, j] <- fdadensity:::trapzRcpp(X = Qin, Y = QCen[[i]][j, ]*QCen[[i]][k, ])
      }
    }
  }
  
# Construct overall WC Matrix, unless UseExp = TRUE
  
  WC <- NULL
  if(UseExp == FALSE) {
    WC <- matrix(NA, nrow = p, ncol = p)
    for(j in 1:p) {
      for(k in j:p) {
        tmp <- sapply(1:n, function(i) WC_ind[[i]][j, k])
        WC[j, k] <- WC[k, j] <- mean(tmp)
      }
    }
  }
  
  return(list(WC_ind = WC_ind, WC = WC))
}


