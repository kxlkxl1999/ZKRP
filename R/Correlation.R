symbolic.r <- function(model) {
    x = model$predictor
    y = model$response
    haty = model$fitted.values
    
    # sample deviation
    sample.dev <- function(lower, upper) {
        n = length(lower)
        s.cov = (sum(lower^2 + lower*upper + upper^2) / (3*n)) - (sum(lower + upper)^{2} / (4*n^{2}))
        
        return(s.cov)
    }
    
    n = nrow(x)
    meany = sum(y[, 1:2]) / (2*n)
    meanhaty = sum(haty[, 1:2]) / (2*n)
    
    # symbolic covariance
    sym.c = ( sum(2*(y[, 1] - meany)*(haty[, 1] - meanhaty) +
                      (y[, 1] - meany)*(haty[, 2] - meanhaty) +
                      (y[, 2] - meany)*(haty[, 1] - meanhaty) +
                      2*(y[, 2] - meany)*(haty[, 2] - meanhaty)) ) / ( 6*n )
    
    sym.r = sym.c / (sqrt(sample.dev(y[, 1], y[, 2])) * sqrt(sample.dev(haty[, 1], haty[, 2])))
    
    return(sym.r)
}