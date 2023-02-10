library(quantreg)

HF1_sub <- function(x,y)
{
    model = rq(y~x, method = "pfn")
    return(as.numeric(model$coefficients))
}
