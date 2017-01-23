#' Fitting a hawkes (temporal) point process model
#'
#' Fits a simple univariate point process model
#'
#' @param params a numeric vector of length three of starting values
#' for the paramerers mu, alpha, beta respectively.
#' @param times a numeric vector of times  
#' @return estimates parameter values
#' 
#' @export

fit.hawkes <- function(params, times){
    fit<- optim(params,loglik.hawkes, times = times)
    pars <- fit$par
    out <- matrix(NA,ncol = 2, nrow = 3)
    out[,1] <- pars
    rownames(out) <- c("mu","alpha","beta")
    colnames(out) <- c("param","se")
    out
}
