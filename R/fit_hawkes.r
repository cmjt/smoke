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

fit.hawkes <- function(params, times, marks = NULL){
    if(!is.null(marks)){
        fit <- optim(params,loglik.marked.hawkes, times = times, marks = marks)
        rownames <- c("mu","alpha","beta","nu","delta")
    }else{
        fit<- optim(params,loglik.hawkes, times = times)
        rownames <- c("mu","alpha","beta")
    }
    pars <- fit$par
    nrow <- ifelse(!is.null(marks),5,3)
    out <- matrix(NA,ncol = 2, nrow = nrow)
    out[,1] <- pars
    rownames(out) <- rownames
    colnames(out) <- c("param","se")
    out
}
