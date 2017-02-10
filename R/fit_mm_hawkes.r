#' Functions from the package HiddenMarkov
#'
#'




BW <- function (object, control = bwcontrol(), ...) {
    tau <- object$tau[-1] - object$tau[-length(object$tau)] ## times minus the inital zero an scaled
    Q <- object$Q ## transition matrix vals
    delta <- object$delta ## delta vals
    lambda <- object$lambda
    ## param.mu <- object$param.mu ##MY ADDITION initial mu param vals
    ## param.alpha <- object$param.alpha ##MY ADDITION initial alpha param vals
    ## param.beta <- object$param.beta ##MY ADDITION initial beta param values
    tol <- control$tol
    m <- nrow(Q) ## number states
    n <- length(tau) ## number of events
    oldLL <- -Inf ## setting initial log likelihood
    for (iter in 1:control$maxiter) { ## number of iterations to do
        cond <- Estep.mmpp(tau, Q, delta, lambda) ##retuns conditional expectations and log likelihood
        LL <- cond$LL ## extract log likelihood 
        diff <- LL - oldLL ## calculate difference
        if (control$prt) {
            cat("iter =", iter, "\n")
            cat("LL =", formatC(LL, digits = log10(1/tol) + 2, 
                                format = "f"), "\n")
            cat("diff =", diff, "\n\n")
        } ## if you want vals to be printed for each iteration
        if (diff < 0 & control$posdiff) 
            stop("Worse log-likelihood on last iteration") ## stop if log likelihood heading in wrong direction
        if (eval(control$converge)) 
            break
        Q <- Q * (diag(1/diag(cond$A)) %*% cond$A) ## Ck matrix art 320goes on about??
        diag(Q) <- 0
        diag(Q) <- -apply(Q, MARGIN = 1, FUN = sum) ## update Q?
        lambda <- cond$B/diag(cond$A) ## update baselines?
        if (object$nonstat) 
            delta <- exp(cond$logalpha[1, ] + cond$logbeta[1, ])
        else delta <- compdelta(solve(diag(lambda) - Q) %*% diag(lambda)) ## update Delta
        oldLL <- LL ## reset log likelihood
    }
    x <- list(delta = delta, Q = Q, lambda = lambda, LL = LL, 
              iter = iter, diff = diff) ## return results
    nms <- names(x)
    k <- length(nms)
    for (i in 1:k) object[[nms[i]]] <- x[[nms[i]]]
    return(object)
}
###
Q <- matrix(c(-2, 2,1, -1),byrow=TRUE, nrow=2)/10
library(HiddenMarkov)
x <- mmpp(NULL, Q, delta=c(0, 1), lambda=c(5, 1))
object <- simulate(x, nsim=50, seed=5)
object$param.mu <- object$lambda
object$param.beta <- c(0.8,0.8)
object$param.alpha <- c(0.6,0.6)


forwardback.mmpp <- function (tau, Q, delta, lambda) {
    m <- nrow(Q)
    n <- length(tau)
    Lambda <- diag(lambda)
    decomp <- eigen(Q - Lambda, symmetric = FALSE)
    if (any(duplicated(decomp$values))) 
        stop("repeated eigenvalues")
    S <- decomp$vectors
    Sinv <- solve(S)
    eigenval <- decomp$values
    phi <- as.double(delta)
    logalpha <- matrix(as.double(rep(0, m * (n + 1))), nrow = (n + 
        1))
    logalpha[1, ] <- log(phi)
    scalefac <- as.double(rep(0, n))
    post <- Sinv %*% Lambda
    psi <- array(as.double(0), dim = c(n, m, m))
 
    for (i in 1:n) {
        psi[i, , ] <- S %*% diag(exp(eigenval * tau[i])) %*% 
            post
        phi <- phi %*% psi[i, , ]
        sumphi <- sum(phi)
            scalefac[i] <- log(sumphi)
        phi <- phi/sumphi
        logalpha[i + 1, ] <- log(phi)
    }
    logbeta <- matrix(rep(as.double(0), m * (n + 1)), nrow = (n + 
        1))
    phi <- rep(as.double(1/m), m)
 
        lscale <- log(m)
        logck <- 0
        for (i in seq(n, 1, -1)) {
            phi <- psi[i, , ] %*% phi
            logck <- logck + scalefac[i]
            logbeta[i, ] <- log(phi) + lscale - logck
            sumphi <- sum(phi)
            phi <- phi/sumphi
            lscale <- lscale + log(sumphi)
        }
  
    return(list(logalpha = logalpha, logbeta = logbeta, eigenval = eigenval, 
        S = S, Sinv = Sinv, scalefac = scalefac, LL = sum(scalefac)))
}

Estep.mmpp <- function (tau, Q, delta, lambda, fortran = TRUE) 
{
    m <- ncol(Q)
    n <- length(tau)
    x <- forwardback.mmpp(tau, Q, delta, lambda)
    logbeta <- x$logbeta
    logalpha <- x$logalpha
    alpha <- exp(logalpha)
    beta <- exp(logbeta)
    d <- x$eigenval
    S <- x$S
    Sinv <- x$Sinv
    diff <- outer(d, d, FUN = "-") + diag(m)
    post0 <- Sinv %*% diag(lambda)
    A <- matrix(as.double(0), nrow = m, ncol = m)
    TT <- array(as.double(0), dim = c(n, m, m))
    
    for (i in 1:n) {
        expdtau <- exp(d * tau[i])
        difftau <- outer(expdtau, expdtau, FUN = "-")
            TT[i, , ] <- (difftau + diag(tau[i] * expdtau))/diff/exp(x$scalefac[i])
    }
    
    
    for (k in 1:m) {
        pre <- S %*% diag(Sinv[, k])
        for (j in 1:m) {
                post <- diag(S[j, ]) %*% post0
                for (i in 1:n) {
                    A[k, j] <- A[k, j] + alpha[i, ] %*% pre %*% 
                        TT[i, , ] %*% post %*% beta[(i + 1), ]
                }
        }
    }
   
    if (n == 1) 
        B <- exp(logalpha[-1, ] + logbeta[-1, ])
    else B <- apply(exp(logalpha[-1, ] + logbeta[-1, ]), MARGIN = 2, 
                    FUN = sum)
    return(list(A = A, B = B, logalpha = logalpha, logbeta = logbeta, 
                LL = x$LL))
}
