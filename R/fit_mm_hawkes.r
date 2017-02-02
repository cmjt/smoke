#' Functions from the package HiddenMarkov
#'
#'

FB.mine <- function(times, mu,  alpha, beta, Q, pi){
    r <- nrow(Q)
    n <- length(times)
    Lambda <- diag(mu)

}

LL1 <- function(pi, D, Q, w){
    sum(log(pi)) - sum(D*q) + sum(w*log(q))
}
    


FB <- function (tau, Q, delta, lambda) {
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
    logalpha <- matrix(as.double(rep(0, m * (n + 1))), nrow = (n + 1))
    logalpha[1, ] <- log(phi)
    scalefac <- as.double(rep(0, n))
    post <- Sinv %*% Lambda
    psi <- array(as.double(0), dim = c(n, m, m))
    for (i in 1:n) {
        psi[i, , ] <- S %*% diag(exp(eigenval * tau[i])) %*% post
        phi <- phi %*% psi[i, , ]
        sumphi <- sum(phi)
        scalefac[i] <- log(sumphi)
        phi <- phi/sumphi
        logalpha[i + 1, ] <- log(phi)
    }
    logbeta <- matrix(rep(as.double(0), m * (n + 1)), nrow = (n +  1))
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

ES <- function (tau, Q, delta, lambda) {
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




BW <- function (object, control = bwcontrol(), ...) {
    tau <- object$tau[-1] - object$tau[-length(object$tau)]
    Q <- object$Q
    delta <- object$delta
    lambda <- object$lambda
    tol <- control$tol
    m <- nrow(Q)
    n <- length(tau)
    oldLL <- -Inf
    for (iter in 1:control$maxiter) {
        cond <- Estep.mmpp(tau, Q, delta, lambda)
        LL <- cond$LL
        diff <- LL - oldLL
        if (control$prt) {
            cat("iter =", iter, "\n")
            cat("LL =", formatC(LL, digits = log10(1/tol) + 2, 
                                format = "f"), "\n")
            cat("diff =", diff, "\n\n")
        }
        if (diff < 0 & control$posdiff) 
            stop("Worse log-likelihood on last iteration")
        if (eval(control$converge)) 
            break
        Q <- Q * (diag(1/diag(cond$A)) %*% cond$A)
        diag(Q) <- 0
        diag(Q) <- -apply(Q, MARGIN = 1, FUN = sum)
        lambda <- cond$B/diag(cond$A)
        if (object$nonstat) 
            delta <- exp(cond$logalpha[1, ] + cond$logbeta[1, 
                                                           ])
        else delta <- compdelta(solve(diag(lambda) - Q) %*% diag(lambda))
        oldLL <- LL
    }
    x <- list(delta = delta, Q = Q, lambda = lambda, LL = LL, 
              iter = iter, diff = diff)
    nms <- names(x)
    k <- length(nms)
    for (i in 1:k) object[[nms[i]]] <- x[[nms[i]]]
    return(object)
}
