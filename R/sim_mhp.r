#' Function to simulate a Markov-modulated Hawkes process with stepwise decay
#'
#'
#' @export
sim.mhp <- function(n = NULL,initial.state = NULL, Q = matrix(c(-0.5,0.2,0.5,-0.5),nrow = 2,byrow=TRUE),
                    mu = c(1,2), alpha = c(1,0.5), beta = c(2,0.8),seed = 1234,plot.proc = FALSE){
    set.seed(seed)
    n.points <- n
    n.states <- nrow(Q)
    #initialise in first state 
    j <- initial.state
    times <- numeric()
    states <- numeric()
    states <- c(states,j)
    t <- s <- 0
    times <- c(times,t)
    while(length(times) < n.points){
        q_j <- -Q[j,j]
        r_i <- q_j + mu[j] + alpha[j] * sum(exp(-beta[j]*(t -times))[times<t])
        tau_i <- rexp(1,r_i)
        s <- s + tau_i
        U <- runif(1)
        if(U > q_j/r_i) {t <- s; times <- c(times,t); states <- c(states,j)}
        if(U <= q_j/r_i) { k <- sample(1:n.states,1); y2 <- k ; j <- y2; s <- s+1; times <- times}
        
    }
    if(plot.proc){plot.hawkes(times, mu, alpha, beta, states = states)}else{return(cbind(times,states))}
}
