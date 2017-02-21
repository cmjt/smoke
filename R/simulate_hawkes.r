#' Function to simulate a simple univariate Hawkes process
#'
#' basic code based on simulation from http://radhakrishna.typepad.com/simulation-of-univariate-hawkes-via-thinning-1.pdf
#'
#' @param mu numeric integer specifying base rate of the hawkes process
#' @param alpha numeric integer specifying intensity jump after an event occurence
#' @param beta numeric integer specifying exponential intensity decay
#' @param maxT numeric specifying end of the time line within which to simulat the proces, by default = 10
#' @param plot.proc logical if TRUE the process is plotted along with the intensity
#' @param seed by default = 1
#'
#' @export
 
simulate.hawkes <- function(mu = NULL,alpha = NULL,beta = NULL,maxT = 10, plot.proc = FALSE,seed = 1){
    if(alpha >= beta){stop("alpha must be less then beta to ensure intensity decreases quicker than new events increase it")}
    set.seed(seed)
    times <- numeric()
    s <- 0
    t <- 0
    lambda_star <- mu
    s <- s -log(runif(1))/lambda_star
    t <- s
    dlambda <- alpha
    times <- c(times, t)
    while(s < maxT) {
        U <- runif(1)
        s <- s -log(U)/lambda_star
        u <- runif(1)
        if(u <= (mu + dlambda*exp(-beta*(s-t)))/lambda_star){
            dlambda <- alpha + dlambda*exp(-beta*(s-t))
            lambda_star <- lambda_star + alpha
            t <- s
            times <- c(times, t)
        }
    }
    if(plot.proc){
        plot.hawkes(times, mu, alpha, beta)
        return{times}
    }else{return(times)}
}
