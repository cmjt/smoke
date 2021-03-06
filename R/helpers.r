## simple extraction of scottish county

county <- function(sp = NULL, name = NULL){
    county <- sp[sp$NAME == name, ]
    county
}
    
##extract values of random field at specified locations

extract <- function(x = NULL, locs = NULL,mesh = NULL, t = NULL, sd = FALSE){
    fields <- summary(x)$random.names[summary(x)$random.model == 
                                                "SPDE2 model" | summary(x)$random.model == "Copy"]
    idx <- which(summary(x)$random.model == "SPDE2 model" | summary(x)$random.model == 
                                                                      "Copy")
    n <- length(fields)
    A <- inla.mesh.project(mesh = mesh,loc = locs)$A
    if (!is.null(t)) {
        means <- list()
        for (i in idx[1]:idx[n]) {
            means[[i - idx[1] + 1]] <- lapply(1:t, function(j) {
                meshNodes <- x$summary.random[[i]]$mean[1:mesh$n + (i-1)*mesh$n]
                means[[i]] <- drop(A%*%meshNodes)
            })
        }
        sds <- list()
        for (i in idx[1]:idx[n]) {
            sds[[i - idx[1] + 1]] <- lapply(1:t, function(j) {
                meshNodes <- x$summary.random[[i]]$mean[1:mesh$n + (i-1)*mesh$n]
                sds[[i]] <- drop(A%*%meshNodes)
            })
        }
    } else {
        means <- list()
        for (i in idx[1]:idx[n]) {
            meshNodes <- x$summary.random[[i - idx[1] + 1]]$mean
            means[[i - idx[1] + 1]] <-  drop(A%*%meshNodes)
         
        }
        sds <- list()
        for (i in idx[1]:idx[n]) {
            meshNodes <- x$summary.random[[i - idx[1] + 1]]$sd
            sds[[i - idx[1] + 1]] <-  drop(A%*%meshNodes)
         
        }    
    }
    names(sds) <- names(means) <- fields
    ifelse(sd, return(sds), return(means))
}

   
## fit spatial only  marked model

fit.smoke.spatial.joint <- function(mesh = NULL, locs=NULL,  mark = NULL, verbose = FALSE,
                                    hyper = list(theta=list(prior='normal', param = c(0,10))),
                                    control.inla=list(strategy='gaussian',int.strategy = 'eb'), ...){
    n.marks <- length(table(mark))
    ## Error if there are more than two marks
    if (n.marks != 2){
        stop("Number of marks should be 2")
    }
    spde <- inla.spde2.matern(mesh, alpha=2)
    ##create sub point patterns from the mark 
    locs.s <- subset(locs,mark == 1)
    locs.c <-subset(locs,mark == 0)
   
    Ast.s <- inla.spde.make.A(mesh = mesh, loc = locs.s)
    Ast.c <- inla.spde.make.A(mesh = mesh, loc = locs.c)
    
    ##number of spatial mesh verticies
    m <- spde$n.spde                 
    n.s <- nrow(locs.s)                                  
    n.c <- nrow(locs.c)
    ## response for the LGCP model and ``effect'' size for poisson likelihood
                                       
    y1 <- rep(0:1, c( m, n.s))
    expected.s <-c(diag(spde$param.inla$M0), rep(0,n.s))
    y2 <- rep(0:1, c( m, n.c))
    expected.c <- c(diag(spde$param.inla$M0), rep(0,n.c))
                                       
    stk0 <- inla.stack(data=list(y = cbind(y1,NA),expect = expected.s),
                       A = list(rBind(Diagonal(n = m), Ast.s),1),tag="smoke",
                       effects=list(field1 = 1:m, a0 = rep(1,n.s + m)))

    stk1 <- inla.stack(data=list(y=cbind(NA,y2), expect=expected.c),
                       A=list(rBind(Diagonal(n = m), Ast.c),
                             rBind(Diagonal(n = m), Ast.c),1),tag = "crave",
                       effects=list(b_field1 = 1:m,field2 = 1:m, b0 = rep(1,n.c+m)))

    ##cojoin stacks
    stack = inla.stack(stk0,stk1)

    formula = y ~ -1  + a0 + b0 +
        f(field1, model = spde) +
        f(b_field1, copy = "field1", fixed = FALSE, hyper = hyper) 

    result = inla(formula, data=inla.stack.data(stack),
                  family=c("poisson","poisson"),only.hyperparam = FALSE,
                  control.predictor=list(A=inla.stack.A(stack),compute=TRUE),
                  control.compute = list(config = TRUE),
                  control.inla = control.inla,
                  verbose = verbose,
                  ...)
}


## estimated intensity for the simple univariate hawkes process
## note for intensity as suggested by http://radhakrishna.typepad.com/simulation-of-univariate-hawkes-via-thinning-1.pdf
## times <= p

## included hawkes process with general immigration

hawke.intensity <- function(mu, alpha, beta,times,p = NULL, n.s = NULL,gamma = NULL,states = NULL){
    if(is.null(p)) p <- times
    lam <- function(p){
        mu + alpha*sum(exp(-beta*(p-times))[times<p])
    }
    if(is.null(states)){
        lam.p <- rep(0,length(p))
        for(i in 1:length(p)){
            lam.p[i] <- lam(p[i])
        }
    }
    if(!is.null(states)){
        mu.t <- mu; alpha.t <- alpha; beta.t <- beta; times.t <- times ## temp storage
        lam.p <- rep(0,length(times))
        for(i in 1:length(times.t)){
            s <- states[i]
            mu <- mu.t[s]; alpha <- alpha.t[s]; beta <- beta.t[s]; times <- times.t[states==s]
            lam.p[i] <- lam(p = times.t[i])
        }
    }
    if(!is.null(n.s)){
        if(is.null(gamma))stop("intensity jump, gamma, due to external event must be supplied for a general immigrant process")
        if(n.s > length(p))stop("general immigration cannot exceed time span")
        lam.II <- function(s){
            gamma*sum(exp(-beta*(s-times))[times<s])
        }
        s <- sample(p,n.s)
        lam.s <- rep(0, length(s))
        for(i in 1:length(s)){
            lam.s[i] <- lam.II(s[i])
        }
        lam.p[which(p%in%s)]<- lam.p[which(p%in%s)] + lam.s
    }
    lam.p
}
    
## hawke intensity for marked process (same as ETAS model i.e., exponential mark)
hawke.mark.intensity <- function(mu, alpha, beta, nu, delta, times, marks, p = NULL, k = NULL){
    if(is.null(p)) p <- times
    if(is.null(k)) k <- marks
    if(length(p)!=length(k)) stop("k must be the same length as p")
    lam <- function(p,k){
        (mu + alpha*sum(exp(nu*marks)[times < p]*exp(-beta*(p-times))[times < p]))*delta*exp(-delta*k)
    }
    lam.p <- rep(0, length(p))
    for(i in 1:length(p)){
        lam.p[i] <- lam(p[i],k[i])
    }
    lam.p
}

## hawke intensity for marked process (same as ETAS model i.e., exponential mark)
hawke.population.intensity <- function(mu, alpha, beta, times, population, p = NULL, k = NULL){
    if(is.null(p)) p <- times
    if(is.null(k)) k <- population
    if(length(p)!=length(k)) stop("k must be the same length as p")
    lam <- function(p,k){
        (mu + alpha*sum(exp(-beta*(p-times))[times<p]))/k
    }
    lam.p <- rep(0, length(p))
    for(i in 1:length(p)){
        lam.p[i] <- lam(p[i],k[i])
    }
    lam.p
}


## hawke log likelihood
loglik.hawkes <- function(params, times){ 
    mu_i <- params[1]
    alpha_i <- params[2] 
    beta_i <- params[3]
    n <- length(times)
    term_1 <- -mu_i*times[n]
    term_2 <- sum(alpha_i/beta_i*(exp( -beta_i * (times[n] - times)) - 1))
    Ai <- c(0, sapply(2:n, function(z) {
        sum(exp( -beta_i * (times[z]- times[1:(z - 1)])))
    }))
    term_3 <- sum(log( mu_i + alpha_i * Ai)) 
    return(-term_1- term_2 -term_3)
}

## hawke log likelihood multiplied by population
loglik.hawkes.population <- function(params, times, population){ 
    mu_i <- params[1]
    alpha_i <- params[2] 
    beta_i <- params[3]
    n <- length(times)
    term_1 <- -mu_i*times[n]
    term_2 <- sum(alpha_i/beta_i*(exp( -beta_i * (times[n] - times)) - 1))
    Ai <- c(0, sapply(2:n, function(z) {
        sum(exp( -beta_i * (times[z]- times[1:(z - 1)])))
    }))
    term_3 <- sum(log((mu_i + alpha_i * Ai)/population)) 
    return(-term_1- term_2 -term_3)
}


## marked hawke log likelihood
loglik.marked.hawkes <- function(params, times,marks){ 
    mu_i <- params[1]
    alpha_i <- params[2] 
    beta_i <- params[3]
    nu_i <- params[4]
    delta_i <- params[5]
    n <- length(times)
    term_1 <- -mu_i*times[n]
    term_2 <- sum(alpha_i/beta_i*(exp(nu_i*marks)*(exp( -beta_i * (times[n] - times)) - 1)))
    Ai <- c(0, sapply(2:n, function(z) {
       sum(exp( nu_i* marks[z] -beta_i * (times[z]- times[1:(z - 1)])))
    }))
    marks.t <- c(0,marks)
    term_3 <- sum(log(mu_i + alpha_i*Ai*delta_i*exp(-delta_i*marks.t)))
    return(-term_1- term_2 -term_3)
}

## Haversine
hdist <- function (lon1, lat1, lon2, lat2){
    rad <- pi/180
    a1 <- lat1 * rad
    a2 <- lon1 * rad
    b1 <- lat2 * rad
    b2 <- lon2 * rad
    dlon <- b2 - a2
    dlat <- b1 - a1
    a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
    c <- 2 * atan2(sqrt(a), sqrt(1 - a))
    R <- 6378.145
    d <- R * c
    return(d)
}

int.hawkes <- function(mu,alpha,beta, times,t0, T){
    n <- length(times)
    result <- sapply(1:n, function(z){
        mu*times[z] + sum(alpha/beta*(1-exp(-beta*(times[z]-times[1:z]))))
    })
    lower <- which(times==t0)
    upper <- which(times==T)
    return(result[upper] - result[lower])
}
