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
                meshNodes <- x$summary.random$field.pp$mean[1:mesh$n + (i-1)*mesh$n]
                means[[i]] <- drop(A%*%meshNodes)
                return(r)
            })
        }
        sds <- list()
        for (i in idx[1]:idx[n]) {
            sds[[i - idx[1] + 1]] <- lapply(1:t, function(j) {
                meshNodes <- x$summary.random$field.pp$mean[1:mesh$n + (i-1)*mesh$n]
                sds[[i]] <- drop(A%*%meshNodes)
                return(r)
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
                                    control.inla=list(strategy='gaussian',int.strategy = 'eb')){
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
                       A = list(rBind(Diagonal(n = m), Ast.s)),tag="smoke",
                       effects=list(field1 = 1:m))

    stk1 <- inla.stack(data=list(y=cbind(NA,y2), expect=expected.c),
                       A=list(rBind(Diagonal(n = m), Ast.c),
                              rBind(Diagonal(n = m), Ast.c)),tag = "crave",
                       effects=list(b_field1 = 1:m,field2 = 1:m))

    ##cojoin stacks
    stack = inla.stack(stk0,stk1)

    formula = y ~ -1 +
        f(field1, model = spde) +
        f(field2, model = spde)+
        f(b_field1, copy = "field1", fixed = FALSE, hyper = hyper ) 

    result = inla(formula, data=inla.stack.data(stack),
                  family=c("poisson","poisson"),only.hyperparam = FALSE,
                  control.predictor=list(A=inla.stack.A(stack),compute=TRUE),
                  control.compute = list(config = TRUE),
                  control.inla = control.inla,
                  verbose = verbose)
}


## estimated intensity for the simple univariate hawkes process
## note for intensity as suggested by http://radhakrishna.typepad.com/simulation-of-univariate-hawkes-via-thinning-1.pdf
## times <= p


hawke.intensity <- function(mu, alpha, beta,times,p = NULL){
    if(is.null(p)) p <- times
    lam <- function(p){
        mu + alpha*sum(exp(-beta*(p-times))[times<p])
    }
    lam.p <- rep(0, length(times))
    for(i in 1:length(p)){
        lam.p[i] <- lam(p[i])
    }
    lam.p
}
