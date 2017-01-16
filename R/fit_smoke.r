#' Fitting a marked point process model
#'
#' Fits a spatio-temporal marked point process model
#' using INLA coupled with the SPDE approch 
#'
#' @return An inla model fit object, 
#' @param mesh delauney triangulation of area, an object returned by \link{make.mesh} is suitable.
#' @param locs a matrix of \code{nrow} locations in \code{ncol} dimesions.
#' @param mark a vector of length \code{nrow} of marks refering to each point location must be zeros and ones
#' @param verbose Logical, if \code{TRUE}, model fitting is output
#' the console.
#'
#' @importMethodsFrom Matrix diag
#' 
#' @export
fit.smoke <- function(mesh = NULL, locs=NULL,  mark = NULL, verbose = FALSE,
                      hyper = list(theta=list(prior='normal', param = c(0,10)))){
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
                       A = list(rBind(Diagonal(n = m), Ast.s), 1),tag="smoke",
                       effects=list(field1 = 1:m,
                                    list(b0=rep(1,m + n.s))))

    stk1 <- inla.stack(data=list(y=cbind(NA,y2), expect=expected.c),
                       A=list(rBind(Diagonal(n = m), Ast.c),
                              rBind(Diagonal(n = m), Ast.c),1),tag = "crave",
                       effects=list(b_field1 = 1:m,field2 = 1:m, list(a0 = rep(1, m + n.c))))

    ##cojoin stacks
    stack = inla.stack(stk0,stk1)

    formula = y ~ -1 + b0 +  a0  +
        f(field1, model = spde) +
        f(field2, model = spde)+
        f(b_field1, copy = "field1", fixed=FALSE, hyper = hyper ) 

    result = inla(formula, data=inla.stack.data(stack),
                  family=c("poisson","poisson"),only.hyperparam = FALSE,
                  control.predictor=list(A=inla.stack.A(stack),compute=TRUE),
                  control.compute = list(config = TRUE),
                  verbose = verbose)
}
