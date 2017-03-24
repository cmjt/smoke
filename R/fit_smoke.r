#' Fitting a point process model
#'
#' Fits a spatio-temporal marked point process model
#' using INLA coupled with the SPDE approch 
#'
#' @return An inla model fit object, 
#' @param mesh delauney triangulation of area, an object returned by \link{make.mesh} is suitable.
#' @param locs a matrix of \code{nrow} locations in \code{ncol} dimesions.
#' @param mark a vector of length \code{nrow} of marks refering to each point location must be zeros and ones
#' @param t a numeric vector specifying a temporal index for each observation (starting at 1.....T).
#' @param verbose Logical, if \code{TRUE}, model fitting is output
#' the console.
#' @param hyper prior given to the joint interaction parameter, for a marked model only. By default
#' this is Normal(0,10)
#' @param prior.rho prior for the temporal correlation coefficient, by default a \code{pcprior} is used
#' with \code{param=c(0-0.9)}.
#' @param control.inla a list as supplied by call to \code{inla} see http://www.r-inla.org/
#' @importMethodsFrom Matrix diag
#' 
#' @export
fit.smoke <- function(mesh = NULL, locs=NULL,  mark = NULL, t = NULL, verbose = FALSE,
                      hyper = list(theta=list(prior='normal', param = c(0,10))),
                      prior.rho = list(theta = list(prior='pccor1', param = c(0, 0.9))),
                      control.inla=list(strategy='gaussian',int.strategy = 'eb'),
                      ...){
    
   if(!is.null(mark)){
       result <- fit.smoke.spatial.joint(mesh = mesh, locs = locs,  mark = mark, verbose = verbose,
                                    hyper = hyper, control.inla = control.inla, ...)
           }
    if(!is.null(t)){
        result <- lgcpSPDE:::fit.lgcp(mesh = mesh, locs = locs, temp = t, prior.rho = prior.rho,
                                      verbose = verbose, control.inla = control.inla)
        }
   result
}


#' fitting a multiple marked model, specific to Duke data
#' 
#' export
#'
fit.smoke.multi <- function(mesh = NULL, locs = NULL, verbose = FALSE,
                            h.b1 = list(theta=list(prior='normal', param = c(0,10))),
                            h.b2 = list(theta=list(prior='normal', param = c(0,10))),
                            h.g1 = list(theta=list(prior='normal', param = c(0,10))),
                            h.g2 = list(theta=list(prior='normal', param = c(0,10))),
                            prior.rho = list(theta = list(prior='pccor1', param = c(0, 0.9))),
                            control.inla=list(strategy='gaussian',int.strategy = 'eb'),
                            ...){
    spde <- inla.spde2.matern(mesh, alpha=2) ## maybe update to include priors on spatial field later on
    m <- spde$n.spde
    locs.s <- locs[["smoker"]]
    n.s <- nrow(locs.s)
    locs.ns <- locs[["non.smoker"]]
    n.ns <- nrow(locs.ns)
    locs.t <- locs[["tobacco"]]
    n.t <- nrow(locs.t)
    locs.nt <- locs[["non.tobacco"]]
    n.nt <- nrow(locs.nt)
    A.s <- inla.spde.make.A(mesh = mesh, loc = locs.s)
    A.ns <- inla.spde.make.A(mesh = mesh, loc = locs.ns)
    A.t <- inla.spde.make.A(mesh = mesh, loc = locs.t)
    A.nt <- inla.spde.make.A(mesh = mesh, loc = locs.nt)
    y.s <- rep(0:1, c( m, n.s))
    expected.s <-c(diag(spde$param.inla$M0), rep(0,n.s))
    y.ns <- rep(0:1, c( m, n.ns))
    expected.ns <- c(diag(spde$param.inla$M0), rep(0,n.ns))
    y.t <- rep(0:1, c( m, n.t))
    expected.t <-c(diag(spde$param.inla$M0), rep(0,n.t))
    y.nt <- rep(0:1, c( m, n.nt))
    expected.nt <-c(diag(spde$param.inla$M0), rep(0,n.nt))
    stk.s <- inla.stack(data=list(y = cbind(y.s,NA,NA,NA),expect = expected.s),
                       A = list(rBind(Diagonal(n = m), A.s),rBind(Diagonal(n = m), A.s)),tag="smoke",
                       effects=list(field1 = 1:m,s.copy.field2 = 1:m))
    stk.ns <- inla.stack(data=list(y = cbind(NA,y.ns,NA,NA),expect = expected.ns),
                       A = list(rBind(Diagonal(n = m), A.ns),rBind(Diagonal(n = m), A.ns)),tag="non-smoke",
                       effects=list(copy.field1 = 1:m,ns.copy.field2 = 1:m))
    stk.t <- inla.stack(data=list(y = cbind(NA,NA,y.t,NA),expect = expected.t),
                       A = list(rBind(Diagonal(n = m), A.t)),tag="tobacco",
                       effects=list(field2 = 1:m))
    stk.nt <- inla.stack(data=list(y = cbind(NA,NA,NA,y.nt),expect = expected.nt),
                       A = list(rBind(Diagonal(n = m), A.nt)),tag="non-tobacco",
                       effects=list(copy.field2 = 1:m))
    stack <- inla.stack(stk.s, stk.ns, stk.t, stk.nt)
    formula <- y ~ 0 + f(field1, model = spde) +
        f(copy.field1,copy = "field1",fixed = FALSE,hyper = h.b1) + 
        f(field2, model = spde) + f(copy.field2, copy = "field2",fixed = FALSE, hyper = h.b2) +
        f(s.copy.field2, copy = "field2",fixed = FALSE, hyper = h.g1) +
        f(ns.copy.field2, copy = "field2",fixed = FALSE, hyper = h.g2)
    result = inla(formula, data=inla.stack.data(stack),
                  family=rep("poisson",4),only.hyperparam = FALSE,
                  control.predictor=list(A=inla.stack.A(stack),compute=TRUE),
                  control.compute = list(config = TRUE),
                  control.inla = control.inla,
                  verbose = verbose,
                  ...)
}
        
        
