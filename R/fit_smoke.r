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
                      control.inla=list(strategy='gaussian',int.strategy = 'eb')){
    
   if(!is.null(mark)){
       result <- fit.smoke.spatial.joint(mesh = mesh, locs = locs,  mark = mark, verbose = verbose,
                                    hyper = hyper, control.inla = control.inla)
           }
    if(!is.null(time)){
        result <- lgcpSPDE:::fit.lgcp(mesh = mesh, locs = locs, temp = t, prior.rho = prior.rho,
                                      verbose = verbose, control.inla = control.inla)
        }
   result
}
