#' simulate two point patterns with a shared stochastic structure using a spde model for the laten field
#'
#' @return a list of length two matrices containing the point pattern locations
#'
#'
#'
#' @export
#'
sim.joint <- function(mesh = NULL, field.shared = list(mu = 1, kappa = 1, sigma2 = 1,seed = 1), beta = 1){
    mu1 <- field.shared[["mu"]]
    k1 <- field.shared[["kappa"]]
    s21 <- field.shared[["sigma2"]]
    seed1 <- field.shared[["seed"]]
    spde <- inla.spde2.matern(mesh = mesh, alpha=2)
    proj<-inla.mesh.projector(mesh = mesh)
    theta1 <- c(-0.5 * log(4 * pi * s21 * k1^2), log(k1))
    Q1 <- inla.spde2.precision(spde, theta1)
    ## the samples at the mesh nodes
    x1 <-  mu1 + inla.qsample(Q = Q1, seed = seed1, constr = spde$f$extraconstr)
    x2 <-  beta*x1 #+ x2
    x1 <- rpoispp(as.im(exp(matrix(inla.mesh.project(proj, x1),
                                   length(proj$x),length(proj$y)))/area(ripras(mesh$loc))))
    x2 <- rpoispp(as.im(exp(matrix(inla.mesh.project(proj, x2),
                                   length(proj$x),length(proj$y)))/area(ripras(mesh$loc))))
    return(list(pp1 = cbind(x1$x,x1$y), pp2 = cbind(x2$x,x2$y)))
}
