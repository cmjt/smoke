#' simulate smoking logging locations, times and satisfaction
#' *************time and space are independent ATM********************
#' @return a 2xn matrix of values with logging locations, hour of day (1:23), day of week (1:7), satisfaction scores (1:20)
#' and craving scores (1:20)
#' *****************need to add in parameters for tod, wday, etc. trends********************
#' @param x  either 1) a matrix of locations, 2) a spatial density object as \code{class} \code{im},
#' or 3) a model fit object as returned by \code{fit_smoke} (yet to add in).
#' 
#'
#' @export

sim.smoke <- function(x = NULL, sp = NULL, tod.prob = rep(1/23, times = 23), dow.prob = rep(1/7,times = 7)){
    if(class(x)=="im"){
        im <- x
    }else{
        if(is.matrix(x)){
            if(!is.null(sp)){w <- as.owin(sp)}else{w <- ripras(x)}
            x <- suppressWarnings(ppp(x[,1],x[,2],window = w))
            im <- density(x)
        }else{
            print("input data must be either a matrix, a object of class im or a model fit")
        }
    }
    out <- rpoispp(im)
    tod <- sample(1:23, out$n,replace = TRUE, prob = tod.prob)
    dow <- sample(1:7, out$n,replace = TRUE, prob = dow.prob)
    c <- sample(1:20, out$n,replace = TRUE)
    out <- cbind(out$x, out$y,tod,dow,c)
    colnames(out) <- c("x","y","Hour","Day","craving")
    out
}
    
