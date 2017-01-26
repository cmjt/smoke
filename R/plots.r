#' Plot point process data on Google earth
#'
#' @return an explorartory plot of the supplied data. markers will popup with factor values of each location which
#' can also be subset to plot only a single factor value.
#' 
#' @param x  a matrix of three named colums Longitude, Latitude, factor respectively,
#' 
#' @param id an optional argument specifying a factor of interest supplied to \code{locs} to plot
#'
#' @export

plot.smoke <- function(x = NULL, id = NULL){
    if(!is.null(id)){
        id <- id
        x <- subset(x, x$factor == id)
    }
    locs <- cbind(x[,'Longitude'],x[,'Latitude'])
    x <- SpatialPointsDataFrame(locs
                               ,data.frame(type = x[,'factor']))
    if(length(table(x$type)) == 2){ # if there are two factors plot a legend and corresponding coloured dots
        pal <- colorFactor(c( "red","navy"), domain = c("0", "1"))
        map <- leaflet(x) %>% addTiles() %>%
            setView(mean(locs[,1]),mean(locs[,2]), zoom = 4) %>%
            addCircleMarkers(data = x, color = ~pal(x$type),
                             radius = ~ifelse(type=="0", 3, 3),
                             stroke = TRUE, fillOpacity = 0.5) %>%
            addLegend(pal = pal, values = na.omit(x$type), opacity = 1,title="Factor",na.label="")
    }else{
        map <- leaflet(x) %>% addTiles() %>%
            setView(mean(locs[,1]),mean(locs[,2]), zoom = 4) %>%
            addCircleMarkers(data = x,stroke = FALSE, fillOpacity = 0.5,popup = as.character(x$type),
                             clusterOptions = markerClusterOptions())
    }
    map
}
        


#' Plotting simple intensity of a hawkes process
#' that is where intensity  = mu + sum(alpha * b(t))
#' where b(t) = beta*exp(-beta*t)
#' @param times numeric vector of univariate point locations
#' @param mu numeric immigrant intensity
#' @param alpha numeric offsping intensity
#' @param beta numeric normalised offspring intensity
#' @param n.s optional numeric value giving the extected number of extraneous exents
#' @param gamma optional numeric value giving the jump in intensity post an extraneous event s. must
#' @param mark a named list to be supplied if the process is marked where rge marks follow an exponential
#' distribution. list must contain \code{delta} a numeric rate parameter of mark distribution, \code{nu}
#' the exponential growth of the marks contribution to the self excitement of the process, \code{marks}
#' a numeric vector of marks at each time in \code{time}. NOTE marks are assumed to follow an exponential distribution.
#' be supplied if n.s is
#'
#' @export
#'
plot.hawkes <- function(times = NULL, mu = NULL, alpha = NULL, beta = NULL,n.s = NULL, gamma = NULL, mark = NULL){
    n <- length(times)
    max <- max(times)
    p <- seq(0,max,length.out=500)
    if(!is.null(n.s)){
        lam.p <- hawke.intensity(mu = mu, alpha = alpha, beta = beta, times = times,
                                 p = p, n.s = n.s, gamma = gamma)
        ylab <- expression(lambda(t))
    }
    if(!is.null(mark)){
        nu <- mark[["nu"]]
        delta <- mark[["delta"]]
        marks <- mark[["marks"]]
        k <- seq(min(marks),max(marks),length.out=500)
        lam.p <- hawke.mark.intensity(mu = mu, alpha = alpha, beta = beta, nu = nu,
                                      delta = delta, times = times, marks = marks, p = p, k = k)
        ylab <- expression(lambda(t,m))
        }else{
            lam.p <- hawke.intensity(mu = mu, alpha = alpha, beta = beta, times = times, p = p)
            ylab <- expression(lambda(t))
    }
    lmax <- max(lam.p)
    lmin <- min(lam.p)
    plot(times,rep(lmin-0.5,n),ylim = c(lmin-1.5,lmax),xlab="time",ylab = ylab)
    lines(p,lam.p,col="grey")
}
