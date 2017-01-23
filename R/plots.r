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
#'
#' @export
#'
plot.hawkes <- function(times = NULL, mu = NULL, alpha = NULL, beta = NULL){
    n <- length(times)
    max <- max(times)
    p <- seq(0,max,length.out=500)
    lam.p <- hawke.intensity(mu = mu, alpha = alpha, beta = beta, times = times, p = p)
    lmax <- max(lam.p)
    plot(times,rep(1,n),ylim = c(-0.5,lmax),xlab="time",ylab = expression(lambda(t)))
    lines(p,lam.p,col=2)
}
