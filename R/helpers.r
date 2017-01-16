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
        
    
   
