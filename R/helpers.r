## simple extraction of scottish county

county <- function(sp = NULL, name = NULL){
    county <- sp[sp$NAME == name, ]
    county
}
    
