moveData<-function (data) 
{
    if (is.null(data$ID) | is.null(data$step) | is.null(data$x) | 
        is.null(data$y)) 
        stop("Can't construct moveData object: fields are missing")
    obj <- data
    class(obj) <- append("moveData", class(obj))
    return(obj)
}