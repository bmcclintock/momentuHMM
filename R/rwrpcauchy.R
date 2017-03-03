#' @importFrom stats rcauchy
# wrapped cauchy distribution with rho in (-1,1)
rwrpcauchy<-function (n, location = 0, rho = exp(-1)) 
{
  if(rho == 0) 
    result <- runif(n, 0, 2 * pi)
  else if(abs(rho)==1) 
    result <- rep(location - pi * (rho<0), n)%%(2*pi)
  else {
    if(rho < 0){ 
      rho <- abs(rho)
      location <- location - pi
    }
    scale <- -log(rho)
    result <- rcauchy(n, location, scale)%%(2 * pi)
  }
  result
}