# # Setting lower and upper bounds of the prior parameter pdf, and starting point of the chain
# Try wide priors
#-------------------------------------------------------------------------------------
param.Y <- matrix(c(0.2,0.3,0.4) , nrow=1, ncol=3, byrow=T)

#-------------------------------------------------------------------------------------
if (no.param > 1) {
    param.Y <- rbind(param.Y, c(-(param.Y[3]-param.Y[1])/2/nrow(data), 0, (param.Y[3]-param.Y[1])/2/nrow(data)))
  }
  
if (no.param > 2) {
    param.Y <- rbind(param.Y, c((param.Y[1,1]-param.Y[1,3]-param.Y[2,3]*nrow(data))/(2*nrow(data)^2), 0, (param.Y[1,3]-param.Y[1,1]-param.Y[2,1]*nrow(data))/(2*nrow(data)^2)))
  }
  
if (no.param > 3) {
    param.Y <- rbind(param.Y, c((param.Y[1,1]-param.Y[1,3]-param.Y[2,3]*nrow(data)-param.Y[3,3]*(nrow(data)^2))/(2*nrow(data)^3), 0, (param.Y[1,3]-param.Y[1,1]-param.Y[2,1]*nrow(data)-param.Y[3,1]*(nrow(data)^2))/(2*nrow(data)^3)))
  }
  
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
if (with.storage==F) {
  param.1 = data.frame(param.Y)
  names(param) <- c("Y_min","Y","Y_max")
  pMinima <- param[ ,c("Y_min")]
  pMaxima <- param[ ,c("Y_max")]
  pValues <- param[ ,c("Y")] # Starting point of the chain
} else { # (with.storage==T) 
  param.k <- matrix(c(0,0.15,0.25) , nrow=1, ncol=3, byrow=T)
    if (no.param > 1) {
      param.k <- rbind(param.k, c(-(param.k[3]-param.k[1])/nrow(data), 0, (param.k[3]-param.k[1])/nrow(data)))
    }
    if (no.param > 2) {
      param.k <- rbind(param.k, c((param.k[1,1]-param.k[1,3]-param.k[2,3]*nrow(data))/(nrow(data)^2), 0, (param.k[1,3]-param.k[1,1]-param.k[2,1]*nrow(data))/(nrow(data)^2)))
    }
    if (no.param > 3) {
      param.k <- rbind(param.k, c((param.k[1,1]-param.k[1,3]-param.k[2,3]*nrow(data)-param.k[3,3]*(nrow(data)^2))/(2*nrow(data)^3), 0, (param.k[1,3]-param.k[1,1]-param.k[2,1]*nrow(data)-param.k[3,1]*(nrow(data)^2))/(2*nrow(data)^3)))
    }
  param.1 = data.frame(param.k,param.Y)
  param.1 = rbind(param.1,param.1)
  names(param.1) <- c("k_min","k","k_max","Y_min","Y","Y_max")
  pMinima.1 <- param.1[ ,c("k_min","Y_min")]
  pMaxima.1 <- param.1[ ,c("k_max","Y_max")]
  pValues.1 <- param.1[ ,c("k","Y")] # Starting point of the chain
}
# pChain <- matrix(0, nrow=chainLength, ncol = no.param*no.var+1) # Initialising the chain
#-------------------------------------------------------------------------------------

