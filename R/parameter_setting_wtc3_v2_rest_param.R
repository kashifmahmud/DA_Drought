# # Setting lower and upper bounds of the prior parameter pdf, and starting point of the chain
# Try wide priors
#-------------------------------------------------------------------------------------
param.af <- matrix(c(0,0.25,0.5) , nrow=1, ncol=3, byrow=T)
param.as <- matrix(c(0.3,0.65,0.9) , nrow=1, ncol=3, byrow=T)
param.sf <- matrix(c(0,0.001,0.002) , nrow=1, ncol=3, byrow=T)
param.sr <- matrix(c(0,0.01,0.02) , nrow=1, ncol=3, byrow=T)

#-------------------------------------------------------------------------------------
if (no.param > 1) {
  param.af <- rbind(param.af, c(-(param.af[3]-param.af[1])/nrow(data), 0, (param.af[3]-param.af[1])/nrow(data)))
  param.as <- rbind(param.as, c(-(param.as[3]-param.as[1])/2/nrow(data), 0, (param.as[3]-param.as[1])/2/nrow(data)))
  param.sf <- rbind(param.sf, c(-(param.sf[3]-param.sf[1])/nrow(data), 0, (param.sf[3]-param.sf[1])/nrow(data)))
  param.sr <- rbind(param.sr, c(-(param.sr[3]-param.sr[1])/nrow(data), 0, (param.sr[3]-param.sr[1])/nrow(data)))
}

if (no.param > 2) {
  param.af <- rbind(param.af, c((param.af[1,1]-param.af[1,3]-param.af[2,3]*nrow(data))/(2*nrow(data)^2), 0, (param.af[1,3]-param.af[1,1]-param.af[2,1]*nrow(data))/(nrow(data)^2)))
  param.as <- rbind(param.as, c((param.as[1,1]-param.as[1,3]-param.as[2,3]*nrow(data))/(2*nrow(data)^2), 0, (param.as[1,3]-param.as[1,1]-param.as[2,1]*nrow(data))/(2*nrow(data)^2)))
  param.sf <- rbind(param.sf, c((param.sf[1,1]-param.sf[1,3]-param.sf[2,3]*nrow(data))/(2*nrow(data)^2), 0, (param.sf[1,3]-param.sf[1,1]-param.sf[2,1]*nrow(data))/(nrow(data)^2)))
  param.sr <- rbind(param.sr, c((param.sr[1,1]-param.sr[1,3]-param.sr[2,3]*nrow(data))/(2*nrow(data)^2), 0, (param.sr[1,3]-param.sr[1,1]-param.sr[2,1]*nrow(data))/(nrow(data)^2)))
}
if (no.param > 3) {
  param.af <- rbind(param.af, c((param.af[1,1]-param.af[1,3]-param.af[2,3]*nrow(data)-param.af[3,3]*(nrow(data)^2))/(2*nrow(data)^3), 0, (param.af[1,3]-param.af[1,1]-param.af[2,1]*nrow(data)-param.af[3,1]*(nrow(data)^2))/(2*nrow(data)^3)))
  param.as <- rbind(param.as, c((param.as[1,1]-param.as[1,3]-param.as[2,3]*nrow(data)-param.as[3,3]*(nrow(data)^2))/(2*nrow(data)^3), 0, (param.as[1,3]-param.as[1,1]-param.as[2,1]*nrow(data)-param.as[3,1]*(nrow(data)^2))/(2*nrow(data)^3)))
  param.sf <- rbind(param.sf, c((param.sf[1,1]-param.sf[1,3]-param.sf[2,3]*nrow(data)-param.sf[3,3]*(nrow(data)^2))/(2*nrow(data)^3), 0, (param.sf[1,3]-param.sf[1,1]-param.sf[2,1]*nrow(data)-param.sf[3,1]*(nrow(data)^2))/(2*nrow(data)^3)))
  param.sr <- rbind(param.sr, c((param.sr[1,1]-param.sr[1,3]-param.sr[2,3]*nrow(data)-param.sr[3,3]*(nrow(data)^2))/(2*nrow(data)^3), 0, (param.sr[1,3]-param.sr[1,1]-param.sr[2,1]*nrow(data)-param.sr[3,1]*(nrow(data)^2))/(2*nrow(data)^3)))
}
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
param.2 = data.frame(param.af,param.as,param.sf,param.sr)
names(param.2) <- c("af_min","af","af_max","as_min","as","as_max","sf_min","sf","sf_max","sr_min","sr","sr_max")
pMinima.2 <- param.2[ ,c("af_min","as_min","sf_min","sr_min")]
pMaxima.2 <- param.2[ ,c("af_max","as_max","sf_max","sr_max")]
pValues.2 <- param.2[ ,c("af","as","sf","sr")] # Starting point of the chain

pChain <- matrix(0, nrow=chainLength, ncol = no.param*no.var+1) # Initialising the chain
#-------------------------------------------------------------------------------------

