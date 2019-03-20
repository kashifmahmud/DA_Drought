# # Setting lower and upper bounds of the prior parameter pdf, and starting point of the chain
# Try wide priors
#-------------------------------------------------------------------------------------
# param.Y <- matrix(c(0.2,0.3,0.4) , nrow=1, ncol=3, byrow=T)
# param.af <- matrix(c(0,0.25,1) , nrow=1, ncol=3, byrow=T) # Forcing af not to be negetive
# param.as <- matrix(c(0,0.65,1) , nrow=1, ncol=3, byrow=T)
# param.sf <- matrix(c(0,0.0005,0.0015) , nrow=1, ncol=3, byrow=T) # All Groups having same sf
# param.sr <- matrix(c(0,0.0001,0.0002) , nrow=1, ncol=3, byrow=T) # All Groups having same sr

param.Y <- matrix(c(0.2,0.3,0.4) , nrow=1, ncol=3, byrow=T)
param.af <- matrix(c(0,0.25,0.4) , nrow=1, ncol=3, byrow=T)
param.as <- matrix(c(0.4,0.6,0.8) , nrow=1, ncol=3, byrow=T)
param.sf <- matrix(c(0,0.001,0.002) , nrow=1, ncol=3, byrow=T)
param.sr <- matrix(c(0,0.01,0.02) , nrow=1, ncol=3, byrow=T)

# # wide ranges (didn't work)
# param.Y = matrix(c(0.2,0.3,0.4) , nrow=1, ncol=3, byrow=T)
# param.af = matrix(c(0,0.2,1) , nrow=1, ncol=3, byrow=T)
# param.as = matrix(c(0,0.5,1) , nrow=1, ncol=3, byrow=T)
# param.sf = matrix(c(0,0.01,0.02) , nrow=1, ncol=3, byrow=T)
# param.sr = matrix(c(0,0.01,0.02) , nrow=1, ncol=3, byrow=T)

# # initialize 'sf' prior differently for grouped treatments
# if (v[[1]]==1 || v[[1]]==2 || v[[1]]==3) {
#   param.sf <- matrix(c(0,0.01,0.02) , nrow=1, ncol=3, byrow=T) # Group 1
# } else if (v[[1]]==4 || v[[1]]==5 || v[[1]]==6) {
#   param.sf <- matrix(c(0,0.005,0.01) , nrow=1, ncol=3, byrow=T) # Group 2
# } else if (v[[1]]==7) {
#   param.sf <- matrix(c(0,0.0025,0.005) , nrow=1, ncol=3, byrow=T) # Free seedling
# }

# # initialize 'sf' prior differently for grouped treatments
# if (v[[1]]==1 || v[[1]]==2 || v[[1]]==3) {
#   param.sf <- matrix(c(0,0.0125,0.025) , nrow=1, ncol=3, byrow=T) # Group 1
# } else if (v[[1]]==4 || v[[1]]==5 || v[[1]]==6) {
#   param.sf <- matrix(c(0,0.0075,0.015) , nrow=1, ncol=3, byrow=T) # Group 2
# } else if (v[[1]]==7) {
#   param.sf <- matrix(c(0,0.005,0.01) , nrow=1, ncol=3, byrow=T) # Free seedling
# }

# if (length(vol.group)==1) {
#   param.sf <- matrix(c(0,0.01,0.02) , nrow=1, ncol=3, byrow=T) # 1 Group
# } else if (length(vol.group)==2) {
#   if (v1==1) {
#     param.sf <- matrix(c(0,0.01,0.02) , nrow=1, ncol=3, byrow=T) # Group 1
#   } else if (v1==2) {
#     param.sf <- matrix(c(0,0.0025,0.005) , nrow=1, ncol=3, byrow=T) # Group 2
#   }
# } else if (length(vol.group)==3) {
#   if (v1==1) {
#     param.sf <- matrix(c(0,0.01,0.02) , nrow=1, ncol=3, byrow=T) # Group 1
#   } else if (v1==2) {
#     param.sf <- matrix(c(0,0.005,0.01) , nrow=1, ncol=3, byrow=T) # Group 2
#   } else if (v1==3) {
#     param.sf <- matrix(c(0,0.0025,0.005) , nrow=1, ncol=3, byrow=T) # Free seedling
#   }
# } else if (length(vol.group)>3) {
#   if (v[[1]]==1 || v[[1]]==2 || v[[1]]==3) {
#     param.sf <- matrix(c(0,0.01,0.02) , nrow=1, ncol=3, byrow=T) # Group 1
#   } else if (v[[1]]==4 || v[[1]]==5 || v[[1]]==6) {
#     param.sf <- matrix(c(0,0.005,0.01) , nrow=1, ncol=3, byrow=T) # Group 2
#   } else if (v[[1]]==7) {
#     param.sf <- matrix(c(0,0.0025,0.005) , nrow=1, ncol=3, byrow=T) # Free seedling
#   }
# }

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
if (no.param > 1) {
  param.Y <- rbind(param.Y, c(-(param.Y[3]-param.Y[1])/2/nrow(data), 0, (param.Y[3]-param.Y[1])/2/nrow(data)))
  param.af <- rbind(param.af, c(-(param.af[3]-param.af[1])/nrow(data), 0, (param.af[3]-param.af[1])/nrow(data)))
  param.as <- rbind(param.as, c(-(param.as[3]-param.as[1])/2/nrow(data), 0, (param.as[3]-param.as[1])/2/nrow(data)))
  param.sf <- rbind(param.sf, c(-(param.sf[3]-param.sf[1])/nrow(data), 0, (param.sf[3]-param.sf[1])/nrow(data)))
  param.sr <- rbind(param.sr, c(-(param.sr[3]-param.sr[1])/nrow(data), 0, (param.sr[3]-param.sr[1])/nrow(data)))
}

if (no.param > 2) {
  param.Y <- rbind(param.Y, c((param.Y[1,1]-param.Y[1,3]-param.Y[2,3]*nrow(data))/(2*nrow(data)^2), 0, (param.Y[1,3]-param.Y[1,1]-param.Y[2,1]*nrow(data))/(2*nrow(data)^2)))
  param.af <- rbind(param.af, c((param.af[1,1]-param.af[1,3]-param.af[2,3]*nrow(data))/(2*nrow(data)^2), 0, (param.af[1,3]-param.af[1,1]-param.af[2,1]*nrow(data))/(nrow(data)^2)))
  param.as <- rbind(param.as, c((param.as[1,1]-param.as[1,3]-param.as[2,3]*nrow(data))/(2*nrow(data)^2), 0, (param.as[1,3]-param.as[1,1]-param.as[2,1]*nrow(data))/(2*nrow(data)^2)))
  param.sf <- rbind(param.sf, c((param.sf[1,1]-param.sf[1,3]-param.sf[2,3]*nrow(data))/(2*nrow(data)^2), 0, (param.sf[1,3]-param.sf[1,1]-param.sf[2,1]*nrow(data))/(nrow(data)^2)))
  param.sr <- rbind(param.sr, c((param.sr[1,1]-param.sr[1,3]-param.sr[2,3]*nrow(data))/(2*nrow(data)^2), 0, (param.sr[1,3]-param.sr[1,1]-param.sr[2,1]*nrow(data))/(nrow(data)^2)))
}
if (no.param > 3) {
  param.Y <- rbind(param.Y, c((param.Y[1,1]-param.Y[1,3]-param.Y[2,3]*nrow(data)-param.Y[3,3]*(nrow(data)^2))/(2*nrow(data)^3), 0, (param.Y[1,3]-param.Y[1,1]-param.Y[2,1]*nrow(data)-param.Y[3,1]*(nrow(data)^2))/(2*nrow(data)^3)))
  param.af <- rbind(param.af, c((param.af[1,1]-param.af[1,3]-param.af[2,3]*nrow(data)-param.af[3,3]*(nrow(data)^2))/(2*nrow(data)^3), 0, (param.af[1,3]-param.af[1,1]-param.af[2,1]*nrow(data)-param.af[3,1]*(nrow(data)^2))/(2*nrow(data)^3)))
  param.as <- rbind(param.as, c((param.as[1,1]-param.as[1,3]-param.as[2,3]*nrow(data)-param.as[3,3]*(nrow(data)^2))/(2*nrow(data)^3), 0, (param.as[1,3]-param.as[1,1]-param.as[2,1]*nrow(data)-param.as[3,1]*(nrow(data)^2))/(2*nrow(data)^3)))
  param.sf <- rbind(param.sf, c((param.sf[1,1]-param.sf[1,3]-param.sf[2,3]*nrow(data)-param.sf[3,3]*(nrow(data)^2))/(2*nrow(data)^3), 0, (param.sf[1,3]-param.sf[1,1]-param.sf[2,1]*nrow(data)-param.sf[3,1]*(nrow(data)^2))/(2*nrow(data)^3)))
  param.sr <- rbind(param.sr, c((param.sr[1,1]-param.sr[1,3]-param.sr[2,3]*nrow(data)-param.sr[3,3]*(nrow(data)^2))/(2*nrow(data)^3), 0, (param.sr[1,3]-param.sr[1,1]-param.sr[2,1]*nrow(data)-param.sr[3,1]*(nrow(data)^2))/(2*nrow(data)^3)))
}
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
if (with.storage==F) {
  param = data.frame(param.Y,param.af,param.as,param.sf,param.sr)
  names(param) <- c("Y_min","Y","Y_max","af_min","af","af_max","as_min","as","as_max","sf_min","sf","sf_max","sr_min","sr","sr_max")
  pMinima <- param[ ,c("Y_min","af_min","as_min","sf_min","sr_min")]
  pMaxima <- param[ ,c("Y_max","af_max","as_max","sf_max","sr_max")]
  pValues <- param[ ,c("Y","af","as","sf","sr")] # Starting point of the chain
} else { # (with.storage==T) 
  # param.k <- matrix(c(0,0.25,1) , nrow=1, ncol=3, byrow=T)
  param.k <- matrix(c(0,0.15,0.3) , nrow=1, ncol=3, byrow=T)
  
  if (no.param > 1) {
    param.k <- rbind(param.k, c(-(param.k[3]-param.k[1])/nrow(data), 0, (param.k[3]-param.k[1])/nrow(data)))
  } 
  if (no.param > 2) {
    param.k <- rbind(param.k, c((param.k[1,1]-param.k[1,3]-param.k[2,3]*nrow(data))/(nrow(data)^2), 0, (param.k[1,3]-param.k[1,1]-param.k[2,1]*nrow(data))/(nrow(data)^2)))
  }
  if (no.param > 3) {
    param.k <- rbind(param.k, c((param.k[1,1]-param.k[1,3]-param.k[2,3]*nrow(data)-param.k[3,3]*(nrow(data)^2))/(2*nrow(data)^3), 0, (param.k[1,3]-param.k[1,1]-param.k[2,1]*nrow(data)-param.k[3,1]*(nrow(data)^2))/(2*nrow(data)^3)))
  }
  param = data.frame(param.k,param.Y,param.af,param.as,param.sf,param.sr)
  names(param) <- c("k_min","k","k_max","Y_min","Y","Y_max","af_min","af","af_max","as_min","as","as_max","sf_min","sf","sf_max","sr_min","sr","sr_max")
  pMinima <- param[ ,c("k_min","Y_min","af_min","as_min","sf_min","sr_min")]
  pMaxima <- param[ ,c("k_max","Y_max","af_max","as_max","sf_max","sr_max")]
  pValues <- param[ ,c("k","Y","af","as","sf","sr")] # Starting point of the chain
}
pChain <- matrix(0, nrow=chainLength, ncol = no.param*no.var+1) # Initialising the chain
#-------------------------------------------------------------------------------------

