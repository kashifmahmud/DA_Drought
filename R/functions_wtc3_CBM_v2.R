#----------------------------------------------------------------------------------------------------------------
#- This is a collection of many functions that run the CBM. These functions
#    are called by just a few lines of code in "CentralScript.R" to recreate the analyses and figures.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

# This script calcualtes LogLikelihood to find the most accurate model
logLikelihood.wtc3 <- function (no.param.par.var,data.set,output,with.storage,model.comparison) {
  logLi <- matrix(0, nrow=nrow(data.set), ncol = 1) # Initialising the logLi
  # data_count = sum(!is.na(data.set$LM)) + sum(!is.na(data.set$WM)) + sum(!is.na(data.set$RM)) + sum(!is.na(data.set$Ra)) + sum(!is.na(data.set$TNC_leaf)) + sum(!is.na(data.set$litter))
  data_count = sum(!is.na(data.set$LM)) + sum(!is.na(data.set$WM)) + sum(!is.na(data.set$RM)) + sum(!is.na(data.set$TNC_leaf)) + sum(!is.na(data.set$litter))
  # data_count = sum(!is.na(data.set$LM)) + sum(!is.na(data.set$WM)) + sum(!is.na(data.set$RM)) + sum(!is.na(data.set$litter))
  
  for (i in 1:nrow(data.set)) {
    if (!is.na(data.set$LM[i])) {
      logLi[i] = - 2*((data_count/sum(!is.na(data.set$LM)))*(0.5*((output$Mleaf[i] - data.set$LM[i])/data.set$LM_SE[i])^2 - log(data.set$LM_SE[i]) - log(2*pi)^0.5))
      # logLi[i] = - 0.5*((output$Mleaf[i] - data.set$LM[i])/data.set$LM_SE[i])^2 - log(data.set$LM_SE[i]) - log(2*pi)^0.5
    }
    if (!is.na(data.set$WM[i])) {
      logLi[i] = logLi[i] - ((data_count/sum(!is.na(data.set$WM)))*(0.5*((output$Mwood[i] - data.set$WM[i])/data.set$WM_SE[i])^2 - log(data.set$WM_SE[i]) - log(2*pi)^0.5))
      # logLi[i] = logLi[i] - 0.5*((output$Mwood[i] - data.set$WM[i])/data.set$WM_SE[i])^2 - log(data.set$WM_SE[i]) - log(2*pi)^0.5
    }
    if (!is.na(data.set$RM[i])) {
      logLi[i] = logLi[i] - ((data_count/sum(!is.na(data.set$RM)))*(0.5*((output$Mroot[i] - data.set$RM[i])/data.set$RM_SE[i])^2 - log(data.set$RM_SE[i]) - log(2*pi)^0.5)) # multiplied by 2 to give extra weight
      # logLi[i] = logLi[i] - ((output$Mroot[i] - data.set$RM[i])/data.set$RM_SE[i])^2 - log(data.set$RM_SE[i]) - log(2*pi)^0.5 # multiplied by 20 to give extra weight
    }
    # if (i > 1) {
    #   if (!is.na(data.set$Ra[i])) {
    #     logLi[i] = logLi[i] - ((data_count/sum(!is.na(data.set$Ra)))*(0.5*((output$Rabove[i] - data.set$Ra[i]) / data.set$Ra_SE[i])^2 - log(data.set$Ra_SE[i]) - log(2*pi)^0.5))
    #   }
    # }
    if (!is.null(data.set$litter)) {
      if (!is.na(data.set$litter[i])) {
        logLi[i] = logLi[i] - ((data_count/sum(!is.na(data.set$litter)))*(0.5*((output$Mlit[i] - data.set$litter[i])/data.set$litter_SE[i])^2 - log(data.set$litter_SE[i]) - log(2*pi)^0.5))
        # logLi[i] = logLi[i] - 0.5*((output$Mlit[i] - data.set$litter[i])/data.set$litter_SE[i])^2 - log(data.set$litter_SE[i]) - log(2*pi)^0.5
      }
    }
    # if (model.comparison==F) {
    if (with.storage==T) {
      if (!is.null(data.set$TNC_leaf)) {
        if (!is.na(data.set$TNC_leaf[i])) {
          logLi[i] = logLi[i] - ((data_count/sum(!is.na(data.set$TNC_leaf)))*(0.5*((output$Sleaf[i] - data.set$TNC_leaf[i])/data.set$TNC_leaf_SE[i])^2 - log(data.set$TNC_leaf_SE[i]) - log(2*pi)^0.5))
          # logLi[i] = logLi[i] - (0.5*((output$Sleaf[i] - data.set$TNC_leaf[i])/data.set$TNC_leaf_SE[i])^2 - log(data.set$TNC_leaf_SE[i]) - log(2*pi)^0.5)
        }
      }
    }
    # }
  }
  return(sum(logLi)/data_count)
  # return(sum(logLi))
}

# This script calcualtes LogLikelihood to find the most accurate model
logLikelihood.wtc3.final <- function (no.param.par.var,data.set,output,with.storage,model.comparison) {
  logLi <- matrix(0, nrow=nrow(data.set), ncol = 1) # Initialising the logLi
  
  for (i in 1:nrow(data.set)) {
    if (!is.na(data.set$LM[i])) {
      logLi[i] = - 0.5*((output$Mleaf[i] - data.set$LM[i])/data.set$LM_SE[i])^2 - log(data.set$LM_SE[i]) - log(2*pi)^0.5
    }
    if (!is.na(data.set$WM[i])) {
      logLi[i] = logLi[i] - 0.5*((output$Mwood[i] - data.set$WM[i])/data.set$WM_SE[i])^2 - log(data.set$WM_SE[i]) - log(2*pi)^0.5
    }
    if (!is.na(data.set$RM[i])) {
      # logLi[i] = logLi[i] - 10*(0.5*((output$Mroot[i] - data.set$RM[i])/data.set$RM_SE[i])^2 - log(data.set$RM_SE[i]) - log(2*pi)^0.5) # multiplied by 20 to give extra weight
      logLi[i] = logLi[i] - 0.5*((output$Mroot[i] - data.set$RM[i])/data.set$RM_SE[i])^2 - log(data.set$RM_SE[i]) - log(2*pi)^0.5 # multiplied by 20 to give extra weight
    }
    # if (i > 1) {
    #   if (!is.na(data.set$Ra[i])) {
    #     # if (no.param.par.var < 4) {
    #     # logLi[i] = logLi[i] - 0.5*(((data.set$Rd.foliage.mean[i]*output$Mleaf[i] + data.set$Rd.stem.mean[i]*output$Mwood[i]*data.set$SMratio[i] + data.set$Rd.branch.mean[i]*output$Mwood[i]*data.set$BMratio[i] +
    #     #                                Y[i]*(output$Mleaf[i]-output$Mleaf[i-1] + output$Mwood[i]-output$Mwood[i-1]) ) -
    #     #                               data.set$Ra[i]) / data.set$Ra_SE[i])^2 - log(data.set$Ra_SE[i]) - log(2*pi)^0.5
    #     logLi[i] = logLi[i] - (0.5*((output$Rabove[i] - data.set$Ra[i]) / data.set$Ra_SE[i])^2 - log(data.set$Ra_SE[i]) - log(2*pi)^0.5)
    #     # logLi[i] = logLi[i] - 0.1*(0.5*((output$Rabove[i] - data.set$Ra[i]) / data.set$Ra_SE[i])^2 - log(data.set$Ra_SE[i]) - log(2*pi)^0.5)
    #     # } else { # no.param.par.var > 4; monthly parameter setting)
    #     # logLi[i] = logLi[i] - 0.5*(((data.set$Rd.foliage.mean[i]*output$Mleaf[i] + data.set$Rd.stem.mean[i]*output$Mwood[i]*data.set$SMratio[i] + data.set$Rd.branch.mean[i]*output$Mwood[i]*data.set$BMratio[i] +
    #     #                                Y[(i-1)-(j[i-1])]*(output$Mleaf[i]-output$Mleaf[i-1] + output$Mwood[i]-output$Mwood[i-1]) ) -
    #     #                               data.set$Ra[i]) / data.set$Ra_SE[i])^2 - log(data.set$Ra_SE[i]) - log(2*pi)^0.5
    #     # logLi[i] = logLi[i] - (0.5*((output$Rabove[i] - data.set$Ra[i]) / data.set$Ra_SE[i])^2 - log(data.set$Ra_SE[i]) - log(2*pi)^0.5)
    #     # logLi[i] = logLi[i] - 0.1*(0.5*((output$Rabove[i] - data.set$Ra[i]) / data.set$Ra_SE[i])^2 - log(data.set$Ra_SE[i]) - log(2*pi)^0.5)
    #     # }
    #   }
    # }
    if (!is.null(data.set$litter)) {
      if (!is.na(data.set$litter[i])) {
        logLi[i] = logLi[i] - 0.5*((output$Mlit[i] - data.set$litter[i])/data.set$litter_SE[i])^2 - log(data.set$litter_SE[i]) - log(2*pi)^0.5
      }
    }
    
    if (model.comparison==F) {
      if (with.storage==T) {
        if (!is.null(data.set$TNC_leaf)) {
          if (!is.na(data.set$TNC_leaf[i])) {
            # logLi[i] = logLi[i] - 1.5*(0.5*((output$Sleaf[i] - data.set$TNC_leaf[i])/data.set$TNC_leaf_SE[i])^2 - log(data.set$TNC_leaf_SE[i]) - log(2*pi)^0.5)
            logLi[i] = logLi[i] - (0.5*((output$Sleaf[i] - data.set$TNC_leaf[i])/data.set$TNC_leaf_SE[i])^2 - log(data.set$TNC_leaf_SE[i]) - log(2*pi)^0.5)
          }
        }
      }
    }
  }
  return(sum(logLi))
}
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Carbon balance model (CBM)
# Developed by Kashif Mahmud and Belinda Medlyn (March 2017)
# k.mahmud@westernsydney.edu.au

# This version tries to group various treatments according to their similarities to have a trend in paramter settings

# This code carries out Bayesian calibration for 5 variables (allocation fractions: "k","Y",af","as","sf") on 
# various temporal scales (e.g. 1,2,...,121 days) to estimate Carbon pools (Cstorage,Cleaf,Cwood,Croot) and fluxes
#-------------------------------------------------------------------------------------
CBM.wtc3 <- function(chainLength, no.param.par.var, treat.group, with.storage, model.comparison, model.optimization) {
  
  source("R/load_packages_wtc3.R")
  
  # Assign inputs for MCMC
  bunr_in = chainLength * 0.1 # Discard the first 10% iterations for Burn-IN in MCMC (According to Oijen, 2008)
  if (with.storage==T) {
    no.var = 6 # variables to be modelled are: k,Y,af,as,sf,sr
  } else {
    no.var = 5 # variables to be modelled are: Y,af,as,sf,sr
  }
  
  treat = as.factor(c("ambient","elevated","amb+ele","amb+ele+const")) # Assign all treatments
  treat = factor(treat, levels=c(sort(setdiff(unique(treat), 'amb+ele'),decreasing=FALSE), 'amb+ele'))
  
  param.mean = data.frame(matrix(ncol = no.var+1, nrow = length(no.param.par.var)*length(treat.group)))
  if (with.storage==T) {
    names(param.mean) = c("k","Y","af","as","ar","sf","sr")
  } else {
    names(param.mean) = c("Y","af","as","ar","sf","sr")
  }
  aic.bic = data.frame(matrix(ncol = 6, nrow = length(no.param.par.var)*length(treat.group)))
  names(aic.bic) <- c("logLi","aic","bic","time","no.param","treatment")
  time = data.frame(no.param=rep(no.param.par.var,length(treat.group)),
                    start.time=numeric(length(no.param.par.var)*length(treat.group)),
                    end.time=numeric(length(no.param.par.var)*length(treat.group)),
                    time.taken=numeric(length(no.param.par.var)*length(treat.group)))
  
  q = 0 # Indicates the iteration number
  # set.seed(3)
  # set.seed(15) # final seed for reproducible results
  # set.seed(18) 
  
  # Start the iteration for different treatment group and number of parameters
  for (v1 in 1:length(treat.group)) {
    v = unlist(treat.group[v1])
    # # This script take the subset of processed data for particular treatment group
    # source("R/data_processing_wtc3.R", local=TRUE)
    data = subset(data.all,(treatment.no %in% treat.group[v]))
    
    for (z in 1:length(no.param.par.var)) {
      # Initialize few output data files
      q = q + 1
      time$start.time[q] <- Sys.time()
      # param.vary = 31 # monthly = 31
      param.vary = ceiling(nrow(data)/no.param.par.var[z]) # How many days the parameter set remain unchanged (weekly = 8; monthly = 31; just one parameter = nrow(data))
      # no.param = ceiling(nrow(data.set)/param.vary) # number of parameter set for the whole duration of experiment (121 days)
      no.param = ceiling(nrow(data)/param.vary) # number of parameter set for the whole duration of experiment (121 days)
      
      if (no.param.par.var < 5) {
        # This script initializes the parameter setting
        if (!treat[v1] %in% as.factor("amb+ele+const")) {
          source("R/parameter_setting_wtc3.R", local=TRUE) # initialize 'sf' prior differently for grouped treatments
        } else { # (treat[v1] %in% as.factor("amb+ele+const")) 
          source("R/parameter_setting_wtc3_v2_k_Y.R", local=TRUE) # k, Y - Dont vary; Rest of the parameter - Vary over time
          for (j in 1:length(v)) {
            source("R/parameter_setting_wtc3_v2_rest_param.R", local=TRUE)
            if (j == 1) {
              output = output.set
            }
            if (j > 1) {
              output = rbind(output,output.set)
            }
          }
        }
      } else { # no.param.par.var > 5; monthly parameter setting) 
        j = c()
        j[1] = 0
        i = seq(1,nrow(data),1)
        j[i] = i - ceiling(i/param.vary)*1  # j is for parameter settings for various time frames
        
        # Setting lower and upper bounds of the prior parameter pdf, and starting point of the chain
        if (with.storage==T) {
          param.k <- matrix(c(0,0.15,0.25) , nrow=no.param, ncol=3, byrow=T)
          # param.k <- matrix(c(0,0.15,1) , nrow=no.param, ncol=3, byrow=T)
        }
        param.Y <- matrix(c(0.2,0.3,0.4) , nrow=no.param, ncol=3, byrow=T)
        param.af <- matrix(c(0,0.25,0.5) , nrow=no.param, ncol=3, byrow=T)
        param.as <- matrix(c(0.3,0.65,0.9) , nrow=no.param, ncol=3, byrow=T)
        param.sf <- matrix(c(0,0.0005,0.0015) , nrow=no.param, ncol=3, byrow=T)
        param.sr <- matrix(c(0,0.0001,0.0002) , nrow=no.param, ncol=3, byrow=T)
        
        # param.Y <- matrix(c(0.2,0.3,0.4) , nrow=no.param, ncol=3, byrow=T)
        # param.af <- matrix(c(0,0.25,1) , nrow=no.param, ncol=3, byrow=T) # Forcing af not to be negetive
        # param.as <- matrix(c(0,0.65,1) , nrow=no.param, ncol=3, byrow=T)
        # param.sf <- matrix(c(0,0.0005,0.001) , nrow=no.param, ncol=3, byrow=T) # All Groups having same sf
        # param.sr <- matrix(c(0,0.0001,0.0002) , nrow=no.param, ncol=3, byrow=T) # All Groups having same sr
        
        if (with.storage==T) {
          param = data.frame(param.k,param.Y,param.af,param.as,param.sf,param.sr)
          names(param) <- c("k_min","k","k_max","Y_min","Y","Y_max","af_min","af","af_max","as_min","as","as_max","sf_min","sf","sf_max","sr_min","sr","sr_max")
          pMinima <- param[ ,c("k_min","Y_min","af_min","as_min","sf_min","sr_min")]
          pMaxima <- param[ ,c("k_max","Y_max","af_max","as_max","sf_max","sr_max")]
          pValues <- param[ ,c("k","Y","af","as","sf","sr")] # Starting point of the chain
        } else { # (with.storage==F) 
          param = data.frame(param.Y,param.af,param.as,param.sf,param.sr)
          names(param) <- c("Y_min","Y","Y_max","af_min","af","af_max","as_min","as","as_max","sf_min","sf","sf_max","sr_min","sr","sr_max")
          pMinima <- param[ ,c("Y_min","af_min","as_min","sf_min","sr_min")]
          pMaxima <- param[ ,c("Y_max","af_max","as_max","sf_max","sr_max")]
          pValues <- param[ ,c("Y","af","as","sf","sr")] # Starting point of the chain
        }
        pChain <- matrix(0, nrow=chainLength, ncol = no.param*no.var+1) # Initialising the chain
      }
      
      # Defining the variance-covariance matrix for proposal generation
      # vcov = (0.006*(pMaxima-pMinima))^2
      vcov = (0.05*(pMaxima-pMinima))^2
      # vcov = (0.0025*(pMaxima-pMinima))^2
      # vcov = (0.001*(pMaxima-pMinima))^2
      vcovProposal =  vcov # The higher the coefficient, the higher the deviations in parameter time series
      
      
      # Find the Prior probability density
      prior.dist = vector("list", no.var)
      for (i in 1:no.var) {
        prior.dist[i] = list(log(dnorm(pValues[ , i], (pMinima[ , i] + pMaxima[ , i])/2, (pMaxima[ , i] - pMinima[ , i])/3))) # Prior normal gaussian distribution
      }
      logPrior0 <- sum(unlist(prior.dist))
      
      # Calculating model outputs for the starting point of the chain
      for (j in 1:length(v)) {
        data.set = subset(data,(treatment.no %in% treat.group[v[j]]))
        if (no.param.par.var < 5) {
          if (with.storage==T) {
            # output.set = model(no.param,data.set,tnc.partitioning,Y=pValues$Y,k=pValues$k,af=pValues$af,as=pValues$as,sf=pValues$sf,sr=pValues$sr)
            output.set = model_v2(no.param,data.set,tnc.partitioning,Y=pValues$Y,k=pValues$k,af=pValues$af,as=pValues$as,sf=pValues$sf,sr=pValues$sr,v1,treat) # Model for different parameter variation
          } else {
            output.set = model.without.storage(no.param,data.set,Y=pValues$Y,af=pValues$af,as=pValues$as,sf=pValues$sf,sr=pValues$sr)
          }
        } else { # no.param.par.var > 5; monthly parameter setting)
          if (with.storage==T) {
            output.set = model.monthly(data.set,j,tnc.partitioning,Y=pValues$Y,k=pValues$k,af=pValues$af,as=pValues$as,sf=pValues$sf,sr=pValues$sr)
          } else {
            output.set = model.without.storage.monthly(data.set,j,Y=pValues$Y,af=pValues$af,as=pValues$as,sf=pValues$sf,sr=pValues$sr)
          }
        }
        # output.set$treatment = as.factor(treat[v1])
        output.set$treatment = as.factor(treat[v[j]])
        # output = output.set
        if (j == 1) {
          output = output.set
        }
        if (j > 1) {
          output = rbind(output,output.set)
        }
      }
      # # Consider the input uncertainty with u=10 model runs
      # listofdfs <- list()
      # # for (j in 1:length(v)) {
      #   data.set = subset(data,(treatment.no %in% treat.group[v[j]]))
      #   for (u in 1:10) {
      #     if (no.param.par.var < 5) {
      #       if (with.storage==T) {
      #         output.set = model(no.param,data.set,tnc.partitioning,Y=pValues$Y,k=pValues$k,af=pValues$af,as=pValues$as,sf=pValues$sf,sr=pValues$sr)
      #         output.set = model_v2(no.param,data.set,tnc.partitioning,Y=pValues$Y,k=pValues$k,af=pValues$af,as=pValues$as,sf=pValues$sf,sr=pValues$sr)
      #       } else {
      #         output.set = model.without.storage(no.param,data.set,Y=pValues$Y,af=pValues$af,as=pValues$as,sf=pValues$sf,sr=pValues$sr)
      #       }
      #     } else { # no.param.par.var > 5; monthly parameter setting)
      #       if (with.storage==T) {
      #         output.set = model.monthly(data.set,j,tnc.partitioning,Y=pValues$Y,k=pValues$k,af=pValues$af,as=pValues$as,sf=pValues$sf,sr=pValues$sr)
      #       } else {
      #         output.set = model.without.storage.monthly(data.set,j,Y=pValues$Y,af=pValues$af,as=pValues$as,sf=pValues$sf,sr=pValues$sr)
      #       }
      #     }
      #     listofdfs[[u]] <- output.set
      #   }
      #   output.set = aaply(laply(listofdfs, as.matrix), c(2, 3), mean)
      # 
      #   output.set$treatment = as.factor(treat[v1])
      #   output = output.set
      #   # if (j == 1) {
      #   #   output = output.set
      #   # }
      #   # if (j > 1) {
      #   #   output = rbind(output,output.set)
      #   # }
      # # }
      # output = as.data.frame(output)
      
      
      # output$treatment = as.factor(treat.group[v])
      # data.set = data.set[order(data.set$treatment),]
      logL0 <- logLikelihood.wtc3(no.param.par.var,data,output,with.storage,model.comparison) # Calculate log likelihood of starting point of the chain
      
      if (with.storage==T) {
        pChain[1,] <- c(pValues$k,pValues$Y,pValues$af,pValues$as,pValues$sf,pValues$sr,logL0) # Assign the first parameter set with log likelihood
      } else {
        pChain[1,] <- c(pValues$Y,pValues$af,pValues$as,pValues$sf,pValues$sr,logL0) # Assign the first parameter set with log likelihood
      }
      
      
      # Calculating the next candidate parameter vector, as a multivariate normal jump away from the current point
      for (c in (2 : chainLength)) {
        candidatepValues = matrix(ncol = no.var, nrow = no.param)
        for (i in 1:no.var) {
          candidatepValues[,i] = rmvnorm(n=1, mean=pValues[,i],
                                         sigma=diag(vcovProposal[,i],no.param)) 
        }
        candidatepValues = data.frame(candidatepValues)
        if (with.storage==T) {
          names(candidatepValues) <- c("k","Y","af","as","sf","sr")
        } else {
          names(candidatepValues) <- c("Y","af","as","sf","sr")
        }
        
        # Reflected back to generate another candidate value
        reflectionFromMin = pmin( unlist(matrix(0,nrow=no.param,ncol=no.var)), unlist(candidatepValues-pMinima) )
        reflectionFromMax = pmax( unlist(list(rep(0, no.param))), unlist(candidatepValues-pMaxima) )
        candidatepValues = candidatepValues - 2 * reflectionFromMin - 2 * reflectionFromMax 
        
        
        # Calculating the prior probability density for the candidate parameter vector
        if (all(candidatepValues>=pMinima) && all(candidatepValues<=pMaxima)){
          uni.dist = vector("list", no.var)
          for (i in 1:no.var) {
            # uni.dist[i] = list(log(dunif(candidatepValues[ , i], (pMinima[ , i] + pMaxima[ , i])/2, (pMaxima[ , i] - pMinima[ , i])/3))) # Prior normal gaussian distribution
            uni.dist[i] = list(log(dnorm(candidatepValues[ , i], (pMinima[ , i] + pMaxima[ , i])/2, (pMaxima[ , i] - pMinima[ , i])/3))) # Prior normal gaussian distribution
          }
          logPrior1 <- sum(unlist(uni.dist))
          Prior1 = 1
        } else {
          Prior1 <- 0
        }
        
        
        # Calculating the outputs for the candidate parameter vector and then log likelihood
        if (Prior1 > 0) {
          for (j in 1:length(v)) {
            data.set = subset(data,(treatment.no %in% treat.group[v[j]]))
            #   Mleaf = Mwood = Mroot = c()
            #   Mleaf[1] <- data.set$Mleaf[1]
            #   Mwood[1] <- data.set$Mwood[1]
            #   Mroot[1] <- data.set$Mroot[1]
            
            if (no.param.par.var < 5) {
              if (with.storage==T) {
                # out.cand.set = model(no.param,data.set,tnc.partitioning,candidatepValues$Y,
                #                      candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf,candidatepValues$sr)
                out.cand.set = model_v2(no.param,data.set,tnc.partitioning,candidatepValues$Y, 
                                        candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf,candidatepValues$sr,v1,treat)
              } else {
                out.cand.set = model.without.storage(no.param,data.set,candidatepValues$Y,
                                                     candidatepValues$af,candidatepValues$as,candidatepValues$sf,candidatepValues$sr)
              }
            } else { # no.param.par.var > 5; monthly parameter setting)
              if (with.storage==T) {
                out.cand.set = model.monthly(data.set,j,tnc.partitioning,candidatepValues$Y,
                                             candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf,candidatepValues$sr)
              } else {
                out.cand.set = model.without.storage.monthly(data.set,j,candidatepValues$Y,
                                                             candidatepValues$af,candidatepValues$as,candidatepValues$sf,candidatepValues$sr)
              }
            }
            out.cand.set$treatment = as.factor(treat[v[j]])
            # out.cand = out.cand.set
            if (j == 1) {
              out.cand = out.cand.set
            }
            if (j > 1) {
              out.cand = rbind(out.cand,out.cand.set)
            }
          }
          
          # # Consider the input uncertainty with u=10 model runs
          # listofdfs <- list()
          # # for (j in 1:length(v)) {
          #   data.set = subset(data,(treatment.no %in% treat.group[v[j]]))
          #   for (u in 1:10) {
          #     if (no.param.par.var < 5) {
          #       if (with.storage==T) {
          #         out.cand.set = model(no.param,data.set,tnc.partitioning,candidatepValues$Y,
          #                              candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf,candidatepValues$sr)
          #       } else {
          #         out.cand.set = model.without.storage(no.param,data.set,candidatepValues$Y,
          #                                              candidatepValues$af,candidatepValues$as,candidatepValues$sf,candidatepValues$sr)
          #       }
          #     } else { # no.param.par.var > 5; monthly parameter setting)
          #       if (with.storage==T) {
          #         out.cand.set = model.monthly(data.set,j,tnc.partitioning,candidatepValues$Y,
          #                                      candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf,candidatepValues$sr)
          #       } else {
          #         out.cand.set = model.without.storage.monthly(data.set,j,candidatepValues$Y,
          #                                                      candidatepValues$af,candidatepValues$as,candidatepValues$sf,candidatepValues$sr)
          #       }
          #     }
          #     listofdfs[[u]] <- out.cand.set
          #   }
          #   out.cand.set = aaply(laply(listofdfs, as.matrix), c(2, 3), mean)
          # 
          #   out.cand = as.data.frame(out.cand.set)
          #   out.cand.set$treatment = as.factor(treat[v1])
          #   out.cand = out.cand.set
          #   # if (j == 1) {
          #   #   out.cand = out.cand.set
          #   # }
          #   # if (j > 1) {
          #   #   out.cand = rbind(out.cand,out.cand.set)
          #   # }
          # # }
          
          # data = data[order(data$volume),]
          logL1 <- logLikelihood.wtc3(no.param.par.var,data,out.cand,with.storage,model.comparison) # Calculate log likelihood
          
          # Calculating the logarithm of the Metropolis ratio
          logalpha <- (logPrior1+logL1) - (logPrior0+logL0) 
          
          # Accepting or rejecting the candidate vector
          if (!treat[v1] %in% as.factor("amb+ele+const")) {
            if (no.param.par.var == 4) {
              if ( log(runif(1, min = 0, max =1)) < logalpha && (candidatepValues$sr[1] + candidatepValues$sr[2]*(nrow(data.set)) + candidatepValues$sr[3]*(nrow(data.set))^2 + candidatepValues$sr[4]*(nrow(data.set))^3) >= 0
                   && (candidatepValues$sf[1] + candidatepValues$sf[2]*(nrow(data.set)) + candidatepValues$sf[3]*(nrow(data.set))^2 + candidatepValues$sf[4]*(nrow(data.set))^3) >= 0
                   && (candidatepValues$k[1] + candidatepValues$k[2]*(nrow(data.set)) + candidatepValues$k[3]*(nrow(data.set))^2 + candidatepValues$k[4]*(nrow(data.set))^3) >= 0
                   && (candidatepValues$as[1] + candidatepValues$as[2]*(nrow(data.set)) + candidatepValues$as[3]*(nrow(data.set))^2 + candidatepValues$as[4]*(nrow(data.set))^3) >= 0
                   && (candidatepValues$af[1] + candidatepValues$af[2]*(nrow(data.set)) + candidatepValues$af[3]*(nrow(data.set))^2 + candidatepValues$af[4]*(nrow(data.set))^3) >= 0
                   # && (1-candidatepValues$af[1]+candidatepValues$as[1]) >= 0
                   && (1 - (candidatepValues$af[1] + candidatepValues$af[2]*(nrow(data.set)) + candidatepValues$af[3]*(nrow(data.set))^2 + candidatepValues$af[4]*(nrow(data.set))^3) - 
                       (candidatepValues$as[1] + candidatepValues$as[2]*(nrow(data.set)) + candidatepValues$as[3]*(nrow(data.set))^2 + candidatepValues$as[4]*(nrow(data.set))^3)) >= 0 ) {
                # && ((1-candidatepValues$af[1]-candidatepValues$as[1])*(nrow(data.set)) + (1-candidatepValues$af[2]-candidatepValues$as[2])*(nrow(data.set))^2 + (1-candidatepValues$af[3]-candidatepValues$as[3])*(nrow(data.set))^3 + (1-candidatepValues$af[4]-candidatepValues$as[4])*(nrow(data.set))^4) >= 0) {
                pValues <- candidatepValues
                logPrior0 <- logPrior1
                logL0 <- logL1
              }
            } else if (no.param.par.var == 3) {
              if ( log(runif(1, min = 0, max =1)) < logalpha && (candidatepValues$sr[1]*(nrow(data.set)) + candidatepValues$sr[2]*(nrow(data.set))^2 + candidatepValues$sr[3]*(nrow(data.set))^3) >= 0
                   && (candidatepValues$sf[1] + candidatepValues$sf[2]*(nrow(data.set)) + candidatepValues$sf[3]*(nrow(data.set))^2) >= 0
                   && (candidatepValues$k[1] + candidatepValues$k[2]*(nrow(data.set)) + candidatepValues$k[3]*(nrow(data.set))^2) >= 0
                   && (candidatepValues$as[1] + candidatepValues$as[2]*(nrow(data.set)) + candidatepValues$as[3]*(nrow(data.set))^2) >= 0
                   && (candidatepValues$af[1] + candidatepValues$af[2]*(nrow(data.set)) + candidatepValues$af[3]*(nrow(data.set))^2) >= 0
                   && (1 - (candidatepValues$af[1] + candidatepValues$af[2]*(nrow(data.set)) + candidatepValues$af[3]*(nrow(data.set))^2) - 
                       (candidatepValues$as[1] + candidatepValues$as[2]*(nrow(data.set)) + candidatepValues$as[3]*(nrow(data.set))^2)) >= 0 ) {
                pValues <- candidatepValues
                logPrior0 <- logPrior1
                logL0 <- logL1
              }
            } else if (no.param.par.var == 2) {
              if ( log(runif(1, min = 0, max =1)) < logalpha && (candidatepValues$sr[1]*(nrow(data.set)) + candidatepValues$sr[2]*(nrow(data.set))^2) >= 0
                   && (candidatepValues$sf[1] + candidatepValues$sf[2]*(nrow(data.set))) >= 0
                   # && (candidatepValues$k[1] + candidatepValues$k[2]*(nrow(data.set))) >= 0
                   && (candidatepValues$k[1]) >= 0
                   && (candidatepValues$as[1] + candidatepValues$as[2]*(nrow(data.set))) >= 0
                   && (candidatepValues$af[1] + candidatepValues$af[2]*(nrow(data.set))) >= 0
                   && (1 - (candidatepValues$af[1] + candidatepValues$af[2]*(nrow(data.set))) - 
                       (candidatepValues$as[1] + candidatepValues$as[2]*(nrow(data.set)))) >= 0 ) {
                pValues <- candidatepValues
                logPrior0 <- logPrior1
                logL0 <- logL1
              }
            } else if (no.param.par.var == 1) {
              if ( log(runif(1, min = 0, max =1)) < logalpha && (candidatepValues$sr[1]) >= 0
                   && (candidatepValues$sf[1]) >= 0
                   && (candidatepValues$k[1]) >= 0
                   && (candidatepValues$as[1]) >= 0
                   && (candidatepValues$af[1]) >= 0
                   && (1 - candidatepValues$af[1] - candidatepValues$as[1]) >= 0 ) {
                pValues <- candidatepValues
                logPrior0 <- logPrior1
                logL0 <- logL1
              }
            } else {
              if ( log(runif(1, min = 0, max =1)) < logalpha && candidatepValues >= 0 && all((candidatepValues$af + candidatepValues$as) < 1)) {
                pValues <- candidatepValues
                logPrior0 <- logPrior1
                logL0 <- logL1
              }
            }
          } else if (treat[v1] %in% as.factor("amb+ele+const")) {
            if (no.param.par.var == 4) {
              if ( log(runif(1, min = 0, max =1)) < logalpha && (candidatepValues$sr[1] + candidatepValues$sr[2]*(nrow(data.set)) + candidatepValues$sr[3]*(nrow(data.set))^2 + candidatepValues$sr[4]*(nrow(data.set))^3) >= 0
                   && (candidatepValues$sf[1] + candidatepValues$sf[2]*(nrow(data.set)) + candidatepValues$sf[3]*(nrow(data.set))^2 + candidatepValues$sf[4]*(nrow(data.set))^3) >= 0
                   && (candidatepValues$k[1] + candidatepValues$k[2]*(nrow(data.set)) + candidatepValues$k[3]*(nrow(data.set))^2 + candidatepValues$k[4]*(nrow(data.set))^3) >= 0
                   && (candidatepValues$as[1] + candidatepValues$as[2]*(nrow(data.set)) + candidatepValues$as[3]*(nrow(data.set))^2 + candidatepValues$as[4]*(nrow(data.set))^3) >= 0
                   && (candidatepValues$af[1] + candidatepValues$af[2]*(nrow(data.set)) + candidatepValues$af[3]*(nrow(data.set))^2 + candidatepValues$af[4]*(nrow(data.set))^3) >= 0
                   # && (1-candidatepValues$af[1]+candidatepValues$as[1]) >= 0
                   && (1 - (candidatepValues$af[1] + candidatepValues$af[2]*(nrow(data.set)) + candidatepValues$af[3]*(nrow(data.set))^2 + candidatepValues$af[4]*(nrow(data.set))^3) - 
                       (candidatepValues$as[1] + candidatepValues$as[2]*(nrow(data.set)) + candidatepValues$as[3]*(nrow(data.set))^2 + candidatepValues$as[4]*(nrow(data.set))^3)) >= 0 ) {
                # && ((1-candidatepValues$af[1]-candidatepValues$as[1])*(nrow(data.set)) + (1-candidatepValues$af[2]-candidatepValues$as[2])*(nrow(data.set))^2 + (1-candidatepValues$af[3]-candidatepValues$as[3])*(nrow(data.set))^3 + (1-candidatepValues$af[4]-candidatepValues$as[4])*(nrow(data.set))^4) >= 0) {
                pValues <- candidatepValues
                logPrior0 <- logPrior1
                logL0 <- logL1
              }
            } else if (no.param.par.var == 3) {
              if ( log(runif(1, min = 0, max =1)) < logalpha && (candidatepValues$sr[1]*(nrow(data.set)) + candidatepValues$sr[2]*(nrow(data.set))^2 + candidatepValues$sr[3]*(nrow(data.set))^3) >= 0
                   && (candidatepValues$sf[1] + candidatepValues$sf[2]*(nrow(data.set)) + candidatepValues$sf[3]*(nrow(data.set))^2) >= 0
                   && (candidatepValues$k[1] + candidatepValues$k[2]*(nrow(data.set)) + candidatepValues$k[3]*(nrow(data.set))^2) >= 0
                   && (candidatepValues$as[1] + candidatepValues$as[2]*(nrow(data.set)) + candidatepValues$as[3]*(nrow(data.set))^2) >= 0
                   && (candidatepValues$af[1] + candidatepValues$af[2]*(nrow(data.set)) + candidatepValues$af[3]*(nrow(data.set))^2) >= 0
                   && (1 - (candidatepValues$af[1] + candidatepValues$af[2]*(nrow(data.set)) + candidatepValues$af[3]*(nrow(data.set))^2) - 
                       (candidatepValues$as[1] + candidatepValues$as[2]*(nrow(data.set)) + candidatepValues$as[3]*(nrow(data.set))^2)) >= 0 ) {
                pValues <- candidatepValues
                logPrior0 <- logPrior1
                logL0 <- logL1
              }
            } else if (no.param.par.var == 2) {
              if ( log(runif(1, min = 0, max =1)) < logalpha && (candidatepValues$sr[1]*(nrow(data.set)) + candidatepValues$sr[2]*(nrow(data.set))^2) >= 0
                   && (candidatepValues$sf[1] + candidatepValues$sf[2]*(nrow(data.set))) >= 0
                   # && (candidatepValues$k[1] + candidatepValues$k[2]*(nrow(data.set))) >= 0
                   && (candidatepValues$k[1]) >= 0
                   && (candidatepValues$as[1] + candidatepValues$as[2]*(nrow(data.set))) >= 0
                   && (candidatepValues$af[1] + candidatepValues$af[2]*(nrow(data.set))) >= 0
                   && (1 - (candidatepValues$af[1] + candidatepValues$af[2]*(nrow(data.set))) - 
                       (candidatepValues$as[1] + candidatepValues$as[2]*(nrow(data.set)))) >= 0 ) {
                pValues <- candidatepValues
                logPrior0 <- logPrior1
                logL0 <- logL1
              }
            } else if (no.param.par.var == 1) {
              if ( log(runif(1, min = 0, max =1)) < logalpha && (candidatepValues$sr[1]) >= 0
                   && (candidatepValues$sf[1]) >= 0
                   && (candidatepValues$k[1]) >= 0
                   && (candidatepValues$as[1]) >= 0
                   && (candidatepValues$af[1]) >= 0
                   && (1 - candidatepValues$af[1] - candidatepValues$as[1]) >= 0 ) {
                pValues <- candidatepValues
                logPrior0 <- logPrior1
                logL0 <- logL1
              }
            } else {
              if ( log(runif(1, min = 0, max =1)) < logalpha && candidatepValues >= 0 && all((candidatepValues$af + candidatepValues$as) < 1)) {
                pValues <- candidatepValues
                logPrior0 <- logPrior1
                logL0 <- logL1
              }
            }
          }
        }
        if (with.storage==T) {
          pChain[c,] <- c(pValues$k,pValues$Y,pValues$af,pValues$as,pValues$sf,pValues$sr,logL0) # Assign the first parameter set with log likelihood
        } else {
          pChain[c,] <- c(pValues$Y,pValues$af,pValues$as,pValues$sf,pValues$sr,logL0) # Assign the first parameter set with log likelihood
        }
      }
      
      # Discard the first 500 iterations for Burn-IN in MCMC
      pChain <- pChain[(bunr_in+1):nrow(pChain),]
      pChain = as.data.frame(pChain)
      if (no.param.par.var < 5) {
        if (with.storage==T) {
          if (no.param.par.var[z]==1) {
            names(pChain) <- c("k1","Y1","af1","as1","sf1","sr1","logli")
          } else if (no.param.par.var[z]==2) {
            names(pChain) <- c("k1","k2","Y1","Y2","af1","af2","as1","as2","sf1","sf2","sr1","sr2","logli")
          } else if (no.param.par.var[z]==3) {
            names(pChain) <- c("k1","k2","k3","Y1","Y2","Y3","af1","af2","af3","as1","as2","as3","sf1","sf2","sf3","sr1","sr2","sr3","logli")
          }
        } else {
          if (no.param.par.var[z]==1) {
            names(pChain) <- c("Y1","af1","as1","sf1","sr1","logli")
          } else if (no.param.par.var[z]==2) {
            names(pChain) <- c("Y1","Y2","af1","af2","as1","as2","sf1","sf2","sr1","sr2","logli")
          } else if (no.param.par.var[z]==3) {
            names(pChain) <- c("Y1","Y2","Y3","af1","af2","af3","as1","as2","as3","sf1","sf2","sf3","sr1","sr2","sr3","logli")
          }
        }
      }
      
      # Store the final parameter set values
      param.set = colMeans(pChain[ , 1:(no.param*no.var)])
      param.SD = apply(pChain[ , 1:(no.param*no.var)], 2, sd)
      param.final = data.frame(matrix(ncol = (no.var)*2, nrow = no.param))
      if (with.storage==T) {
        names(param.final) <- c("k","Y","af","as","sf","sr","k_SD","Y_SD","af_SD","as_SD","sf_SD","sr_SD")
        # param.final$k = param.set[1:no.param]
        # param.final$Y = param.set[(1+no.param):(2*no.param)]
        param.final$k = mean(param.set[1:no.param])
        param.final$Y = mean(param.set[(1+no.param):(2*no.param)])
        param.final$af = param.set[(1+2*no.param):(3*no.param)]
        param.final$as = param.set[(1+3*no.param):(4*no.param)]
        param.final$sf = param.set[(1+4*no.param):(5*no.param)]
        param.final$sr = param.set[(1+5*no.param):(6*no.param)]
        
        # param.final$k_SD = param.SD[1:no.param]
        # param.final$Y_SD = param.SD[(1+no.param):(2*no.param)]
        param.final$k_SD = mean(param.SD[1:no.param])
        param.final$Y_SD = mean(param.SD[(1+no.param):(2*no.param)])
        param.final$af_SD = param.SD[(1+2*no.param):(3*no.param)]
        param.final$as_SD = param.SD[(1+3*no.param):(4*no.param)]
        param.final$sf_SD = param.SD[(1+4*no.param):(5*no.param)]
        param.final$sr_SD = param.SD[(1+5*no.param):(6*no.param)]
      } else {
        names(param.final) <- c("Y","af","as","sf","sr","Y_SD","af_SD","as_SD","sf_SD","sr_SD")
        param.final$Y = param.set[1:no.param]
        param.final$af = param.set[(1+no.param):(2*no.param)]
        param.final$as = param.set[(1+2*no.param):(3*no.param)]
        param.final$sf = param.set[(1+3*no.param):(4*no.param)]
        param.final$sr = param.set[(1+4*no.param):(5*no.param)]
        
        param.final$Y_SD = param.SD[1:no.param]
        param.final$af_SD = param.SD[(1+no.param):(2*no.param)]
        param.final$as_SD = param.SD[(1+2*no.param):(3*no.param)]
        param.final$sf_SD = param.SD[(1+3*no.param):(4*no.param)]
        param.final$sr_SD = param.SD[(1+4*no.param):(5*no.param)]
      }
      
      # Calculate final output set from the predicted parameter set
      for (j in 1:length(v)) {
        data.set = subset(data,(treatment.no %in% treat.group[v[j]]))
        #   Mleaf = Mwood = Mroot = c()
        #   Mleaf[1] <- data.set$Mleaf[1]
        #   Mwood[1] <- data.set$Mwood[1]
        #   Mroot[1] <- data.set$Mroot[1]
        
        if (no.param.par.var < 5) {
          if (with.storage==T) {
            # output.final.set = model(no.param,data.set,tnc.partitioning,param.final$Y,
            #                          param.final$k,param.final$af,param.final$as,param.final$sf,param.final$sr)
            output.final.set = model_v2(no.param,data.set,tnc.partitioning,param.final$Y,
                                        param.final$k,param.final$af,param.final$as,param.final$sf,param.final$sr)
          } else {
            output.final.set = model.without.storage(no.param,data.set,param.final$Y,
                                                     param.final$af,param.final$as,param.final$sf,param.final$sr)
          }
        } else { # no.param.par.var > 5; monthly parameter setting) 
          if (with.storage==T) {
            output.final.set = model.monthly(data.set,j,tnc.partitioning,param.final$Y,
                                             param.final$k,param.final$af,param.final$as,param.final$sf,param.final$sr)
          } else {
            output.final.set = model.without.storage.monthly(data.set,j,param.final$Y,
                                                             param.final$af,param.final$as,param.final$sf,param.final$sr)
          }
        }
        
        output.final.set$treatment = as.factor(treat[v[j]])
        # output.final = output.final.set
        if (j == 1) {
          output.final = output.final.set
        }
        if (j > 1) {
          output.final = rbind(output.final,output.final.set)
        }
      }
      
      # #----------------------------------------------------------------------------------------------------------------
      # if (with.storage==T) {
      #   output.final$Sleaf = output.final$Sleaf / output.final$Mleaf * 100
      # }
      # #----------------------------------------------------------------------------------------------------------------
      
      if (no.param.par.var < 5) {
        # Calculate daily parameter values with SD
        Days <- seq(1,nrow(data.set), length.out=nrow(data.set))
        param.daily = param.final[1,]
        
        ###################
        for (i in 2:length(Days)) {
          param.daily[i,c("k","Y","k_SD","Y_SD")] = param.final[1,c("k","Y","k_SD","Y_SD")]
        }
        if (no.param == 2) {
          for (i in 2:length(Days)) {
            param.daily[i,3:no.var] = param.final[1,3:no.var] + param.final[2,3:no.var] * i
          }
          for (i in (no.var+3):(2*no.var)) {
            param.daily[,i] = ((param.final[1,i]^2 + param.final[2,i]^2)/2)^0.5
          }
        }
        ###################
        
        if (no.param == 1) {
          for (i in 2:length(Days)) {
            param.daily[i,] = param.final[1,]
          }
        }
        # if (no.param == 2) {
        #   for (i in 2:length(Days)) {
        #     param.daily[i,1:no.var] = param.final[1,1:no.var] + param.final[2,1:no.var] * i
        #   }
        #   for (i in (no.var+1):(2*no.var)) {
        #     param.daily[,i] = ((param.final[1,i]^2 + param.final[2,i]^2)/2)^0.5
        #   }
        # }
        if (no.param == 3) {
          for (i in 2:length(Days)) {
            param.daily[i,1:no.var] = param.final[1,1:no.var] + param.final[2,1:no.var] * i + param.final[3,1:no.var] * i^2
          }
          for (i in (no.var+1):(2*no.var)) {
            param.daily[,i] = ((param.final[1,i]^2 + param.final[2,i]^2 + param.final[3,i]^2)/3)^0.5
          }
        }
        if (no.param == 4) {
          for (i in 2:length(Days)) {
            param.daily[i,1:no.var] = param.final[1,1:no.var] + param.final[2,1:no.var] * i + param.final[3,1:no.var] * i^2 + param.final[4,1:no.var] * i^3
          }
          for (i in (no.var+1):(2*no.var)) {
            param.daily[,i] = ((param.final[1,i]^2 + param.final[2,i]^2 + param.final[3,i]^2 + param.final[4,i]^3)/4)^0.5
          }
        }
        param.daily$ar = 1 - param.daily$af - param.daily$as
        param.daily$ar_SD = with(param.daily, ((af_SD*af_SD + as_SD*as_SD)/2)^0.5)
        param.daily$Date = as.Date(data.set$Date)
        
        
        # Plotting the parameter sets over time
        if (with.storage==T) { 
          melted.param1 = melt(param.daily[,c("k","Y","af","as","ar","sf","sr","Date")], id.vars="Date")
          melted.param2 = melt(param.daily[,c("k_SD","Y_SD","af_SD","as_SD","ar_SD","sf_SD","sr_SD","Date")], id.vars="Date")
        } else {
          melted.param1 = melt(param.daily[,c("Y","af","as","ar","sf","sr","Date")], id.vars="Date")
          melted.param2 = melt(param.daily[,c("Y_SD","af_SD","as_SD","ar_SD","sf_SD","sr_SD","Date")], id.vars="Date")
        }
        melted.param = data.frame(melted.param1$Date, melted.param1$variable, melted.param1$value, melted.param2$value)
        names(melted.param) = c("Date","variable","Parameter","Parameter_SD")
        melted.param$Date = as.Date(melted.param$Date)
        melted.param$treatment = as.factor(treat[v1])
        # melted.param$treat.type = as.factor(v1)
        melted.param$no.param = as.factor(no.param.par.var[z])
      } else { # no.param.par.var > 5; monthly parameter setting) 
        
        # Parameter sets over time
        param.final$ar = 1 - param.final$af - param.final$as
        param.final$ar_SD = with(param.final, (af_SD*af_SD + as_SD*as_SD)^0.5)
        
        param.final$Date = data.set$Date[seq(1,nrow(data.set),param.vary)]
        melted.param1 = melt(param.final[,c("k","Y","af","as","ar","sf","sr","Date")], id.vars="Date")
        melted.param2 = melt(param.final[,c("k_SD","Y_SD","af_SD","as_SD","ar_SD","sf_SD","sr_SD","Date")], id.vars="Date")
        melted.param = data.frame(melted.param1$Date, melted.param1$variable, melted.param1$value, melted.param2$value)
        names(melted.param) = c("Date","variable","Parameter","Parameter_SD")
        melted.param$Date = as.Date(melted.param$Date)
        melted.param$treatment = as.factor(treat[v1])
        melted.param$treat.type = as.factor(length(v))
        melted.param$no.param = as.factor(ceiling(no.param.par.var[z]))
      }
      
      # Plotting the Measured (data) vs Modelled Plant Carbon pools for plotting and comparison
      #----------------------------------------------------------------------------------------------------------------
      # lm.daily = read.csv("processed_data/Cleaf_daily_data.csv") # Unit gC per gC plant
      #----------------------------------------------------------------------------------------------------------------
      
      for (j in 1:length(v)) {
        data.set = subset(data,(treatment.no %in% treat.group[v[j]]))
        
        output.final.set = subset(output.final,(treatment %in% treat[v[j]]))
        output.final.set$Date = data.set$Date
        
        if (with.storage==T) {
          names(output.final.set) = c("Cstorage.modelled","Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Sleaf.modelled","Rm","Rabove","Treatment","Date")
          melted.output = melt(output.final.set[,c("Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Sleaf.modelled","Rm","Rabove","Date")], id.vars="Date")
          melted.Cstorage = output.final.set[,c("Cstorage.modelled","Treatment","Date")]
          melted.data = melt(data.set[ , c("LM","WM","RM","litter","TNC_leaf","Ra","Date")], id.vars="Date")
          melted.error = melt(data.set[ , c("LM_SE","WM_SE","RM_SE","litter_SE","TNC_leaf_SE","Ra_SE","Date")], id.vars="Date")
          melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
          melted.Cstorage$treatment = as.factor(treat[v[j]])
          melted.Cstorage$treat.type = as.factor(length(v))
          melted.Cstorage$no.param = as.factor(ceiling(no.param.par.var[z]))
        } else {
          names(output.final.set) = c("Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Rm","Rabove","Treatment","Date")
          melted.output = melt(output.final.set[,c("Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Rm","Rabove","Date")], id.vars="Date")
          melted.data = melt(data.set[ , c("LM","WM","RM","litter","Ra","Date")], id.vars="Date")
          melted.error = melt(data.set[ , c("LM_SE","WM_SE","RM_SE","litter_SE","Ra_SE","Date")], id.vars="Date")
        }
        melted.output$Date = as.Date(melted.output$Date)
        melted.output$treatment = as.factor(treat[v[j]])
        melted.output$treat.type = as.factor(length(v))
        melted.output$no.param = as.factor(ceiling(no.param.par.var[z]))
        
        melted.data$Date = as.Date(melted.data$Date)
        melted.data$treatment = as.factor(treat[v[j]])
        melted.data$treat.type = as.factor(length(v))
        
        melted.error$Date = as.Date(melted.error$Date)
        melted.error$treatment = as.factor(treat[v[j]])
        melted.error$treat.type = as.factor(length(v))
        melted.error$parameter = melted.data$value
        melted.error$no.param = as.factor(ceiling(no.param.par.var[z]))
        
        # if (v1 < 8){
        #   melted.output$volume.group = as.factor(1)
        #   if (with.storage==T) {
        #     melted.Cstorage$volume.group = as.factor(1)
        #   }
        #   melted.error$volume.group = as.factor(1)
        # }
        # if (v1 == 8){
        #   melted.output$volume.group = as.factor(2)
        #   if (with.storage==T) {
        #     melted.Cstorage$volume.group = as.factor(2)
        #   }
        #   melted.error$volume.group = as.factor(2)
        # }
        
        
        # Storing the summary of this volume group of data, outputs, Cstorage (Parameter is same for the group, will be stored later)
        if (j == 1) {
          summary.data.set = melted.data
          summary.error.set = melted.error
          summary.output.set = melted.output
          if (with.storage==T) {
            summary.Cstorage.set = melted.Cstorage
          }
        }
        if (j > 1) {
          summary.output.set = rbind(summary.output.set,melted.output)
          if (with.storage==T) {
            summary.Cstorage.set = rbind(summary.Cstorage.set,melted.Cstorage)
          }
          summary.error.set = rbind(summary.error.set,melted.error)
          if (z == 1) {
            summary.data.set = rbind(summary.data.set,melted.data)
          }
        }
      }
      
      # # Measured (data) vs Modelled Plant Carbon pools for plotting and comparison
      # output.final$Date = data.set$Date
      # if (with.storage==T) { 
      #   names(output.final) = c("Cstorage.modelled","Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Sleaf.modelled","Rm","Rabove","Date")
      #   melted.output = melt(output.final[,c("Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Sleaf.modelled","Rm","Rabove","Date")], id.vars="Date")
      #   melted.Cstorage = output.final[,c("Cstorage.modelled","Date")]
      #   melted.data = melt(data.set[ , c("LM","WM","RM","litter","TNC_leaf","Ra","Date")], id.vars="Date")
      #   melted.error = melt(data.set[ , c("LM_SE","WM_SE","RM_SE","litter_SE","TNC_leaf_SE","Ra_SE","Date")], id.vars="Date")
      #   melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
      #   melted.Cstorage$treatment = as.factor(treat.group[v])
      #   melted.Cstorage$no.param = as.factor(ceiling(no.param.par.var[z]))
      # } else {
      #   names(output.final) = c("Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Rm","Rabove","Date")
      #   melted.output = melt(output.final[,c("Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Rm","Rabove","Date")], id.vars="Date")
      #   melted.data = melt(data.set[ , c("LM","WM","RM","litter","Ra","Date")], id.vars="Date")
      #   melted.error = melt(data.set[ , c("LM_SE","WM_SE","RM_SE","litter_SE","Ra_SE","Date")], id.vars="Date")
      # }
      # melted.data$Date = as.Date(melted.data$Date)
      # melted.error$Date = as.Date(melted.error$Date)
      # melted.error$treatment = as.factor(treat.group[v])
      # melted.error$parameter = melted.data$value
      # melted.output$Date = as.Date(melted.output$Date)
      # melted.data$treatment = as.factor(treat.group[v])
      # melted.output$treatment = as.factor(treat.group[v])
      # melted.output$no.param = as.factor(ceiling(no.param.par.var[z]))
      # 
      # # Storing the summary of data, outputs, Cstorage, parameters
      # if (q == 1) {
      #   summary.data = melted.data
      #   summary.error = melted.error
      #   summary.output = melted.output
      #   if (with.storage==T) {
      #     summary.Cstorage = melted.Cstorage
      #   }
      #   summary.param = melted.param
      # }
      # if (q > 1) {
      #   summary.output = rbind(summary.output,melted.output)
      #   if (with.storage==T) {
      #     summary.Cstorage = rbind(summary.Cstorage,melted.Cstorage)
      #   }
      #   summary.param = rbind(summary.param,melted.param)
      #   summary.error = rbind(summary.error,melted.error)
      #   if (z == 1) {
      #     summary.data = rbind(summary.data,melted.data)
      #   }
      # }
      
      # Storing the summary of all volume group's data, outputs, Cstorage, parameters
      if (q == 1) {
        summary.data = summary.data.set
        summary.error = summary.error.set
        summary.output = summary.output.set
        if (with.storage==T) {
          summary.Cstorage = summary.Cstorage.set
        }
        summary.param = melted.param
      }
      if (q > 1) {
        summary.output = rbind(summary.output,summary.output.set)
        if (with.storage==T) {
          summary.Cstorage = rbind(summary.Cstorage,summary.Cstorage.set)
        }
        summary.param = rbind(summary.param,melted.param)
        # if (v1 < length(treat.group)) {
        summary.error = rbind(summary.error,summary.error.set)
        if (z == 1) {
          summary.data = rbind(summary.data,summary.data.set)
        }
        # }
      }
      
      
      # Display the Acceptance rate of the chain
      nAccepted = length(unique(pChain[,1]))
      acceptance = (paste("Treatment =",treat[v1],", Total Parameter number =",ceiling(no.param.par.var[z]),": ", nAccepted, "out of ", chainLength-bunr_in, "candidates accepted ( = ",
                          round(100*nAccepted/chainLength), "%)"))
      print(acceptance)
      
      
      # Plotting all parameter whole iterations for Day 1 only to check the convergance
      png(file = paste("output/Parameter_iterations_day1_treatment_",treat.group[v],"_par_",no.param.par.var[z], ".png", sep = ""))
      par(mfrow=c(3,3),oma = c(0, 0, 2, 0))
      if (with.storage==T) { 
        plot(pChain[,1],col="red",main="Utilization coefficient at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="k",ylim=c(param.k[1,1],param.k[1,3]))
        plot(pChain[,1+no.param],col="green",main="Alloc frac to Biomass at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="Y",ylim=c(param.Y[1,1],param.Y[1,3]))
        plot(pChain[,1+2*no.param],col="magenta",main="Alloc frac to foliage at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="af",ylim=c(param.af[1,1],param.af[1,3]))
        plot(pChain[,1+3*no.param],col="blue",main="Alloc frac to wood at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="as",ylim=c(param.as[1,1],param.as[1,3]))
        plot(pChain[,1+4*no.param],col="green",main="Foliage turnover at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="sf",ylim=c(param.sf[1,1],param.sf[1,3]))
        plot(pChain[,1+5*no.param],col="magenta",main="Root turnover at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="sr",ylim=c(param.sr[1,1],param.sr[1,3]))
        plot(pChain[,1+6*no.param],col="cyan",main="Log-likelihood",cex.lab = 1.5,xlab="Iterations",ylab="Log-likelihood")
      } else {
        plot(pChain[,1],col="green",main="Alloc frac to Biomass at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="Y",ylim=c(param.Y[1,1],param.Y[1,3]))
        plot(pChain[,1+no.param],col="magenta",main="Alloc frac to foliage at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="af",ylim=c(param.af[1,1],param.af[1,3]))
        plot(pChain[,1+2*no.param],col="blue",main="Alloc frac to wood at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="as",ylim=c(param.as[1,1],param.as[1,3]))
        plot(pChain[,1+3*no.param],col="green",main="Foliage turnover at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="sf",ylim=c(param.sf[1,1],param.sf[1,3]))
        plot(pChain[,1+4*no.param],col="magenta",main="Root turnover at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="sr",ylim=c(param.sr[1,1],param.sr[1,3]))
        plot(pChain[,1+5*no.param],col="cyan",main="Log-likelihood",cex.lab = 1.5,xlab="Iterations",ylab="Log-likelihood")
      }
      title(main = paste("First day Parameter iterations for treatment group",treat.group[v],"with par",no.param.par.var[z]), outer=TRUE, cex = 1.5)
      dev.off()
      
      # Calculate LogLi, AIC, BIC, Time to find the most accurate model for best balance between model fit and complexity
      output.final1 = output.final
      if (with.storage==T) { 
        names(output.final1) = c("Cstorage","Mleaf","Mwood","Mroot","Mlit","Sleaf","Rm","Rabove","Treatment") # Rename for the logLikelihood function
      } else {
        names(output.final1) = c("Mleaf","Mwood","Mroot","Mlit","Rm","Rabove","Treatment")
      }
      
      # data = data[with(data, order(volume)), ]
      # row.names(data) = c(1:nrow(data))
      # aic.bic[q,1] <- logLikelihood.wtc3(no.param.par.var,data.set,output.final1,with.storage,model.comparison) # Calculate logLikelihood
      aic.bic[q,1] <- logLikelihood.wtc3.final(no.param.par.var,data,output.final1,with.storage,model.comparison) # Calculate logLikelihood
      
      k1 = 2 # k = 2 for the usual AIC
      npar = no.param*no.var # npar = total number of parameters in the fitted model
      aic.bic[q,2] = -2*aic.bic[q,1] + k1*npar
      
      if (model.comparison==F) {
        n = sum(!is.na(data.set$TNC_leaf)) + sum(!is.na(data.set$LM)) + sum(!is.na(data.set$WM)) + sum(!is.na(data.set$RM)) + sum(!is.na(data.set$litter))
      } else {
        n = sum(!is.na(data.set$LM)) + sum(!is.na(data.set$WM)) + sum(!is.na(data.set$RM)) + sum(!is.na(data.set$litter))
      }
      k2 = log(n) # n being the number of observations for the so-called BIC
      aic.bic[q,3] = -2*aic.bic[q,1] + k2*npar
      
      time$end.time[q] <- Sys.time()
      time$time.taken[q] <- time$end.time[q] - time$start.time[q]
      aic.bic[q,4] = time$time.taken[q]
      aic.bic[q,5] = no.param
      aic.bic[q,6] = as.character(treat[v1])
    }
  }
  aic.bic$treatment = as.factor(aic.bic$treatment)
  bic = data.frame(aic.bic[,c("bic","no.param","treatment")])
  # melted.aic.bic = melt(aic.bic[,c(1:5)], id.vars=c("no.param"))
  
  # if (model.comparison==T | model.optimization==T) {
  if (model.optimization==T) {
    return(bic)
  } else if (model.optimization==F & with.storage==T) {
    result = list(no.param,summary.param,summary.data,summary.output,summary.error,bic,summary.Cstorage)
    return(result)
  } else {
    result = list(no.param,summary.param,summary.data,summary.output,summary.error,bic)
    return(result)
  }
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#-- Function to run Carbon Balance Model (CBM). 
#-----------------------------------------------------------------------------------------
# This script define the model equations to carry out Bayesian calibration for 
# 5 variables (allocation fractions: "k","Y",af","as","sf") on 
# various temporal scales to estimate Carbon pools (Cstorage,Cleaf,Cwood,Croot)
#-----------------------------------------------------------------------------------------
# Defining the model to iteratively calculate Cstorage, Cleaf, Cwood, Croot, Sleaf, Swood, Sroot
model <- function (no.param,data.set,tnc.partitioning,Y,k,af,as,sf,sr) {
  Mleaf = Mwood = Mroot = Mlit = Y.modelled = Rabove = c()
  Mleaf[1] <- data.set$LM[1]
  Mwood[1] <- data.set$WM[1]
  Mroot[1] <- data.set$RM[1]
  Mlit[1] <- data.set$litter[1]
  
  Cstorage = Sleaf = Swood = Sroot = c()
  
  # From WTC-4 experiment for TNC partitioning to tree organs
  Sleaf[1] = data.set$TNC_leaf[min(which(complete.cases(data.set$TNC_leaf)))]*0.75 # Consider the first available leaf tnc data as the starting point
  Swood[1] = data.set$TNC_wood[min(which(complete.cases(data.set$TNC_wood)))]*0.75 # Consider the first available leaf tnc data as the starting point
  Sroot[1] = data.set$TNC_root[min(which(complete.cases(data.set$TNC_root)))]*0.75 # Consider the first available leaf tnc data as the starting point
  Cstorage[1] = data.set$TNC_tot[min(which(complete.cases(data.set$TNC_tot)))]*0.75 # Consider the first available leaf tnc data as the starting point
  
  # Sleaf[1] = 10
  # Swood[1] = 5
  # Sroot[1] = 5
  # Cstorage[1] = 20
  
  Cleaf <- Croot <- Cwood <- Rm <- c()
  Cleaf[1] <- data.set$LM[1] - Sleaf[1]
  Cwood[1] <- data.set$WM[1] - Swood[1]
  Croot[1] <- data.set$RM[1] - Sroot[1]
  
  Rabove[1] = data.set$Ra[1]
  Rm[1] = data.set$Rd.foliage.mean[1]*Mleaf[1] + data.set$Rd.stem.mean[1]*Mwood[1]*data.set$SMratio[1] +
    data.set$Rd.branch.mean[1]*Mwood[1]*data.set$BMratio[1] +
    data.set$Rd.fineroot.mean[1]*Mroot[1]*data.set$FRratio[1] + data.set$Rd.intermediateroot.mean[1]*Mroot[1]*data.set$IRratio[1] +
    data.set$Rd.coarseroot.mean[1]*Mroot[1]*data.set$CRratio[1] + data.set$Rd.boleroot.mean[1]*Mroot[1]*data.set$BRratio[1]
  
  # Y.modelled = Y[1];
  if (no.param == 1) {
    k.i = k[1]; Y.i = Y[1]; af.i = af[1]; as.i = as[1]; sf.i = sf[1]; sr.i = sr[1]
  }
  for (i in 2:nrow(data.set)) {
    if (no.param == 2) {
      k.i = k[1] + k[2]*i; Y.i = Y[1]+ Y[2]*i; af.i = af[1]+ af[2]*i; as.i = as[1]+ as[2]*i; sf.i = sf[1]+ sf[2]*i; sr.i = sr[1]+ sr[2]*i
    }
    if (no.param == 3) {
      k.i = k[1] + k[2]*i + k[3]*i*i; Y.i = Y[1]+ Y[2]*i + Y[3]*i*i; af.i = af[1]+ af[2]*i + af[3]*i*i;
      as.i = as[1]+ as[2]*i + as[3]*i*i; sf.i = sf[1]+ sf[2]*i + sf[3]*i*i; sr.i = sr[1]+ sr[2]*i + sr[3]*i*i
    }
    if (no.param == 4) {
      k.i = k[1] + k[2]*i + k[3]*i*i + k[4]*i*i*i; Y.i = Y[1]+ Y[2]*i + Y[3]*i*i + Y[4]*i*i*i; af.i = af[1]+ af[2]*i + af[3]*i*i + af[4]*i*i*i;
      as.i = as[1]+ as[2]*i + as[3]*i*i + as[4]*i*i*i; sf.i = sf[1]+ sf[2]*i + sf[3]*i*i + sf[4]*i*i*i; sr.i = sr[1]+ sr[2]*i + sr[3]*i*i + sr[4]*i*i*i
    }
    # if (no.param == 1) {
    #   k.i = k[1]*data.all$Tair[i]; Y.i = Y[1]*data.all$Tair[i]; af.i = af[1]*data.all$Tair[i]; 
    #   as.i = as[1]*data.all$Tair[i]; sf.i = sf[1]*data.all$Tair[i]; sr.i = sr[1]*data.all$Tair[i]
    # }
    # if (no.param == 2) {
    #   k.i = k[1] + k[2]*data.all$Tair[i]; Y.i = Y[1]+ Y[2]*data.all$Tair[i]; af.i = af[1]+ af[2]*data.all$Tair[i]; 
    #   as.i = as[1]+ as[2]*data.all$Tair[i]; sf.i = sf[1]+ sf[2]*data.all$Tair[i]; sr.i = sr[1]+ sr[2]*data.all$Tair[i]
    # }
    # if (no.param == 3) {
    #   k.i = k[1] + k[2]*data.all$Tair[i] + k[3]*data.all$Tair[i]*data.all$Tair[i]; Y.i = Y[1]+ Y[2]*data.all$Tair[i] + Y[3]*data.all$Tair[i]*data.all$Tair[i]; 
    #   af.i = af[1]+ af[2]*data.all$Tair[i] + af[3]*data.all$Tair[i]*data.all$Tair[i]; as.i = as[1]+ as[2]*data.all$Tair[i] + as[3]*data.all$Tair[i]*data.all$Tair[i]; 
    #   sf.i = sf[1]+ sf[2]*data.all$Tair[i] + sf[3]*data.all$Tair[i]*data.all$Tair[i]; sr.i = sr[1]+ sr[2]*data.all$Tair[i] + sr[3]*data.all$Tair[i]*data.all$Tair[i]
    # }
    # if (no.param == 4) {
    #   k.i = k[1] + k[2]*data.all$Tair[i] + k[3]*data.all$Tair[i]*data.all$Tair[i] + k[4]*data.all$Tair[i]*data.all$Tair[i]*data.all$Tair[i]; 
    #   Y.i = Y[1]+ Y[2]*data.all$Tair[i] + Y[3]*data.all$Tair[i]*data.all$Tair[i] + Y[4]*data.all$Tair[i]*data.all$Tair[i]*data.all$Tair[i]; 
    #   af.i = af[1]+ af[2]*data.all$Tair[i] + af[3]*data.all$Tair[i]*data.all$Tair[i] + af[4]*data.all$Tair[i]*data.all$Tair[i]*data.all$Tair[i]; 
    #   as.i = as[1]+ as[2]*data.all$Tair[i] + as[3]*data.all$Tair[i]*data.all$Tair[i] + as[4]*data.all$Tair[i]*data.all$Tair[i]*data.all$Tair[i]; 
    #   sf.i = sf[1]+ sf[2]*data.all$Tair[i] + sf[3]*data.all$Tair[i]*data.all$Tair[i] + sf[4]*data.all$Tair[i]*data.all$Tair[i]*data.all$Tair[i]; 
    #   sr.i = sr[1]+ sr[2]*data.all$Tair[i] + sr[3]*data.all$Tair[i]*data.all$Tair[i] + sr[4]*data.all$Tair[i]*data.all$Tair[i]*data.all$Tair[i]
    # }
    
    Rm[i] = data.set$Rd.foliage.mean[i-1]*Mleaf[i-1] + data.set$Rd.stem.mean[i-1]*Mwood[i-1]*data.set$SMratio[i-1] + 
      data.set$Rd.branch.mean[i-1]*Mwood[i-1]*data.set$BMratio[i-1] + 
      data.set$Rd.fineroot.mean[i-1]*Mroot[i-1]*data.set$FRratio[i-1] + data.set$Rd.intermediateroot.mean[i-1]*Mroot[i-1]*data.set$IRratio[i-1] + 
      data.set$Rd.coarseroot.mean[i-1]*Mroot[i-1]*data.set$CRratio[i-1] + data.set$Rd.boleroot.mean[i-1]*Mroot[i-1]*data.set$BRratio[i-1]
    
    # Cstorage[i] <- Cstorage[i-1] + data.set$GPP[i-1] - Rm[i-1] - k.i*Cstorage[i-1]
    Cstorage[i] <- (Cstorage[i-1] + rnorm(1, data.set$GPP[i], data.set$GPP_SE[i]) - Rm[i]) / (1 + k.i)
    Sleaf[i] <- Cstorage[i] * tnc.partitioning$foliage[i] 
    Swood[i] <- Cstorage[i] * tnc.partitioning$wood[i] 
    Sroot[i] <- Cstorage[i] * tnc.partitioning$root[i] 
    
    # Mlit[i] <- sf.i*Cleaf[i-1] # foliage litter fall
    # Cleaf[i] <- Cleaf[i-1] + k.i*Cstorage[i-1]*af.i*(1-Y.i) - Mlit[i]
    # Cwood[i] <- Cwood[i-1] + k.i*Cstorage[i-1]*as.i*(1-Y.i)
    # Croot[i] <- Croot[i-1] + k.i*Cstorage[i-1]*(1-af.i-as.i)*(1-Y.i) - sr.i*Croot[i-1]
    
    Cleaf[i] <- Cleaf[i-1] + k.i*Cstorage[i]*af.i*(1-Y.i) / (1 + sf.i)
    Cwood[i] <- Cwood[i-1] + k.i*Cstorage[i]*as.i*(1-Y.i)
    Croot[i] <- Croot[i-1] + k.i*Cstorage[i]*(1-af.i-as.i)*(1-Y.i) / (1 + sr.i)
    
    Mleaf[i] <- Cleaf[i] + Sleaf[i]
    Mwood[i] <- Cwood[i] + Swood[i]
    Mroot[i] <- Croot[i] + Sroot[i]
    
    # Y.modelled[i] = Y.i
    # Rabove[i] = data.set$Rd.foliage.mean[i-1]*Mleaf[i-1] + data.set$Rd.stem.mean[i-1]*Mwood[i-1]*data.set$SMratio[i-1] + data.set$Rd.branch.mean[i-1]*Mwood[i-1]*data.set$BMratio[i-1] +
    #   Y.i*(Mleaf[i]-Mleaf[i-1] + Mwood[i]-Mwood[i-1])
    Rabove[i] = data.set$Rd.foliage.mean[i]*Mleaf[i] + data.set$Rd.stem.mean[i]*Mwood[i]*data.set$SMratio[i] + data.set$Rd.branch.mean[i]*Mwood[i]*data.set$BMratio[i] +
      Y.i*(Mleaf[i]-Mleaf[i-1] + Mwood[i]-Mwood[i-1])
    Mlit[i] <- sf.i*Cleaf[i] # foliage litter fall
  }
  # Rm[i] = data.set$Rd.foliage.mean[i]*Mleaf[i] + data.set$Rd.stem.mean[i]*Mwood[i]*data.set$SMratio[i] + 
  #   data.set$Rd.branch.mean[i]*Mwood[i]*data.set$BMratio[i] + 
  #   data.set$Rd.fineroot.mean[i]*Mroot[i]*data.set$FRratio[i] + data.set$Rd.intermediateroot.mean[i]*Mroot[i]*data.set$IRratio[i] + 
  #   data.set$Rd.coarseroot.mean[i]*Mroot[i]*data.set$CRratio[i] + data.set$Rd.boleroot.mean[i]*Mroot[i]*data.set$BRratio[i]
  
  Mlit = cumsum(Mlit)
  output = data.frame(Cstorage,Mleaf,Mwood,Mroot,Mlit,Sleaf,Rm,Rabove)
  
  return(output)
}
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Defining the model to iteratively calculate Cstorage, Cleaf, Cwood, Croot, Sleaf, Swood, Sroot
model_v2 <- function (no.param,data.set,tnc.partitioning,Y,k,af,as,sf,sr,v1,treat) {
  Mleaf = Mwood = Mroot = Mlit = Y.modelled = Rabove = c()
  Mleaf[1] <- data.set$LM[1]
  Mwood[1] <- data.set$WM[1]
  Mroot[1] <- data.set$RM[1]
  Mlit[1] <- data.set$litter[1]
  
  Cstorage = Sleaf = Swood = Sroot = c()
  
  # From WTC-4 experiment for TNC partitioning to tree organs
  Sleaf[1] = data.set$TNC_leaf[min(which(complete.cases(data.set$TNC_leaf)))]*0.75 # Consider the first available leaf tnc data as the starting point
  Swood[1] = data.set$TNC_wood[min(which(complete.cases(data.set$TNC_wood)))]*0.75 # Consider the first available leaf tnc data as the starting point
  Sroot[1] = data.set$TNC_root[min(which(complete.cases(data.set$TNC_root)))]*0.75 # Consider the first available leaf tnc data as the starting point
  Cstorage[1] = data.set$TNC_tot[min(which(complete.cases(data.set$TNC_tot)))]*0.75 # Consider the first available leaf tnc data as the starting point
  
  # Sleaf[1] = 10
  # Swood[1] = 5
  # Sroot[1] = 5
  # Cstorage[1] = 20
  
  Cleaf <- Croot <- Cwood <- Rm <- c()
  Cleaf[1] <- data.set$LM[1] - Sleaf[1]
  Cwood[1] <- data.set$WM[1] - Swood[1]
  Croot[1] <- data.set$RM[1] - Sroot[1]
  
  Rabove[1] = data.set$Ra[1]
  Rm[1] = data.set$Rd.foliage.mean[1]*Mleaf[1] + data.set$Rd.stem.mean[1]*Mwood[1]*data.set$SMratio[1] +
    data.set$Rd.branch.mean[1]*Mwood[1]*data.set$BMratio[1] +
    data.set$Rd.fineroot.mean[1]*Mroot[1]*data.set$FRratio[1] + data.set$Rd.intermediateroot.mean[1]*Mroot[1]*data.set$IRratio[1] +
    data.set$Rd.coarseroot.mean[1]*Mroot[1]*data.set$CRratio[1] + data.set$Rd.boleroot.mean[1]*Mroot[1]*data.set$BRratio[1]
  
  # Y.modelled = Y[1];
  if (treat[v1] %in% as.factor("amb+ele+const")) {
    k.i = k[1]; Y.i = Y[1];
  }
  if (no.param == 1) {
    if (!treat[v1] %in% as.factor("amb+ele+const")) {
      k.i = k[1]; Y.i = Y[1]; 
    }
    af.i = af[1]; as.i = as[1]; sf.i = sf[1]; sr.i = sr[1]
  }
  for (i in 2:nrow(data.set)) {
    if (no.param == 2) {
      if (!treat[v1] %in% as.factor("amb+ele+const")) {
        k.i = k[1] + k[2]*i; Y.i = Y[1]+ Y[2]*i;
      }
      af.i = af[1]+ af[2]*i; as.i = as[1]+ as[2]*i; sf.i = sf[1]+ sf[2]*i; sr.i = sr[1]+ sr[2]*i
    }
    if (no.param == 3) {
      if (!treat[v1] %in% as.factor("amb+ele+const")) {
        k.i = k[1] + k[2]*i + k[3]*i*i; Y.i = Y[1]+ Y[2]*i + Y[3]*i*i;
      }
      af.i = af[1]+ af[2]*i + af[3]*i*i;
      as.i = as[1]+ as[2]*i + as[3]*i*i; sf.i = sf[1]+ sf[2]*i + sf[3]*i*i; sr.i = sr[1]+ sr[2]*i + sr[3]*i*i
    }
    if (no.param == 4) {
      if (!treat[v1] %in% as.factor("amb+ele+const")) {
        k.i = k[1] + k[2]*i + k[3]*i*i + k[4]*i*i*i; Y.i = Y[1]+ Y[2]*i + Y[3]*i*i + Y[4]*i*i*i; 
      }
      af.i = af[1]+ af[2]*i + af[3]*i*i + af[4]*i*i*i;
      as.i = as[1]+ as[2]*i + as[3]*i*i + as[4]*i*i*i; sf.i = sf[1]+ sf[2]*i + sf[3]*i*i + sf[4]*i*i*i; sr.i = sr[1]+ sr[2]*i + sr[3]*i*i + sr[4]*i*i*i
    }
    # if (no.param == 1) {
    #   k.i = k[1]*data.all$Tair[i]; Y.i = Y[1]*data.all$Tair[i]; af.i = af[1]*data.all$Tair[i]; 
    #   as.i = as[1]*data.all$Tair[i]; sf.i = sf[1]*data.all$Tair[i]; sr.i = sr[1]*data.all$Tair[i]
    # }
    # if (no.param == 2) {
    #   k.i = k[1] + k[2]*data.all$Tair[i]; Y.i = Y[1]+ Y[2]*data.all$Tair[i]; af.i = af[1]+ af[2]*data.all$Tair[i]; 
    #   as.i = as[1]+ as[2]*data.all$Tair[i]; sf.i = sf[1]+ sf[2]*data.all$Tair[i]; sr.i = sr[1]+ sr[2]*data.all$Tair[i]
    # }
    # if (no.param == 3) {
    #   k.i = k[1] + k[2]*data.all$Tair[i] + k[3]*data.all$Tair[i]*data.all$Tair[i]; Y.i = Y[1]+ Y[2]*data.all$Tair[i] + Y[3]*data.all$Tair[i]*data.all$Tair[i]; 
    #   af.i = af[1]+ af[2]*data.all$Tair[i] + af[3]*data.all$Tair[i]*data.all$Tair[i]; as.i = as[1]+ as[2]*data.all$Tair[i] + as[3]*data.all$Tair[i]*data.all$Tair[i]; 
    #   sf.i = sf[1]+ sf[2]*data.all$Tair[i] + sf[3]*data.all$Tair[i]*data.all$Tair[i]; sr.i = sr[1]+ sr[2]*data.all$Tair[i] + sr[3]*data.all$Tair[i]*data.all$Tair[i]
    # }
    # if (no.param == 4) {
    #   k.i = k[1] + k[2]*data.all$Tair[i] + k[3]*data.all$Tair[i]*data.all$Tair[i] + k[4]*data.all$Tair[i]*data.all$Tair[i]*data.all$Tair[i]; 
    #   Y.i = Y[1]+ Y[2]*data.all$Tair[i] + Y[3]*data.all$Tair[i]*data.all$Tair[i] + Y[4]*data.all$Tair[i]*data.all$Tair[i]*data.all$Tair[i]; 
    #   af.i = af[1]+ af[2]*data.all$Tair[i] + af[3]*data.all$Tair[i]*data.all$Tair[i] + af[4]*data.all$Tair[i]*data.all$Tair[i]*data.all$Tair[i]; 
    #   as.i = as[1]+ as[2]*data.all$Tair[i] + as[3]*data.all$Tair[i]*data.all$Tair[i] + as[4]*data.all$Tair[i]*data.all$Tair[i]*data.all$Tair[i]; 
    #   sf.i = sf[1]+ sf[2]*data.all$Tair[i] + sf[3]*data.all$Tair[i]*data.all$Tair[i] + sf[4]*data.all$Tair[i]*data.all$Tair[i]*data.all$Tair[i]; 
    #   sr.i = sr[1]+ sr[2]*data.all$Tair[i] + sr[3]*data.all$Tair[i]*data.all$Tair[i] + sr[4]*data.all$Tair[i]*data.all$Tair[i]*data.all$Tair[i]
    # }
    
    Rm[i] = data.set$Rd.foliage.mean[i-1]*Mleaf[i-1] + data.set$Rd.stem.mean[i-1]*Mwood[i-1]*data.set$SMratio[i-1] + 
      data.set$Rd.branch.mean[i-1]*Mwood[i-1]*data.set$BMratio[i-1] + 
      data.set$Rd.fineroot.mean[i-1]*Mroot[i-1]*data.set$FRratio[i-1] + data.set$Rd.intermediateroot.mean[i-1]*Mroot[i-1]*data.set$IRratio[i-1] + 
      data.set$Rd.coarseroot.mean[i-1]*Mroot[i-1]*data.set$CRratio[i-1] + data.set$Rd.boleroot.mean[i-1]*Mroot[i-1]*data.set$BRratio[i-1]
    
    # Cstorage[i] <- Cstorage[i-1] + data.set$GPP[i-1] - Rm[i-1] - k.i*Cstorage[i-1]
    Cstorage[i] <- (Cstorage[i-1] + rnorm(1, data.set$GPP[i], data.set$GPP_SE[i]) - Rm[i]) / (1 + k.i)
    Sleaf[i] <- Cstorage[i] * tnc.partitioning$foliage[i] 
    Swood[i] <- Cstorage[i] * tnc.partitioning$wood[i] 
    Sroot[i] <- Cstorage[i] * tnc.partitioning$root[i] 
    
    # Mlit[i] <- sf.i*Cleaf[i-1] # foliage litter fall
    # Cleaf[i] <- Cleaf[i-1] + k.i*Cstorage[i-1]*af.i*(1-Y.i) - Mlit[i]
    # Cwood[i] <- Cwood[i-1] + k.i*Cstorage[i-1]*as.i*(1-Y.i)
    # Croot[i] <- Croot[i-1] + k.i*Cstorage[i-1]*(1-af.i-as.i)*(1-Y.i) - sr.i*Croot[i-1]
    
    Cleaf[i] <- Cleaf[i-1] + k.i*Cstorage[i]*af.i*(1-Y.i) / (1 + sf.i)
    Cwood[i] <- Cwood[i-1] + k.i*Cstorage[i]*as.i*(1-Y.i)
    Croot[i] <- Croot[i-1] + k.i*Cstorage[i]*(1-af.i-as.i)*(1-Y.i) / (1 + sr.i)
    
    Mleaf[i] <- Cleaf[i] + Sleaf[i]
    Mwood[i] <- Cwood[i] + Swood[i]
    Mroot[i] <- Croot[i] + Sroot[i]
    
    # Y.modelled[i] = Y.i
    # Rabove[i] = data.set$Rd.foliage.mean[i-1]*Mleaf[i-1] + data.set$Rd.stem.mean[i-1]*Mwood[i-1]*data.set$SMratio[i-1] + data.set$Rd.branch.mean[i-1]*Mwood[i-1]*data.set$BMratio[i-1] +
    #   Y.i*(Mleaf[i]-Mleaf[i-1] + Mwood[i]-Mwood[i-1])
    Rabove[i] = data.set$Rd.foliage.mean[i]*Mleaf[i] + data.set$Rd.stem.mean[i]*Mwood[i]*data.set$SMratio[i] + data.set$Rd.branch.mean[i]*Mwood[i]*data.set$BMratio[i] +
      Y.i*(Mleaf[i]-Mleaf[i-1] + Mwood[i]-Mwood[i-1])
    Mlit[i] <- sf.i*Cleaf[i] # foliage litter fall
  }
  # Rm[i] = data.set$Rd.foliage.mean[i]*Mleaf[i] + data.set$Rd.stem.mean[i]*Mwood[i]*data.set$SMratio[i] + 
  #   data.set$Rd.branch.mean[i]*Mwood[i]*data.set$BMratio[i] + 
  #   data.set$Rd.fineroot.mean[i]*Mroot[i]*data.set$FRratio[i] + data.set$Rd.intermediateroot.mean[i]*Mroot[i]*data.set$IRratio[i] + 
  #   data.set$Rd.coarseroot.mean[i]*Mroot[i]*data.set$CRratio[i] + data.set$Rd.boleroot.mean[i]*Mroot[i]*data.set$BRratio[i]
  
  Mlit = cumsum(Mlit)
  output = data.frame(Cstorage,Mleaf,Mwood,Mroot,Mlit,Sleaf,Rm,Rabove)
  
  return(output)
}
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Defining the model to iteratively calculate Cstorage, Cleaf, Cwood, Croot, Sleaf, Swood, Sroot
model.monthly <- function (data.set,j,tnc.partitioning,Y,k,af,as,sf,sr) {
  Mleaf = Mwood = Mroot = Mlit = Rabove = c()
  Mleaf[1] <- data.set$LM[1]
  Mwood[1] <- data.set$WM[1]
  Mroot[1] <- data.set$RM[1]
  Mlit[1] <- data.set$litter[1]
  
  Cstorage = Sleaf = Swood = Sroot = c()
  
  # From WTC-4 experiment for TNC partitioning to tree organs
  # Leaf TNC C / Leaf C =  ; wood TNC C / wood C =  ; Root TNC C / Root C =  
  Sleaf[1] = data.set$TNC_leaf[min(which(complete.cases(data.set$TNC_leaf)))] # Consider the first available leaf tnc data as the starting point
  Swood[1] = data.set$TNC_wood[min(which(complete.cases(data.set$TNC_wood)))] # Consider the first available leaf tnc data as the starting point
  Sroot[1] = data.set$TNC_root[min(which(complete.cases(data.set$TNC_root)))] # Consider the first available leaf tnc data as the starting point
  Cstorage[1] = data.set$TNC_tot[min(which(complete.cases(data.set$TNC_tot)))] # Consider the first available leaf tnc data as the starting point
  
  Cleaf <- Croot <- Cwood <- Rm <- c()
  Cleaf[1] <- data.set$LM[1] - Sleaf[1]
  Cwood[1] <- data.set$WM[1] - Swood[1]
  Croot[1] <- data.set$RM[1] - Sroot[1]
  
  Rabove[1] = data.set$Ra[1]
  Rm[1] = data.set$Rd.foliage.mean[1]*Mleaf[1] + data.set$Rd.stem.mean[1]*Mwood[1]*data.set$SMratio[1] +
    data.set$Rd.branch.mean[1]*Mwood[1]*data.set$BMratio[1] +
    data.set$Rd.fineroot.mean[1]*Mroot[1]*data.set$FRratio[1] + data.set$Rd.intermediateroot.mean[1]*Mroot[1]*data.set$IRratio[1] +
    data.set$Rd.coarseroot.mean[1]*Mroot[1]*data.set$CRratio[1] + data.set$Rd.boleroot.mean[1]*Mroot[1]*data.set$BRratio[1]
  
  for (i in 2:nrow(data.set)) {
    Rm[i] = data.set$Rd.foliage.mean[i-1]*Mleaf[i-1] + data.set$Rd.stem.mean[i-1]*Mwood[i-1]*data.set$SMratio[i-1] + 
      data.set$Rd.branch.mean[i-1]*Mwood[i-1]*data.set$BMratio[i-1] + 
      data.set$Rd.fineroot.mean[i-1]*Mroot[i-1]*data.set$FRratio[i-1] + data.set$Rd.intermediateroot.mean[i-1]*Mroot[i-1]*data.set$IRratio[i-1] + 
      data.set$Rd.coarseroot.mean[i-1]*Mroot[i-1]*data.set$CRratio[i-1] + data.set$Rd.boleroot.mean[i-1]*Mroot[i-1]*data.set$BRratio[i-1]
    
    Cstorage[i] <- Cstorage[i-1] + data.set$GPP[i-1] - Rm[i-1] - k[(i-1)-(j[i-1])]*Cstorage[i-1]
    # Cstorage[i] <- (Cstorage[i-1] + rnorm(1, data.set$GPP[i], data.set$GPP_SE[i]) - Rm[i]) / (1 + k[(i-1)-(j[i-1])])
    
    Sleaf[i] <- Cstorage[i] * tnc.partitioning$foliage[i] # 33% of storage goes to leaf (WTC-4 experiment)
    Swood[i] <- Cstorage[i] * tnc.partitioning$wood[i] # 53% of storage goes to wood (WTC-4 experiment)
    Sroot[i] <- Cstorage[i] * tnc.partitioning$root[i] # 14% of storage goes to root (WTC-4 experiment)
    
    # Sleaf[i] <- Cstorage[i] * 0.33 # 33% of storage goes to leaf (WTC-4 experiment)
    # Swood[i] <- Cstorage[i] * 0.53 # 53% of storage goes to wood (WTC-4 experiment)
    # Sroot[i] <- Cstorage[i] * 0.14 # 14% of storage goes to root (WTC-4 experiment)
    
    # Mlit[i] <- sf[(i-1)-(j[i-1])]*Cleaf[i-1] # foliage litter fall
    # Cleaf[i] <- Cleaf[i-1] + k[(i-1)-(j[i-1])]*Cstorage[i-1]*af[(i-1)-(j[i-1])]*(1-Y[(i-1)-(j[i-1])]) - Mlit[i]
    # Cwood[i] <- Cwood[i-1] + k[(i-1)-(j[i-1])]*Cstorage[i-1]*as[(i-1)-(j[i-1])]*(1-Y[(i-1)-(j[i-1])])
    # Croot[i] <- Croot[i-1] + k[(i-1)-(j[i-1])]*Cstorage[i-1]*(1-af[(i-1)-(j[i-1])]-as[(i-1)-(j[i-1])])*(1-Y[(i-1)-(j[i-1])]) - sr[(i-1)-(j[i-1])]*Croot[i-1]
    
    Cleaf[i] <- Cleaf[i-1] + k[(i-1)-(j[i-1])]*Cstorage[i]*af[(i-1)-(j[i-1])]*(1-Y[(i-1)-(j[i-1])]) / (1 + sf[(i-1)-(j[i-1])])
    Cwood[i] <- Cwood[i-1] + k[(i-1)-(j[i-1])]*Cstorage[i]*as[(i-1)-(j[i-1])]*(1-Y[(i-1)-(j[i-1])])
    Croot[i] <- Croot[i-1] + k[(i-1)-(j[i-1])]*Cstorage[i]*(1-af[(i-1)-(j[i-1])]-as[(i-1)-(j[i-1])])*(1-Y[(i-1)-(j[i-1])]) / (1 + sr[(i-1)-(j[i-1])])
    
    Mleaf[i] <- Cleaf[i] + Sleaf[i]
    Mwood[i] <- Cwood[i] + Swood[i]
    Mroot[i] <- Croot[i] + Sroot[i]
    
    # Rabove[i] = data.set$Rd.foliage.mean[i]*Mleaf[i] + data.set$Rd.stem.mean[i]*Mwood[i]*data.set$SMratio[i] + data.set$Rd.branch.mean[i]*Mwood[i]*data.set$BMratio[i] +
    #   Y[(i-1)-(j[i-1])]*(Mleaf[i]-Mleaf[i-1] + Mwood[i]-Mwood[i-1])
    Rabove[i] = data.set$Rd.foliage.mean[i]*Mleaf[i] + data.set$Rd.stem.mean[i]*Mwood[i]*data.set$SMratio[i] + data.set$Rd.branch.mean[i]*Mwood[i]*data.set$BMratio[i] +
      Y[(i-1)-(j[i-1])]*(Mleaf[i]-Mleaf[i-1] + Mwood[i]-Mwood[i-1])
    Mlit[i] <- sf[(i-1)-(j[i-1])]*Cleaf[i] # foliage litter fall
  }
  # Rm[i] = data.set$Rd.foliage.mean[i]*Mleaf[i] + data.set$Rd.stem.mean[i]*Mwood[i]*data.set$SMratio[i] + 
  #   data.set$Rd.branch.mean[i]*Mwood[i]*data.set$BMratio[i] + 
  #   data.set$Rd.fineroot.mean[i]*Mroot[i]*data.set$FRratio[i] + data.set$Rd.intermediateroot.mean[i]*Mroot[i]*data.set$IRratio[i] + 
  #   data.set$Rd.coarseroot.mean[i]*Mroot[i]*data.set$CRratio[i] + data.set$Rd.boleroot.mean[i]*Mroot[i]*data.set$BRratio[i]
  
  Mlit = cumsum(Mlit)
  output = data.frame(Cstorage,Mleaf,Mwood,Mroot,Mlit,Sleaf,Rm,Rabove)
  
  return(output)
}
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

#-- Function to run Carbon Balance Model (CBM) without considering storage pool. 
#-----------------------------------------------------------------------------------------
# This script define the model equations to carry out Bayesian calibration for 
# 4 variables (allocation fractions: "Y",af","as","sf") on 
# various temporal scales to estimate Carbon pools (Cleaf,Cwood,Croot)
#-----------------------------------------------------------------------------------------
# This version does not consider storage pool

# Defining the model to iteratively calculate Cleaf, Cwood, Croot
model.without.storage <- function (no.param,data.set,Y,af,as,sf,sr) {
  Mleaf = Mwood = Mroot = Mlit = Rm = Y.modelled = Rabove = c()
  Mleaf[1] <- data.set$LM[1]
  Mwood[1] <- data.set$WM[1]
  Mroot[1] <- data.set$RM[1]
  Mlit[1] <- data.set$litter[1]
  
  Rabove[1] = data.set$Ra[1]
  Rm[1] = data.set$Rd.foliage.mean[1]*Mleaf[1] + data.set$Rd.stem.mean[1]*Mwood[1]*data.set$SMratio[1] +
    data.set$Rd.branch.mean[1]*Mwood[1]*data.set$BMratio[1] +
    data.set$Rd.fineroot.mean[1]*Mroot[1]*data.set$FRratio[1] + data.set$Rd.intermediateroot.mean[1]*Mroot[1]*data.set$IRratio[1] +
    data.set$Rd.coarseroot.mean[1]*Mroot[1]*data.set$CRratio[1] + data.set$Rd.boleroot.mean[1]*Mroot[1]*data.set$BRratio[1]
  
  # Y.modelled = Y[1];
  if (no.param == 1) {
    Y.i = Y[1]; af.i = af[1]; as.i = as[1]; sf.i = sf[1]; sr.i = sr[1]
  }
  
  for (i in 2:nrow(data.set)) {
    if (no.param == 2) {
      Y.i = Y[1]+ Y[2]*i; af.i = af[1]+ af[2]*i; as.i = as[1]+ as[2]*i; sf.i = sf[1]+ sf[2]*i; sr.i = sr[1]+ sr[2]*i
    }
    if (no.param == 3) {
      Y.i = Y[1]+ Y[2]*i + Y[3]*i*i; af.i = af[1]+ af[2]*i + af[3]*i*i; as.i = as[1]+ as[2]*i + as[3]*i*i; sf.i = sf[1]+ sf[2]*i + sf[3]*i*i;
      sr.i = sr[1]+ sr[2]*i + sr[3]*i*i
    }
    if (no.param == 4) {
      Y.i = Y[1]+ Y[2]*i + Y[3]*i*i + Y[4]*i*i*i; af.i = af[1]+ af[2]*i + af[3]*i*i + af[4]*i*i*i; as.i = as[1]+ as[2]*i + as[3]*i*i + as[4]*i*i*i; 
      sf.i = sf[1]+ sf[2]*i + sf[3]*i*i + sf[4]*i*i*i; sr.i = sr[1]+ sr[2]*i + sr[3]*i*i + sr[4]*i*i*i
    }
    Rm[i] = data.set$Rd.foliage.mean[i-1]*Mleaf[i-1] + data.set$Rd.stem.mean[i-1]*Mwood[i-1]*data.set$SMratio[i-1] + 
      data.set$Rd.branch.mean[i-1]*Mwood[i-1]*data.set$BMratio[i-1] + 
      data.set$Rd.fineroot.mean[i-1]*Mroot[i-1]*data.set$FRratio[i-1] + data.set$Rd.intermediateroot.mean[i-1]*Mroot[i-1]*data.set$IRratio[i-1] + 
      data.set$Rd.coarseroot.mean[i-1]*Mroot[i-1]*data.set$CRratio[i-1] + data.set$Rd.boleroot.mean[i-1]*Mroot[i-1]*data.set$BRratio[i-1]
    
    Mlit[i] <- sf.i*Mleaf[i-1] # foliage litter fall
    Mleaf[i] <- Mleaf[i-1] + (data.all$GPP[i-1] - Rm[i-1]) * af.i*(1-Y.i) - Mlit[i]
    Mwood[i] <- Mwood[i-1] + (data.all$GPP[i-1] - Rm[i-1]) * as.i*(1-Y.i)
    Mroot[i] <- Mroot[i-1] + (data.all$GPP[i-1] - Rm[i-1]) * (1-af.i-as.i)*(1-Y.i) - sr.i*Mroot[i-1]
    
    # Y.modelled[i] = Y.i
    Rabove[i] = data.set$Rd.foliage.mean[i-1]*Mleaf[i-1] + data.set$Rd.stem.mean[i-1]*Mwood[i-1]*data.set$SMratio[i-1] + data.set$Rd.branch.mean[i-1]*Mwood[i-1]*data.set$BMratio[i-1] +
      Y.i*(Mleaf[i]-Mleaf[i-1] + Mwood[i]-Mwood[i-1])
  }
  Mlit = cumsum(Mlit)
  output = data.frame(Mleaf,Mwood,Mroot,Mlit,Rm,Rabove)
  return(output)
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# This version does not consider storage pool

# Defining the model to iteratively calculate Cleaf, Cwood, Croot
model.without.storage.monthly <- function (data.set,j,Y,af,as,sf,sr) {
  Mleaf = Mwood = Mroot = Mlit = Rm = Rabove = c()
  Mleaf[1] <- data.set$LM[1]
  Mwood[1] <- data.set$WM[1]
  Mroot[1] <- data.set$RM[1]
  Mlit[1] <- data.set$litter[1]
  
  Rabove[1] = data.set$Ra[1]
  Rm[1] = data.set$Rd.foliage.mean[1]*Mleaf[1] + data.set$Rd.stem.mean[1]*Mwood[1]*data.set$SMratio[1] +
    data.set$Rd.branch.mean[1]*Mwood[1]*data.set$BMratio[1] +
    data.set$Rd.fineroot.mean[1]*Mroot[1]*data.set$FRratio[1] + data.set$Rd.intermediateroot.mean[1]*Mroot[1]*data.set$IRratio[1] +
    data.set$Rd.coarseroot.mean[1]*Mroot[1]*data.set$CRratio[1] + data.set$Rd.boleroot.mean[1]*Mroot[1]*data.set$BRratio[1]
  
  if (no.param == 1) {
    Y.i = Y[1]; af.i = af[1]; as.i = as[1]; sf.i = sf[1]; sr.i = sr[1]
  }
  
  for (i in 2:length(GPP)) {
    if (no.param == 2) {
      Y.i = Y[1]+ Y[2]*i; af.i = af[1]+ af[2]*i; as.i = as[1]+ as[2]*i; sf.i = sf[1]+ sf[2]*i; sr.i = sr[1]+ sr[2]*i
    }
    if (no.param == 3) {
      Y.i = Y[1]+ Y[2]*i + Y[3]*i*i; af.i = af[1]+ af[2]*i + af[3]*i*i; 
      as.i = as[1]+ as[2]*i + as[3]*i*i; 
      sf.i = sf[1]+ sf[2]*i + sf[3]*i*i;
      sr.i = sr[1]+ sr[2]*i + sr[3]*i*i
    }
    Rm[i] = data.set$Rd.foliage.mean[i-1]*Mleaf[i-1] + data.set$Rd.stem.mean[i-1]*Mwood[i-1]*data.set$SMratio[i-1] + 
      data.set$Rd.branch.mean[i-1]*Mwood[i-1]*data.set$BMratio[i-1] + 
      data.set$Rd.fineroot.mean[i-1]*Mroot[i-1]*data.set$FRratio[i-1] + data.set$Rd.intermediateroot.mean[i-1]*Mroot[i-1]*data.set$IRratio[i-1] + 
      data.set$Rd.coarseroot.mean[i-1]*Mroot[i-1]*data.set$CRratio[i-1] + data.set$Rd.boleroot.mean[i-1]*Mroot[i-1]*data.set$BRratio[i-1]
    
    Mlit[i] <- sf[(i-1)-(j[i-1])]*Mleaf[i-1] # foliage litter fall
    Mleaf[i] <- Mleaf[i-1] + (data.all$GPP[i-1] - Rm[i-1]) *af[(i-1)-(j[i-1])]*(1-Y[(i-1)-(j[i-1])]) - Mlit[i]
    Mwood[i] <- Mwood[i-1] + (data.all$GPP[i-1] - Rm[i-1]) *as[(i-1)-(j[i-1])]*(1-Y[(i-1)-(j[i-1])])
    Mroot[i] <- Mroot[i-1] + (data.all$GPP[i-1] - Rm[i-1]) *(1-af[(i-1)-(j[i-1])]-as[(i-1)-(j[i-1])])*(1-Y[(i-1)-(j[i-1])]) - sr[(i-1)-(j[i-1])]*Mroot[i-1]
    
    Rabove[i] = data.set$Rd.foliage.mean[i]*Mleaf[i] + data.set$Rd.stem.mean[i]*Mwood[i]*data.set$SMratio[i] + data.set$Rd.branch.mean[i]*Mwood[i]*data.set$BMratio[i] +
      Y[(i-1)-(j[i-1])]*(Mleaf[i]-Mleaf[i-1] + Mwood[i]-Mwood[i-1])
  }
  # Rm[i] = data.set$Rd.foliage.mean[i]*Mleaf[i] + data.set$Rd.stem.mean[i]*Mwood[i]*data.set$SMratio[i] + 
  #   data.set$Rd.branch.mean[i]*Mwood[i]*data.set$BMratio[i] + 
  #   data.set$Rd.fineroot.mean[i]*Mroot[i]*data.set$FRratio[i] + data.set$Rd.intermediateroot.mean[i]*Mroot[i]*data.set$IRratio[i] + 
  #   data.set$Rd.coarseroot.mean[i]*Mroot[i]*data.set$CRratio[i] + data.set$Rd.boleroot.mean[i]*Mroot[i]*data.set$BRratio[i]
  
  Mlit = cumsum(Mlit)
  output = data.frame(Mleaf,Mwood,Mroot,Mlit,Rm,Rabove)
  return(output)
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

################ Figure 4 #####################
# Plot modelled parameters with 3 Grouped treatments and quadratic parameter setting
#-------------------------------------------------------------------------------------
plot.Modelled.parameters.wtc3 <- function(result,with.storage) {
  # listOfDataFrames <- vector(mode = "list", length = 2)
  # for (i in 1:2) {
  #   listOfDataFrames[[i]] <- data.frame(result[[i]][[1]])
  # }
  # no.param.par.var = do.call("rbind", listOfDataFrames)
  # names(no.param.par.var) = "no.param"
  # 
  # listOfDataFrames <- vector(mode = "list", length = 2)
  # for (i in 1:2) {
  #   listOfDataFrames[[i]] <- data.frame(result[[i]][[2]])
  # }
  # summary.param = do.call("rbind", listOfDataFrames)
  
  # cbPalette = c("gray", "skyblue", "orange", "green3", "yellow3", "#0072B2", "#D55E00")
  cbPalette = c("cyan", "firebrick", "darkorange", "deepskyblue3")
  # cbPalette = c("cyan", "darkorange")
  i = 0
  font.size = 10
  plot = list() 
  if (with.storage==T) { 
    var = as.factor(c("k","Y","af","as","ar","sf","sr"))
    title = as.character(c("A","B","C","D","E","F","G"))
  } else {
    var = as.factor(c("Y","af","as","ar","sf","sr"))
    title = as.character(c("A","B","C","D","E","F"))
  }
  pd <- position_dodge(0.5)
  no.param.par.var = result[[1]][[1]]
  summary.param = result[[1]][[2]]
  summary.data = result[[1]][[3]]
  summary.output = result[[1]][[4]]
  summary.error = result[[1]][[5]]
  
  for (p in 1:length(var)) {
    summary.param.set.limit = subset(summary.param, variable %in% var[p])
    for (z in 1:length(no.param.par.var)) {
      summary.param.set = subset(summary.param, variable %in% var[p] & no.param %in% no.param.par.var[z])
      summary.param.set$treatment = unlist(summary.param.set$treatment)
      i = i + 1
      plot[[i]] = ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = treatment, colour=factor(treatment))) +
        geom_ribbon(data = summary.param.set, aes(ymin=Parameter-Parameter_SD, ymax=Parameter+Parameter_SD), linetype=2, alpha=0.1,size=0.1) +
        geom_point(position=pd,size=0.01) +
        geom_line(position=pd,data = summary.param.set, aes(x = Date, y = Parameter,  group = treatment, colour=factor(treatment)),size=1) +
        # ylab(paste(as.character(var[p]),"(fraction)")) +
        ylab(paste(as.character(var[p]))) +
        labs(colour="Treatment") +
        scale_color_manual(values=cbPalette[1:4]) +
        scale_y_continuous(limits = c(min(summary.param.set.limit$Parameter)-2*max(summary.param.set.limit$Parameter_SD),
                                      max(summary.param.set.limit$Parameter)+2*max(summary.param.set.limit$Parameter_SD))) +
        annotate("text", x = min(summary.param.set$Date), y = max(summary.param.set$Parameter) + 2*max(summary.param.set$Parameter_SD), size = font.size-7, label = paste(title[p])) +
        theme_bw() +
        theme(legend.title = element_text(colour="black", size=font.size)) +
        theme(legend.text = element_text(colour="black", size=font.size-3)) +
        theme(legend.position = c(0.65,0.85)) +
        theme(legend.key = element_blank()) +
        theme(text = element_text(size=font.size)) +
        theme(axis.title.x = element_blank()) +
        theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      if (with.storage==T) { 
        if (p==1) {
          # plot[[i]] = plot[[i]] + scale_colour_manual(name="", breaks=c("1", "2", "3"),
          #                                             labels=c("Small", "Large", "Free"), values=cbPalette[2:4]) +
          plot[[i]] = plot[[i]] + scale_colour_manual(name="", values=cbPalette[1:4]) +
            ylab(expression(k~"(g C "*g^"-1"*" C "*d^"-1"*")"))
          plot[[i]] = plot[[i]] + theme(legend.key.height=unit(0.7,"line"))
        } else if (p>1) {
          plot[[i]] = plot[[i]] + guides(colour=FALSE)
        } 
        if (p==2) {
          plot[[i]] = plot[[i]] + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 0.8), units="line"))
        }
        if (p==3) {
          # plot[[i]] = plot[[i]] + ylab(expression(a[f]~"(fraction)")) + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 1), units="line"))
          plot[[i]] = plot[[i]] + ylab(expression(a[f])) + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 1), units="line"))
        }
        if (p==4) {
          plot[[i]] = plot[[i]] + ylab(expression(a[w])) + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 1), units="line"))
        }
        if (p==5) {
          plot[[i]] = plot[[i]] + ylab(expression(a[r])) + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 1), units="line"))
        }
        if (p==6) {
          plot[[i]] = plot[[i]] + ylab(expression(s[f]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
        }
        if (p==7) {
          plot[[i]] = plot[[i]] + ylab(expression(s[r]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
        }
        
      } else {
        if (p==1) {
          plot[[i]] = plot[[i]] + scale_colour_manual(name="", values=cbPalette[1:4])
          plot[[i]] = plot[[i]] + theme(legend.key.height=unit(0.7,"line"))
        } else if (p>1) {
          plot[[i]] = plot[[i]] + guides(colour=FALSE)
        } 
        if (p==2) {
          plot[[i]] = plot[[i]] + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 0.8), units="line"))
          plot[[i]] = plot[[i]] + ylab(expression(a[f])) + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 1), units="line"))
        }
        if (p==3) {
          plot[[i]] = plot[[i]] + ylab(expression(a[w])) + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 1), units="line"))
        }
        if (p==4) {
          plot[[i]] = plot[[i]] + ylab(expression(a[r])) + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 1), units="line"))
        }
        if (p==5) {
          plot[[i]] = plot[[i]] + ylab(expression(s[f]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
        }
        if (p==6) {
          plot[[i]] = plot[[i]] + ylab(expression(s[r]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
        }
      }
    }
  }
  
  png("output/Figure_4_modelled_parameters.png", units="px", width=2000, height=2000, res=250)
  print (do.call(grid.arrange,  plot))
  dev.off()
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

################ Figure 5 #####################
# Plot Daily analysis (lines) with optimum parameter setting and intermittent observations (symbols) of selected carbon stocks
#-------------------------------------------------------------------------------------
plot.Modelled.biomass.wtc3 <- function(result,with.storage) { 
  # listOfDataFrames <- vector(mode = "list", length = 2)
  # for (i in 1:2) {
  #   listOfDataFrames[[i]] <- data.frame(result[[i]][[4]])
  # }
  # summary.output = do.call("rbind", listOfDataFrames)
  # 
  # listOfDataFrames <- vector(mode = "list", length = 2)
  # for (i in 1:2) {
  #   listOfDataFrames[[i]] <- data.frame(result[[i]][[5]])
  # }
  # summary.error = do.call("rbind", listOfDataFrames)
  
  # cbPalette = c("gray", "skyblue", "orange", "green3", "yellow3", "#0072B2", "#D55E00")
  cbPalette = c("cyan", "firebrick", "darkorange", "deepskyblue3")
  # cbPalette = c("cyan", "darkorange")
  i = 0
  font.size = 10
  plot = list() 
  no.param.par.var = result[[1]][[1]]
  summary.param = result[[1]][[2]]
  summary.data = result[[1]][[3]]
  summary.output = result[[1]][[4]]
  summary.error = result[[1]][[5]]
  if (with.storage==T) { 
    meas = as.factor(c("LM","WM","RM","litter","TNC_leaf","Ra"))
    res = as.factor(c("Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Sleaf.modelled","Rabove"))
    error = as.factor(c("LM_SE","WM_SE","RM_SE","litter_SE","TNC_leaf_SE","Ra_SE"))
    title = as.character(c("A","B","C","D","E","F"))
  } else {
    meas = as.factor(c("LM","WM","RM","litter","Ra"))
    res = as.factor(c("Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Rabove"))
    error = as.factor(c("LM_SE","WM_SE","RM_SE","litter_SE","Ra_SE"))
    title = as.character(c("A","B","C","D","E"))
  }
  pd <- position_dodge(2) # move the overlapped errorbars horizontally
  for (p in 1:length(meas)) {
    # summary.data.Cpool = subset(summary.data,variable %in% meas[p])
    summary.output.Cpool = subset(summary.output,variable %in% res[p])
    summary.error.Cpool = subset(summary.error,variable %in% error[p])
    summary.output.Cpool$treat.type = as.factor(ifelse(summary.output.Cpool$treat.type %in% 1, ("Individual"), ("Combined")))
    summary.error.Cpool$treat.type = as.factor(ifelse(summary.error.Cpool$treat.type %in% 1, ("Individual"), ("Combined")))
    
    i = i + 1
    if (meas[p]=="Ra") {
      plot[[i]] = ggplot(summary.error.Cpool, aes(x=Date, y=parameter, group = treatment, colour=treatment)) + 
        geom_point(position=pd,size=0.3) +
        geom_ribbon(data = summary.error.Cpool, aes(ymin=parameter-value, ymax=parameter+value), linetype=2, alpha=0.1,size=0.1) +
        # geom_errorbar(position=pd,aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=0.5) +
        # geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,volume.group,no.param), linetype=volume.group, colour=volume, size=no.param)) +
        # geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,no.param), linetype=no.param, colour=volume)) +
        geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(treatment,treat.type), colour=treatment, linetype=treat.type)) +
        ylab(paste(as.character(meas[p]),"(g C)")) + xlab("Month") +
        # ggtitle("C pools - Measured (points) vs Modelled (lines)") +
        # labs(colour="Soil Volume", linetype="Grouping treatment", size="Total No of Parameter") +
        # labs(colour="Pot Volume (L)", linetype="No. of Parameters") +
        labs(colour="Treatment",linetype="Parameter option") +
        scale_color_manual(values=cbPalette[1:4]) +
        # scale_color_manual(labels = c("Individuals", "One Group"), values = c("blue", "red")) +
        # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g C plant-1)")) +
        theme_bw() +
        annotate("text", x = max(summary.output.Cpool$Date), y = min(summary.output.Cpool$value), size = font.size-7, label = paste(title[p])) +
        theme(legend.title = element_text(colour="black", size=font.size-2)) +
        theme(legend.text = element_text(colour="black", size = font.size-3)) +
        theme(legend.key.height=unit(0.6,"line")) +
        theme(legend.position = c(0.4,0.8)) + theme(legend.box = "horizontal") + 
        theme(legend.key = element_blank()) +
        theme(text = element_text(size=font.size)) +
        theme(axis.title.x = element_blank()) +
        theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
        # theme(plot.title = element_text(hjust = 0)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
      
    } else {
      plot[[i]] = ggplot(summary.error.Cpool, aes(x=Date, y=parameter, group = treatment, colour=treatment)) + 
        geom_point(position=pd) +
        geom_errorbar(position=pd,aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=0.5) +
        # geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,volume.group,no.param), linetype=volume.group, colour=volume, size=no.param)) +
        # geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,no.param), linetype=no.param, colour=volume)) +
        geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(treatment,treat.type), colour=treatment, linetype=treat.type)) +
        ylab(paste(as.character(meas[p]),"(g C)")) + xlab("Month") +
        labs(colour="Treatment",linetype="Parameter option") +
        scale_color_manual(values=cbPalette[1:4]) +
        # scale_color_manual(labels = c("Individuals", "One Group"), values = c("blue", "red")) +
        # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g C plant-1)")) +
        theme_bw() +
        annotate("text", x = max(summary.output.Cpool$Date), y = min(summary.output.Cpool$value), size = font.size-7, label = paste(title[p])) +
        # theme(plot.title = element_text(size = 20, face = "bold")) +
        theme(legend.title = element_text(colour="black", size=font.size-2)) +
        theme(legend.text = element_text(colour="black", size = font.size-3)) +
        theme(legend.key.height=unit(0.6,"line")) +
        theme(legend.position = c(0.4,0.8)) + theme(legend.box = "horizontal") + 
        theme(legend.key = element_blank()) +
        theme(text = element_text(size=font.size)) +
        theme(axis.title.x = element_blank()) +
        theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
        # theme(plot.title = element_text(hjust = 0)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    }
    if (with.storage==T) {
      if (p==1) {
        plot[[i]] = plot[[i]] + ylab(expression(C["t,f"]~"(g C "*plant^"-1"*")"))
        # plot[[i]] = plot[[i]]  + theme(legend.key.height=unit(0.6,"line"))
      } else if (p==2) {
        plot[[i]] = plot[[i]] + ylab(expression(C["t,w"]~"(g C "*plant^"-1"*")"))
        # plot[[i]] = plot[[i]] + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 0.75), units="line"))
      } else if (p==3) {
        plot[[i]] = plot[[i]] + ylab(expression(C["t,r"]~"(g C "*plant^"-1"*")"))
      } else if (p==4) {
        plot[[i]] = plot[[i]] + ylab(expression(C["f,lit"]~"(g C "*plant^"-1"*")"))
      } else if (p==5) {
        plot[[i]] = plot[[i]] + ylab(expression(C["n,f"]~"(g C "*plant^"-1"*")"))
      } else {
        plot[[i]] = plot[[i]] + ylab(expression(R[a]~"(g C "*plant^"-1"*")"))
      }
    } else {
      if (p==1) {
        plot[[i]] = plot[[i]] + ylab(expression(C["t,f"]~"(g C "*plant^"-1"*")"))
        # plot[[i]] = plot[[i]]  + theme(legend.key.height=unit(0.6,"line"))
      } else if (p==2) {
        plot[[i]] = plot[[i]] + ylab(expression(C["t,w"]~"(g C "*plant^"-1"*")"))
        # plot[[i]] = plot[[i]] + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 0.75), units="line"))
      } else if (p==3) {
        plot[[i]] = plot[[i]] + ylab(expression(C["t,r"]~"(g C "*plant^"-1"*")"))
      } else if (p==4) {
        plot[[i]] = plot[[i]] + ylab(expression(C["f,lit"]~"(g C "*plant^"-1"*")"))
      } else {
        plot[[i]] = plot[[i]] + ylab(expression(R[a]~"(g C "*plant^"-1"*")"))
      }
    }
    if (p!=4) {
      plot[[i]] = plot[[i]] + guides(colour=FALSE,linetype=FALSE)
    }
    
    # #----------------------------------------------------------------------------------------------------------------
    # # keeps <- c("Date", "volume", "tnc.conc", "tnc.conc_SE")
    # # tnc.data = tnc.data.processed[ , keeps, drop = FALSE]
    # 
    # if (p == 4) {
    #   plot[[i]] = ggplot(summary.error.Cpool, aes(x=Date, y=parameter, group = volume, colour=volume)) +
    #     geom_point(position=pd) +
    #     geom_errorbar(position=pd,aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=2) +
    #     geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = volume, colour=volume)) +
    #     ylab(paste(as.character(meas[p]),"(g C)")) + xlab("Month") +
    #     labs(colour="Pot Volume (L)") +
    #     theme_bw() +
    #     annotate("text", x = min(summary.output.Cpool$Date), y = max(summary.output.Cpool$value), size = font.size-7, label = paste(title[p])) +
    #     theme(legend.title = element_text(colour="black", size=font.size)) +
    #     theme(legend.text = element_text(colour="black", size = font.size)) +
    #     theme(legend.position = c(0.17,0.7)) +
    #     theme(legend.key = element_blank()) +
    #     theme(text = element_text(size=font.size)) +
    #     theme(axis.title.x = element_blank()) +
    #     theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    #     ylab(expression(S[leaf]~"(% of"~M[leaf]~")")) + guides(colour=FALSE)
    # }
    # #----------------------------------------------------------------------------------------------------------------
    
  }
  
  png("output/Figure_5_modelled_biomass.png", units="px", width=1600, height=1300, res=220)
  print (do.call(grid.arrange,  plot))
  dev.off()
  
  # # #----------------------------------------------------------------------------------------------------------------
  # # # Represent Sleaf as a concentration of Mleaf instead of total mass
  # if (p == 4) {
  #   summary.output.Mleaf = subset(summary.output,variable %in% "Mleaf.modelled")
  #   summary.output.Sleaf = subset(summary.output,variable %in% "Sleaf.modelled")
  #   summary.error.Sleaf = subset(summary.error,variable %in% "Sleaf_SD")
  #   summary.output.Sleaf$value = summary.output.Sleaf$value / summary.output.Mleaf$value * 100
  #   summary.output.Sleaf = summary.output.Sleaf[,-c(5,6)]
  #   
  #   # summary.error.Sleaf$value = summary.error.Sleaf$value / lm.daily.m$leafmass * 100
  #   leafmass.daily = read.csv("processed_data/Cleaf_daily_data.csv") # Unit gC per gC plant
  #   leafmass.daily = leafmass.daily[with(leafmass.daily, order(volume,Date)), ]
  #   summary.error.Sleaf$value = ((summary.error.Sleaf$value*summary.error.Sleaf$value + leafmass.daily$leafmass_SE*leafmass.daily$leafmass_SE)/2)^0.5 / lm.daily.m$leafmass * 100
  #   summary.error.Sleaf$parameter = summary.error.Sleaf$parameter / lm.daily.m$leafmass * 100
  #   summary.error.Sleaf = summary.error.Sleaf[,-c(6,7)]
  #   
  #   pd <- position_dodge(4) # move the overlapped errorbars horizontally
  #   plot[[i]] = ggplot(summary.error.Sleaf, aes(x=Date, y=parameter, group = volume, colour=volume)) +
  #     geom_errorbar(position=pd,aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=0.2) +
  #     geom_line(position=pd,data = summary.output.Sleaf, aes(x = Date, y = value, group = volume, colour=volume)) +
  #     geom_point(position=pd) +
  #     # ylab("Sleaf (g C)") + xlab("Month") +
  #     ylab(paste(as.character(meas[p]),"(g C)")) + xlab("Month") +
  #     labs(colour="Pot Volume (L)") +
  #     theme_bw() +
  #     annotate("text", x = min(summary.output.Sleaf$Date), y = max(summary.output.Sleaf$value), size = font.size-7, label = paste(title[p])) +
  #     theme(legend.title = element_text(colour="black", size=font.size)) +
  #     theme(legend.text = element_text(colour="black", size = font.size)) +
  #     theme(legend.position = c(0.17,0.7)) +
  #     theme(legend.key = element_blank()) +
  #     theme(text = element_text(size=font.size)) +
  #     theme(axis.title.x = element_blank()) +
  #     theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #     ylab(expression(S[leaf]~"(% of"~M[leaf]~")")) + guides(colour=FALSE)
  # }
  # 
  # png("output/Figure_5_modelled_biomass_Sleaf_conc.png", units="px", width=2200, height=1600, res=220)
  # print (do.call(grid.arrange,  plot))
  # dev.off()
  # # #----------------------------------------------------------------------------------------------------------------
  
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
