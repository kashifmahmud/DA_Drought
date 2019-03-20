# Manual calculation of parameters

#-------------------------------------------------------------------------------------
# Clear the workspace (if needed)
rm(list=ls())
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# Load required libraries. There are quite a few, including some non-standard functions that are not on CRAN. 
# This script will check for required libraries and install any that are missing	
source('R/load_packages_wtc3.R')	

# Load the custom analysis and plotting functions that do all of the actual work	
source("R/functions_wtc3.R")	
source("R/functions_wtc3_CBM.R")	

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
CBM.wtc3.manual <- function(no.param.par.var, treat.group, with.storage, model.comparison, model.optimization, param.set) {
  
  if (with.storage==T) {
    no.var = 6 # variables to be modelled are: k,Y,af,as,sf,sr
  } else {
    no.var = 5 # variables to be modelled are: Y,af,as,sf,sr
  }
  
  aic.bic = data.frame(matrix(ncol = 5, nrow = length(no.param.par.var)*length(treat.group)))
  names(aic.bic) <- c("logLi","aic","bic","no.param","treatment")
  q = 0 # Indicates the iteration number
  
  for (v in 1:length(treat.group)) {
    data.set = subset(data.all,(Treatment %in% treat.group[v]))
    
    q = q + 1
    j = c()
    j[1] = 0
    i = seq(1,nrow(data.set),1)
    j[i] = i - ceiling(i/param.vary)*1  # j is for parameter settings for various time frames
    
    param = param.set[[v]]
    
    output.set = model.monthly.manual(data.set,j,tnc.partitioning,param)
    
    # listofdfs <- list()
    # for (u in 1:10) {
    #   output.set = model.monthly.manual(data.set,j,tnc.partitioning,param)
    #   listofdfs[[u]] <- output.set
    # }
    # output.set = aaply(laply(listofdfs, as.matrix), c(2, 3), mean)
    
    output = as.data.frame(output.set)
    # output$treatment = as.factor(treat.group[v])
    
    param$ar = 1 - param$af - param$as
    param$Date = data.set$Date[seq(1,nrow(data.set),param.vary)]
    melted.param = melt(param[,c("k","Y","af","as","ar","sf","sr","Date")], id.vars="Date")
    melted.param = data.frame(melted.param$Date, melted.param$variable, melted.param$value)
    names(melted.param) = c("Date","variable","Parameter")
    melted.param$Date = as.Date(melted.param$Date)
    melted.param$treatment = treat.group[v]
    melted.param$no.param = as.factor(ceiling(no.param.par.var))
    
    # Measured (data) vs Modelled Plant Carbon pools for plotting and comparison
    output$Date = data.set$Date
    if (with.storage==T) { 
      names(output) = c("Cstorage.modelled","Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Sleaf.modelled","Rm","Rabove","Date")
      melted.output = melt(output[,c("Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Sleaf.modelled","Rm","Rabove","Date")], id.vars="Date")
      melted.Cstorage = output[,c("Cstorage.modelled","Date")]
      melted.data = melt(data.set[ , c("LM","WM","RM","litter","TNC_leaf","Ra","Date")], id.vars="Date")
      melted.error = melt(data.set[ , c("LM_SE","WM_SE","RM_SE","litter_SE","TNC_leaf_SE","Ra_SE","Date")], id.vars="Date")
      melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
      melted.Cstorage$treatment = as.factor(treat.group[v])
      melted.Cstorage$no.param = as.factor(ceiling(no.param.par.var))
    } else {
      names(output) = c("Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Rm","Rabove","Date")
      melted.output = melt(output[,c("Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Rm","Rabove","Date")], id.vars="Date")
      melted.data = melt(data.set[ , c("LM","WM","RM","litter","Ra","Date")], id.vars="Date")
      melted.error = melt(data.set[ , c("LM_SE","WM_SE","RM_SE","litter_SE","Ra_SE","Date")], id.vars="Date")
    }
    melted.data$Date = as.Date(melted.data$Date)
    melted.error$Date = as.Date(melted.error$Date)
    melted.error$treatment = as.factor(treat.group[v])
    melted.error$parameter = melted.data$value
    melted.output$Date = as.Date(melted.output$Date)
    melted.data$treatment = as.factor(treat.group[v])
    melted.output$treatment = as.factor(treat.group[v])
    melted.output$no.param = as.factor(ceiling(no.param.par.var))
    
    # Storing the summary of data, outputs, Cstorage, parameters
    if (q == 1) {
      summary.data = melted.data
      summary.error = melted.error
      summary.output = melted.output
      if (with.storage==T) {
        summary.Cstorage = melted.Cstorage
      }
      summary.param = melted.param
    }
    if (q > 1) {
      summary.output = rbind(summary.output,melted.output)
      if (with.storage==T) {
        summary.Cstorage = rbind(summary.Cstorage,melted.Cstorage)
      }
      summary.param = rbind(summary.param,melted.param)
      summary.error = rbind(summary.error,melted.error)
      summary.data = rbind(summary.data,melted.data)
    }
    
    # Calculate LogLi, AIC, BIC, Time to find the most accurate model for best balance between model fit and complexity
    output = output
    if (with.storage==T) { 
      names(output) = c("Cstorage","Mleaf","Mwood","Mroot","Mlit","Sleaf","Rm","Rabove","Date") # Rename for the logLikelihood function
    } else {
      names(output) = c("Mleaf","Mwood","Mroot","Mlit","Rm","Rabove","Date")
    }
    
    aic.bic[q,1] <- logLikelihood.wtc3.final(no.param.par.var,data.set,output,with.storage,model.comparison) # Calculate logLikelihood
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
    
    aic.bic[q,4] = no.param
    aic.bic[q,5] = as.character(treat.group[v])
  }
  
  bic = data.frame(aic.bic[,c("bic","no.param","treatment")])
  
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

# Defining the model to iteratively calculate Cstorage, Cleaf, Cwood, Croot, Sleaf, Swood, Sroot
model.monthly.manual <- function (data.set,j,tnc.partitioning,param) {
  Y = param$Y; k = param$k; af = param$af; as = param$as; sf = param$sf; sr = param$sr
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

################ Figure 4 #####################
# Plot modelled parameters with 3 Grouped treatments and quadratic parameter setting
#-------------------------------------------------------------------------------------
plot.Modelled.parameters <- function(result,with.storage) { 
  # cbPalette = c("gray", "skyblue", "orange", "green3", "yellow3", "#0072B2", "#D55E00")
  cbPalette = c("cyan", "darkorange", "firebrick", "deepskyblue3")
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
  no.param.par.var = result[[1]]
  summary.param = result[[2]]
  summary.data = result[[3]]
  summary.output = result[[4]]
  summary.error = result[[5]]
  
  for (p in 1:length(var)) {
    summary.param.set.limit = subset(summary.param, variable %in% var[p])
    for (z in 1:length(no.param.par.var)) {
      summary.param.set = subset(summary.param, variable %in% var[p] & no.param %in% no.param.par.var[z])
      i = i + 1
      plot[[i]] = ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = treatment, colour=factor(treatment))) +
        # geom_ribbon(data = summary.param.set, aes(ymin=Parameter-Parameter_SD, ymax=Parameter+Parameter_SD), linetype=2, alpha=0.1,size=0.1) +
        geom_point(position=pd,size=0.01) +
        geom_line(position=pd,data = summary.param.set, aes(x = Date, y = Parameter,  group = treatment, colour=factor(treatment)),size=1) +
        # ylab(paste(as.character(var[p]),"(fraction)")) +
        ylab(paste(as.character(var[p]))) +
        labs(colour="Treatment") +
        scale_color_manual(values=cbPalette[1:4]) +
        scale_y_continuous(limits = c(min(summary.param.set.limit$Parameter), max(summary.param.set.limit$Parameter))) +
        annotate("text", x = min(summary.param.set$Date), y = max(summary.param.set$Parameter), size = font.size-7, label = paste(title[p])) +
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
plot.Modelled.biomass <- function(result,with.storage) { 
  # cbPalette = c("gray", "skyblue", "orange", "green3", "yellow3", "#0072B2", "#D55E00")
  cbPalette = c("cyan", "darkorange", "firebrick", "deepskyblue3")
  i = 0
  font.size = 10
  plot = list() 
  no.param.par.var = result[[1]]
  summary.param = result[[2]]
  summary.data = result[[3]]
  summary.output = result[[4]]
  summary.error = result[[5]]
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
    summary.data.Cpool = subset(summary.data,variable %in% meas[p])
    summary.output.Cpool = subset(summary.output,variable %in% res[p])
    summary.error.Cpool = subset(summary.error,variable %in% error[p])
    # if (p==4) {
    #   summary.output.Cpool$value = cumsum(summary.output.Cpool$value)
    # }
    
    i = i + 1
    if (meas[p]=="Ra") {
      plot[[i]] = ggplot(summary.error.Cpool, aes(x=Date, y=parameter, group = treatment, colour=treatment)) + 
        geom_point(position=pd,size=0.3) +
        # geom_errorbar(position=pd,aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=0.5) +
        # geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,volume.group,no.param), linetype=volume.group, colour=volume, size=no.param)) +
        # geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,no.param), linetype=no.param, colour=volume)) +
        geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = treatment, colour=treatment)) +
        ylab(paste(as.character(meas[p]),"(g C)")) + xlab("Month") +
        # ggtitle("C pools - Measured (points) vs Modelled (lines)") +
        # labs(colour="Soil Volume", linetype="Grouping treatment", size="Total No of Parameter") +
        # labs(colour="Pot Volume (L)", linetype="No. of Parameters") +
        labs(colour="Treatment") +
        scale_color_manual(values=cbPalette[1:4]) +
        # scale_color_manual(labels = c("Individuals", "One Group"), values = c("blue", "red")) +
        # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g C plant-1)")) +
        theme_bw() +
        annotate("text", x = max(summary.output.Cpool$Date), y = min(summary.output.Cpool$value), size = font.size-7, label = paste(title[p])) +
        theme(legend.title = element_text(colour="black", size=font.size-2)) +
        theme(legend.text = element_text(colour="black", size = font.size-3)) +
        theme(legend.key.height=unit(0.6,"line")) +
        theme(legend.position = c(0.22,0.8)) +
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
        geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = treatment, colour=treatment)) +
        ylab(paste(as.character(meas[p]),"(g C)")) + xlab("Month") +
        labs(colour="Treatment") +
        scale_color_manual(values=cbPalette[1:4]) +
        # scale_color_manual(labels = c("Individuals", "One Group"), values = c("blue", "red")) +
        # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g C plant-1)")) +
        theme_bw() +
        annotate("text", x = max(summary.output.Cpool$Date), y = min(summary.output.Cpool$value), size = font.size-7, label = paste(title[p])) +
        # theme(plot.title = element_text(size = 20, face = "bold")) +
        theme(legend.title = element_text(colour="black", size=font.size-2)) +
        theme(legend.text = element_text(colour="black", size = font.size-3)) +
        theme(legend.key.height=unit(0.6,"line")) +
        theme(legend.position = c(0.22,0.8)) +
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
    if (p>1) {
      plot[[i]] = plot[[i]] + guides(colour=FALSE)
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

#-------------------------------------------------------------------------------------
# Read data
treat.group = as.factor(c("ambient","elevated")) # Assign all treatments
data.all = read.csv("processed_data/data_all.csv") 
tnc.partitioning = read.csv("processed_data/tnc_partitioning_data.csv")
# data.all[,c("TNC_tot_SE","TNC_leaf_SE","TNC_wood_SE","TNC_root_SE")] = 3 * data.all[,c("TNC_tot_SE","TNC_leaf_SE","TNC_wood_SE","TNC_root_SE")]

# inputs
no.param.par.var = 9
with.storage = T
model.comparison = F
model.optimization = F

param.vary = ceiling(nrow(data.all)/2/no.param.par.var)
no.param = ceiling(nrow(data.all)/2/param.vary)

# Setting parameters values
param.1.k <- matrix(c(0.2,0.2,0.2,0.2,0.25,0.1,0.08,0.1,0.15) , nrow=no.param, ncol=1, byrow=T) 
param.1.Y <- matrix(c(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3) , nrow=no.param, ncol=1, byrow=T) 
param.1.af <- matrix(c(0.3,0.25,0.25,0.2,0.2,0.2,0.2,0.15,0.05) , nrow=no.param, ncol=1, byrow=T) 
param.1.as <- matrix(c(0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.7,0.75) , nrow=no.param, ncol=1, byrow=T) 
param.1.sf <- matrix(c(0.00075,0.00075,0.00075,0.00075,0.0005,0.0005,0.00025,0.00025,0.00025) , nrow=no.param, ncol=1, byrow=T) 
param.1.sr <- matrix(c(0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001) , nrow=no.param, ncol=1, byrow=T)

param.1 = data.frame(param.1.k,param.1.Y,param.1.af,param.1.as,param.1.sf,param.1.sr)
names(param.1) <- c("k","Y","af","as","sf","sr")

param.2.k <- matrix(c(0.2,0.2,0.15,0.1,0.2,0.05,0.05,0.1,0.15) , nrow=no.param, ncol=1, byrow=T) 
param.2.Y <- matrix(c(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3) , nrow=no.param, ncol=1, byrow=T) 
param.2.af <- matrix(c(0.3,0.3,0.3,0.25,0.25,0.15,0.15,0.1,0.025) , nrow=no.param, ncol=1, byrow=T) 
param.2.as <- matrix(c(0.6,0.6,0.65,0.7,0.7,0.75,0.75,0.7,0.7) , nrow=no.param, ncol=1, byrow=T) 
param.2.sf <- matrix(c(0.00075,0.0005,0.0005,0.0005,0.00075,0.00075,0.00095,0.00095,0.00095) , nrow=no.param, ncol=1, byrow=T) 
param.2.sr <- matrix(c(0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001) , nrow=no.param, ncol=1, byrow=T)

param.2 = data.frame(param.2.k,param.2.Y,param.2.af,param.2.as,param.2.sf,param.2.sr)
names(param.2) <- c("k","Y","af","as","sf","sr")

param.set = list(param.1, param.2)

result = CBM.wtc3.manual(no.param.par.var, treat.group, with.storage, model.comparison, model.optimization, param.set) # Quadratic/Cubic parameters

result[[6]]
write.csv(result[[6]], "output/bic.csv", row.names=FALSE) # unit of respiration rates: gC per gC plant per day

# Plot parameters and biomass data fit
plot.Modelled.parameters(result,with.storage)
plot.Modelled.biomass(result,with.storage)


