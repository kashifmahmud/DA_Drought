#----------------------------------------------------------------------------------------------------------------
#- This is a collection of many functions that do the actual work of data analysis and plotting. These functions
#    are called by just a few lines of code in "main script.R" to recreate the analyses and figures.
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
# Demonstration of the short-term temperature effects on photosynthesis and respiration.

plotCUE_conceptual_fig <- function(toexport=T,Tdew=10,Ca=400,Vcmax=100,Jmax=125,Tleaf=20:42,PPFD=1500){
  # load required libraries
  #- set up parameters for model
  #VPD <- RHtoVPD(RH=RH,TdegC=Tleaf)
  VPD <- DewtoVPD(Tdew=Tdew,TdegC=Tleaf)
  
  #- predict photosynthesis and respiration with changing T (and VPD)
  output<- Photosyn(VPD=VPD,Ca=Ca,Vcmax=Vcmax,Jmax=Jmax,Tleaf=Tleaf,Tcorrect=T,Rd0=1,TrefR=20,PPFD=PPFD)
  output$AGROSS <- with(output,ALEAF+Rd)
  output.acclim <-Photosyn(VPD=VPD,Ca=Ca,Vcmax=Vcmax*0.95,Jmax=Jmax*0.95,Tleaf=Tleaf,Tcorrect=T,Rd0=0.65,TrefR=20,delsJ=630,PPFD=PPFD)
  
  output.acclim$AGROSS <- with(output.acclim,ALEAF+Rd)
  
  #- calculate CUE
  CUE <- 1-(output$Rd/output$ALEAF)
  CUEa <- 1-(output.acclim$Rd/output.acclim$ALEAF)
  RtoA <- (output$Rd/output$AGROSS)
  RtoA2 <- (output.acclim$Rd/output.acclim$AGROSS)
  
  #- plot
  # dev.new(width=20,height=40)
  # windows(20,40);
  par(mfrow=c(2,1),mar=c(2,7,1,1),oma=c(3,0,2,0),las=1)
  plot(AGROSS~Tleaf,data=output,type="l",ylab="",cex.lab=1.4,ylim=c(0,25),
       col="forestgreen",lwd=2,xaxt="n",yaxt="n",frame.plot=F)
  magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,majorn=3)
  title(ylab=expression(atop(Leaf~CO[2]~exchange,
                             (mu*mol~CO[2]~m^-2~s^-1))),cex.lab=1.3)
  lines(AGROSS~Tleaf,data=output.acclim,type="l",col="forestgreen",lwd=2,lty=2)
  lines(Rd~Tleaf,data=output,col="brown",lwd=2)
  lines(Rd~Tleaf,data=output.acclim,col="brown",lwd=2,lty=2)
  legend(x=23,y=30,xpd=NA,legend=c("A","R"),lwd=2,col=c("forestgreen","brown"),ncol=2,bty="n")
  legend("topleft","a",cex=1.1,bty="n",inset=-0.05)
  
  plot(RtoA~output$Tleaf,type="l",ylab="",ylim=c(0,0.6),cex.lab=1.3,lwd=2,col="black",xaxt="n",yaxt="n",
       frame.plot=F)
  magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,majorn=3)
  title(ylab=expression(R/A),cex.lab=1.3)
  
  lines(RtoA2~output.acclim$Tleaf,type="l",ylab="CUE (1-R/A)",ylim=c(0.2,1),cex.lab=1.3,lwd=2,lty=2,col="black")
  legend("topleft","b",cex=1.1,bty="n",inset=-0.05)
  
  par(xpd=NA)
  text(x=32,y=-0.15,expression(T[leaf]~(degree*C)),cex=1.4)
  
  #- overlay conceptual points
  
#   Arrows(x0=c(35),y0=c(0.14),x1=c(39),y1=c(0.29),lwd=3,col=c("red"),arr.type="curved")
#   points(x=c(35,40),y=c(0.14,0.32),pch=16,col=c("black","red"),cex=2.5)
#   
  
  Arrows(x0=c(35,35),y0=c(0.14,0.14),x1=c(39,38),y1=c(0.29,0.14),lwd=3,col=c("red","orange"),arr.type="curved")
  points(x=c(35,40,40),y=c(0.14,0.32,0.14),pch=16,col=c("black","red","orange"),cex=2.5)
  
  
  if(toexport==T) dev.copy2pdf(file="./output/Figure1.pdf")
}
#----------------------------------------------------------------------------------------------------------------






#-----------------------------------------------------------------------------------------------------
#- function to return leak parameters of WTC fluxes while in a closed state
return.leaks <- function(plotson=0){

  # I did leak tests for WTC chambers #1-4 on the night of 27 May 2014, and all 12 chambers on 28 May 2014.

  #read in all the minutely data
  leakdat <- read.csv(file="data/WTC_TEMP_CM_WTCFLUX-LEAKS_20140527-20140529_L1.csv")
  dat <- subset(leakdat,chamber<13)
  dat$datetime <- as.POSIXct(paste(dat$date,dat$time,sep=" "),format="%Y-%m-%d %H:%M:%S",tz="GMT")
  
  #set up the times for start and end times for the two nights
  starttime1 <- as.POSIXct("2014-05-27 22:00:10",format="%Y-%m-%d %H:%M:%S",tz="GMT")
  endtime1 <- as.POSIXct("2014-05-28 07:30:00",format="%Y-%m-%d %H:%M:%S",tz="GMT")
  starttime2 <- as.POSIXct("2014-05-28 20:50:00",format="%Y-%m-%d %H:%M:%S",tz="GMT")
  endtime2 <- as.POSIXct("2014-05-29 07:30:10",format="%Y-%m-%d %H:%M:%S",tz="GMT")
  
  
  print("Processing data...")

    #------------------------------------------------------------------------------------------------------------------
  #work out the reference [CO2] data. Note that this is chamber "13", and the data has a different structure.
  ref <- subset(leakdat,chamber==13)
  ref$datetime <- as.POSIXct(paste(ref$date,ref$time,sep=" "),format="%Y-%m-%d %H:%M:%S",tz="GMT")
  
  # just grab the [CO2] data, interpolate it
  ref2 <- ref[,c(13,19)]
  names(ref2)[1] <- "CO2ref"
  times <- data.frame(datetime =seq.POSIXt(from=starttime1,to=endtime2,by="sec"))
  
  ref.i <- merge(times,ref2,all.x=T)
  ref.i <- ref.i[with(ref.i,order(datetime)),]
  
  ref.zoo <- zoo(ref.i$CO2ref)
  ref.zoo <- na.approx(ref.zoo) #gapfill between observaitons
  
  ref.i$CO2ref <- ref.zoo #replace missing values with gapfilled numbers
  #------------------------------------------------------------------------------------------------------------------
  
  
  
  #------------------------------------------------------------------------------------------------------------------
  #subset the data for the first night
  dat1 <- subset(dat,datetime>starttime1&datetime<endtime1)
  dat1 <- suppressWarnings(merge(ref.i,dat1,by="datetime"))
  dat1 <- suppressWarnings(subset(dat1,chamber<=4)) #only the first four chambers were leak tested on the first day
  
  #subset the data for the second night
  dat2 <- subset(dat,datetime>starttime2&datetime<endtime2)
  dat2 <- suppressWarnings(merge(ref.i,dat2,by="datetime"))
  
  
  
  #model the leak, return the squared residual sum for minimiation by optimise(), which is much faster than DEoptim().
  leak.mod <- function(theta, V, Ca, Ci,fit=1){
    
    resid <- rep(NA,length(Ca))
    resid[1] <- 0
    pred <- rep(NA,length(Ca))
    pred[1] <- Ci[1]
    
    for (i in 2:length(Ca)){
      dCO2 <- -theta*(Ci[i-1]-Ca[i-1])/V
      pred[i] <- pred[i-1] + dCO2
      
      resid[i] <- (Ci[i] - pred[i])^2
    }
    resid.sum <- sum(resid)
    if (fit==1) return(resid.sum)
    if (fit==0) return(data.frame(Ci=Ci,Ca=Ca,pred=pred))
  }
  
  
  print("Fitting the model to estimate leak parameters for the first day of measurements...")
  #fit the first day
  dat1.list <- split(dat1,dat1$chamber)
  out1 <- list()
  pred <- list()
  pred[[1]] <- data.frame()
  for (i in 1:length(dat1.list)){
    dat <- dat1.list[[i]]
    
    out1[i] <- suppressWarnings(optimise(f=leak.mod,interval=c(0,0.5),V=53,Ca=dat$CO2ref,Ci=dat$CO2L,fit=1))
    
  }
  
  print("Fitting the model to estimate leak parameters for the second day of measurements...")
  
  #fit the second day
  dat2.list <- split(dat2,dat2$chamber)
  out2 <- list()
  # pred <- list()
  # pred[[1]] <- data.frame()
  for (i in 1:length(dat2.list)){
    dat <- dat2.list[[i]]
    
    out2[i] <- suppressWarnings(optimise(f=leak.mod,interval=c(0,0.5),V=53,Ca=dat$CO2ref,Ci=dat$CO2L,fit=1))
  }
  
  leaks_day1 <- data.frame(theta=do.call(rbind,out1))
  leaks_day1$chamber <- 1:4
  
  
  #compile leak estimates
  leaks <- data.frame(theta=do.call(rbind,out2))
  leaks$chamber <- 1:12
  #leaks$theta2 <- c(do.call(rbind,out1),rep(NA,8))
  leaks <- leaks[,c(2,1)]
  #leaks$theta <- rowMeans(leaks[,2:3],na.rm=T)
  
  
  #-- what was the variation of theta measured across the two nights?
  dtheta <- rbind(leaks_day1,leaks[1:4,])
  leaksd <- summaryBy(theta~chamber,data=dtheta,FUN=sd)
  leakmean <- summaryBy(theta~chamber,data=dtheta,FUN=mean)
  cv <- merge(leaksd,leakmean,by="chamber")
  cv$cv <- with(cv,theta.sd/theta.mean)
  #----------------------------------------------------------------------------------------------
  #plot fits, if plotson==1
  
  if(plotson==1){
    print("Plotting")
    # windows(12,10)
    par(mfrow=c(4,3),mar=c(2,4,2,2))
    
    #fit the second day
    dat2.list <- split(dat2,dat2$chamber)
    out2 <- list()
    pred <- data.frame(Ci=NA,Ca=NA,pred=NA)
    # pred[[1]] <- data.frame()
    for (i in 1:length(dat2.list)){
      dat <- dat2.list[[i]]
      
      #out2[i] <- optimise(f=leak.mod,interval=c(0,0.5),V=30,Ca=dat$CO2ref,Ci=dat$CO2L,fit=1)
      pred <- leak.mod(theta=leaks$theta[i],V=30,Ca=dat$CO2ref,Ci=dat$CO2L,fit=0)
      
      plot(pred$Ci,col="black",lty=1,lwd=2,type="l",ylab="[CO2]",ylim=c(350,1100))
      lines(pred$Ca,col="grey",lty=1,lwd=2)
      lines(pred$pred,col="red",lty=1,lwd=2)
      title(main=paste("Chamber",i,sep=" "))
      if (i==3){
        legend("topright",c("data","pred","Ca"),lty=1,lwd=2,col=c("black","red","grey"),cex=1.5,bty="n")
      }
    }
  }
  
  
  print("Done with leak calculations.")
  return(leaks)
}
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#- code to return WTC flux estimates from closed chambers (whole canopy T-response curves)

return_Rcanopy_closed <- function(){

  # get the leak parameter for each chamber. This will be used later to estimate the rate of respiration
  leaks <- suppressWarnings(return.leaks(plotson=0))
  
  #------------------------------------------------------
  #------------------------------------------------------
  #------------------------------------------------------
  #read in the minutely IRGA data, including a local IRGA (CO2L) and a central LI7000 (CO2CConc)
  Mindat1 <- read.csv("data/WTC_TEMP_CM_WTCFLUX-CLOSED_20140212-20140213_L1.csv")
  Mindat1$DateTime <- paste(Mindat1$date,Mindat1$time,sep=" ")
  Mindat1$DateTime <- as.POSIXct(Mindat1$DateTime,tz="GMT")
  
  
  
  #5 sets of measurements
  #start and stop times are different for this computer relative to the tdl computer
  starttime1 <- as.POSIXct("2014-02-12 19:16:09 GMT",format="%Y-%m-%d %T",tz="GMT")
  stoptime1 <- as.POSIXct("2014-02-12 20:45:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  starttime2 <- as.POSIXct("2014-02-12 21:05:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  stoptime2 <- as.POSIXct("2014-02-12 22:40:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  starttime3 <- as.POSIXct("2014-02-12 23:00:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  stoptime3 <- as.POSIXct("2014-02-13 00:30:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  starttime4 <- as.POSIXct("2014-02-13 01:00:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  stoptime4 <- as.POSIXct("2014-02-13 02:30:00 GMT",format="%Y-%m-%d %T",tz="GMT")
  starttime5 <- as.POSIXct("2014-02-13 03:40:00 GMT",format="%Y-%m-%d %T",tz="GMT") 
  stoptime5 <- as.POSIXct("2014-02-13 04:29:09 GMT",format="%Y-%m-%d %T",tz="GMT")
  
  #Mindat2 <- subset(Mindat,DateTime>starttime1 & DateTime<= stoptime5)
  
  
  #calcualte time elapsed
  dat1.min <- subset(Mindat1,DateTime >= starttime1 & DateTime <= stoptime1);dat1.min$Measurement <- 1
  dat2.min <- subset(Mindat1,DateTime >= starttime2 & DateTime <= stoptime2);dat2.min$Measurement <- 2
  dat3.min <- subset(Mindat1,DateTime >= starttime3 & DateTime <= stoptime3);dat3.min$Measurement <- 3
  dat4.min <- subset(Mindat1,DateTime >= starttime4 & DateTime <= stoptime4);dat4.min$Measurement <- 4
  dat5.min <- subset(Mindat1,DateTime >= starttime5 & DateTime <= stoptime5);dat5.min$Measurement <- 5
  
  Mindat <- rbind(dat1.min,dat2.min,dat3.min,dat4.min,dat5.min)
  Mindat$T_treatment <- ifelse(Mindat$chamber %% 2 ==1,"ambient","elevated")
  
  #------------------------------------------------------------------------------------------------------------------
  #work out the reference [CO2] data. Note that this is chamber "13", and the data has a different structure.
  #note that the actual measured [CO2] by the central irga of the 
  ref <- subset(Mindat,chamber==13)
  
  # just grab the [CO2] data, interpolate it. Note that teh different structure of the reference data
  #  meant that the CO2 data were actually stored in the column labeled "DPCChamb", but ONLY for chamber =13
  ref2 <- ref[,c("DPCChamb","DateTime")]
  names(ref2)[1] <- "CO2ref"
  times <- data.frame(DateTime =seq.POSIXt(from=starttime1,to=stoptime5,by="sec",tz="GMT"))
  
  ref.i <- merge(times,ref2,all=T)
  ref.i <- ref.i[with(ref.i,order(DateTime)),]
  
  ref.zoo <- zoo(ref.i)
  ref.zoo$CO2ref <- na.approx(ref.zoo$CO2ref) #gapfill between observaitons
  ref.df <- fortify.zoo(ref.zoo)
  ref.df$CO2ref <- as.numeric(as.character(ref.df$CO2ref))
  ref.df$DateTime <- as.POSIXct(ref.df$DateTime,tz="GMT")
  

  #------------------------------------------------------------------------------------------------------------------
  
  
  
  #------------------------------------------------------------------------------------------------------------------
  #merge the chamber data with the interpolated reference data
  dat1 <- merge(ref.df,subset(Mindat,chamber<13),by="DateTime")
  
  #-- so that worked great, except for chamber 9, where the local irga died. Ugh.
  # I"ll just fit the first three temperatures for chamber 9, which have good data
  dat9 <- subset(dat1,chamber==9)
  dat92 <- subset(dat9,DateTime <as.POSIXct("2014-02-13 00:35:00",tz="GMT"))
  
  dat1.1 <- subset(dat1,chamber !=9)
  
  dat2 <- rbind(dat1.1,dat92)
  #------------------------------------------------------------------------------------------------------------------
  
  
  
  #----------------------------------------------------------------------------------------------
  # get the chamber temperatures from the inside met data from HIEv
  met <- read.csv("data/WTC_TEMP_CM_WTCMET_20140201-20140228_L1_v1.csv")
  met$DateTime <- as.POSIXct(met$DateTime,tz="GMT")
  #met2 <- subset(met,DateTime>starttime1 & DateTime < stoptime5)
  
  #set up the data for merging with the met
  dat2$chamber <- as.factor(paste0("C",sprintf("%02.0f",dat2$chamber)))
  dat2$DateTime1 <- nearestTimeStep(dat2$DateTime,nminutes=15,align="ceiling")
  
  dat3 <- merge(dat2,met,by.x=c("chamber","DateTime1"),by.y=c("chamber","DateTime"))
  #----------------------------------------------------------------------------------------------
  
  
  #------------------------------------------------------------------------------------------------------------------
  # create a factor for measurement number
  dat3$id <- as.factor(paste(dat3$chamber,dat3$Measurement,sep="-"))
  dat3.l <- split(dat3,dat3$id)
  #------------------------------------------------------------------------------------------------------------------
  
  
  
  
  
  #----------------------------------------------------------------------------------------------
  #model the flux, given the leak rate, return the squared residual sum for minimiation by optimise(), which is much faster than DEoptim().
  closedWTC.mod <- function(R,theta, V, Ca, Ci,fit=1){
    
    resid <- rep(NA,length(Ca))
    resid[1] <- 0
    pred <- rep(NA,length(Ca))
    pred[1] <- Ci[1]
    
    for (i in 2:length(Ca)){
      dCO2 <- (-theta*(Ci[i-1]-Ca[i-1])+R)/V
      pred[i] <- pred[i-1] + dCO2
      
      resid[i] <- (Ci[i] - pred[i])^2
    }
    resid.sum <- sum(resid)
    if (fit==1) return(resid.sum)
    if (fit==0) return(data.frame(Ci=Ci,Ca=Ca,pred=pred))
  }
  #----------------------------------------------------------------------------------------------
  
  
  
  #----------------------------------------------------------------------------------------------
  out1 <- c()
  fits1 <- c()
  #fit the respiration model using the leak estimate measured with an empty chamber
  for (i in 1:length(dat3.l)){
    dat <- dat3.l[[i]]
    theta <- leaks[dat$chamber[1],2]
    out1[[i]] <- optimise(f=closedWTC.mod,interval=c(0,200),theta=theta,V=53,Ca=dat$CO2ref,Ci=dat$CO2L,fit=1)[[1]] # V was 30
    chamber <- dat$chamber[1]
    id <- dat$id[1]
    T_treatment <- dat$T_treatment[1]
    Tair <- median(dat$Tair_al)
    RH <- median(dat$RH_al)
    Measurement <- dat$Measurement[1]
    
    fits1[[i]] <- data.frame(chamber,id,T_treatment,Measurement,theta,Tair,RH,Rcanopy=out1[[i]])
  }
  fits <- do.call(rbind,fits1)
  fits$mol_m3 <- 44.6 
  
  
  
  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  #- Okay, now we have raw respiration data. Let's merge it with tree size to normalize the flux by mass (or leaf area)
  
  #- get an estimate of leaf area for each day of the experiment
  #treeMass <- returnTreeMass()
  treeMass <- read.csv("data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
  treeMass$DateTime <- as.POSIXct(treeMass$DateTime,format="%Y-%m-%d %T",tz="GMT")
  treeMass2 <- subset(treeMass,DateTime==as.POSIXct("2014-02-12 01:00:00",format="%Y-%m-%d %T",tz="GMT"),
                      select=c("chamber","leafArea"))
  
  
  #- merge the respiration data with the tree mass data
  fits.mass <- merge(fits,treeMass2,by=c("chamber"))
  
  #- leaf area of ambient and warmed treatments
  summaryBy(leafArea~T_treatment,data=fits.mass,FUN=mean)
  
  #- unit conversions. Recall that Rcanopy has units of umol CO2 m^3 mol-1 s-1
  fits.mass$Rcanopy_umol <- fits.mass$Rcanopy*fits.mass$mol_m3/60
  
  fits.mass$R_la <- with(fits.mass,Rcanopy_umol/leafArea) # convert from umol CO2 s-1 to umol CO2 m-2 s-1
  fits.mass$mol_m3 <- NULL
  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  
  
  #average by treatment 
  fits.mass$Tbin <- ifelse(fits.mass$Tair<20,18,
                           ifelse(fits.mass$Tair>20 & fits.mass$Tair < 24,22.5,
                                  ifelse(fits.mass$Tair>24 & fits.mass$Tair < 26, 25,
                                         ifelse(fits.mass$Tair>26 & fits.mass$Tair<30,28,32))))
  fits.trt <- summaryBy(Tair+Rcanopy_umol+R_la~T_treatment+Tbin,data=fits.mass,FUN=c(mean,standard.error),keep.names=F)
  
  #----------------------------------------------------------------------------------------------
  
  return(list(fits.mass,fits.trt))
  
}
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#- plot figure 2 (R vs. T curves at leaf and whole canopy scales)
plotRvsT_figure2 <- function(fits.mass=fits.mass,fits.trt=fits.trt,export=T){
  palette(c("black","red"))
  
  #--------------------------------------------------------------------------------------------------------
  #- fit arrheneous models for all three types of respiraiton
  # ANCOVA. No difference in slopes, so interaction removed from model
  fits.mass$logRcanopy <- log(fits.mass$Rcanopy_umol)
  fits.mass$logRarea <- log(fits.mass$R_la)
  
  
  fits.mass$invT <- 1/(fits.mass$Tair+273.15)
  
  lmRcanopy <- lm(logRcanopy~invT+T_treatment,data=fits.mass)
  lmRarea <- lm(logRarea~invT+T_treatment,data=fits.mass)
  #- extract fits for each type of R
  Rgas <- 8.314 #J per mol per K
  
  A_amb1 <- unname(exp(coef(lmRcanopy)[1]))
  A_ele1 <- unname(exp(coef(lmRcanopy)[1]+coef(lmRcanopy)[3]))
  Ea1 <- unname(-1*coef(lmRcanopy)[2]*Rgas)
  #(A_ele1-A_amb1)/A_amb1*100
  
  A_amb2 <- unname(exp(coef(lmRarea)[1]))
  A_ele2 <- unname(exp(coef(lmRarea)[1]+coef(lmRarea)[3]))
  Ea2 <- unname(-1*coef(lmRarea)[2]*Rgas)
  #(A_ele2-A_amb2)/A_amb2*100
  
  
  Q10 <- exp(10*Ea1/(Rgas*(25+273.15)^2))
  # get model predictions of R across all temperatures
  xvals <- seq(18,32.1, length=101)
  predA1 <- A_amb1*exp(-1*Ea1/(Rgas*(xvals+273.15)))
  predE1 <- A_ele1*exp(-1*Ea1/(Rgas*(xvals+273.15)))
  predA2 <- A_amb2*exp(-1*Ea2/(Rgas*(xvals+273.15)))
  predE2 <- A_ele2*exp(-1*Ea2/(Rgas*(xvals+273.15)))
  #--------------------------------------------------------------------------------------------------------
  
  
  
  
  #--------------------------------------------------------------------------------------------------------
  #- process leaf-scale R vs. T data
  rvt <- read.csv("data/WTC_TEMP_CM_GX-RdarkVsT_20140207-20140423_L1.csv")
  rvt$Date <- as.Date(rvt$Date)
  rvt.c <- subset(rvt,Date==as.Date("2014-02-07")) #- just pull out the data measured prior to the drought
  
  rvt.45 <- subset(rvt.c, Tleaf > 18 & Tleaf<=40)
  rvt.45$lnRmass <- log(rvt.45$Rmass)
  rvt.45$lnRarea <- log(rvt.45$Rarea)
  rvt.45$invT <- 1/(rvt.45$Tleaf+273.15)
  rvt.45$Treat <- as.factor(rvt.45$T_treatment)
  rvt.45$datefac <- as.factor(rvt.45$Date)
  
  rvt.45$Tleaf_bin <- cut(rvt.45$Tleaf,breaks=seq(from=18,to=40,length=25))
  rvt.45$Tleaf_bin_mid <- sapply(strsplit(gsub("^\\W|\\W$", "", rvt.45$Tleaf_bin), ","), function(x)sum(as.numeric(x))/2) 
  rvt.45.ch <- summaryBy(Rmass+invT+lnRmass+Rarea+lnRarea~date+chamber+Treat+Tleaf_bin_mid,data=subset(rvt.45, Tleaf_bin_mid<45),
                         keep.names=T,FUN=c(mean))
  
  rvt.treat.bin <- summaryBy(Rmass+Rarea~date+Treat+Tleaf_bin_mid,data=rvt.45.ch,keep.names=F,FUN=c(mean,standard.error))
  rvt.treat.bin$Rmass.high <- with(rvt.treat.bin,Rmass.mean+Rmass.standard.error)
  rvt.treat.bin$Rmass.low <- with(rvt.treat.bin,Rmass.mean-Rmass.standard.error)
  
  #- fit arrhenious. 
  lmRleaf <- lm(lnRarea~invT+Treat,data=rvt.45.ch)
  Rleaf_amb1 <- unname(exp(coef(lmRleaf)[1]))
  Rleaf_ele1 <- unname(exp(coef(lmRleaf)[1]+coef(lmRleaf)[3]))
  Ea_Rleaf <- unname(-1*coef(lmRleaf)[2]*Rgas)

  
  Q10_Rleaf <- exp(10*Ea_Rleaf/(Rgas*(25+273.15)^2))
  # get model predictions of R across all temperatures
  xvals_Rleaf <- seq(18,40, length=101)
  predRleafA <- Rleaf_amb1*exp(-1*Ea_Rleaf/(Rgas*(xvals_Rleaf+273.15)))
  predRleafE <- Rleaf_ele1*exp(-1*Ea_Rleaf/(Rgas*(xvals_Rleaf+273.15)))
  
  #--------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------
  
  
  
  
  
  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  #- plot R vs. T with (a) raw fluxes, (b) fluxes per unit leaf area, (c) fluxes per unit total tree mass
  # windows(22,32);
  par(mfrow=c(3,1),cex.lab=1.5,mar=c(2,7,1,2),oma=c(4,1,0,0),las=1)
  xlims=c(15,42)
  
  # #- plot AREA BASED leaf R over T for the first date of high resolution T-response curves
  plotBy(Rarea.mean~Tleaf_bin_mid|Treat,data=rvt.treat.bin,xaxt="n",yaxt="n",
         ylab=expression(atop(R[leaf],
                              (mu*mol~CO[2]~m^-2~s^-1))),col=c("black","red"),pch=15,
         xlim=xlims,ylim=c(0,4),type="p",lwd=3,cex=1.6,xlab="",cex.lab=1.6,legend=F,
         panel.first=adderrorbars(x=rvt.treat.bin$Tleaf_bin_mid,y=rvt.treat.bin$Rarea.mean,SE=rvt.treat.bin$Rarea.standard.error,direction="updown"))
  magaxis(side=c(1,2,3,4),labels=c(1,1,0,1))
  legend("bottomright","a",bty="n",inset=0.02,cex=1.5)
  lines(x=xvals_Rleaf,y=predRleafA,col="black",lwd=2)
  lines(x=xvals_Rleaf,y=predRleafE,col="red",lwd=2)
  
    #- plot raw T-response curves
  plotBy(Rcanopy_umol.mean~Tair.mean|T_treatment,data=fits.trt,type="p",xlim=xlims,ylim=c(0,30),pch=15,cex=1.6,cex.lab=1.6,xaxt="n",yaxt="n",
         ylab=expression(atop(R[canopy],
                              (mu*mol~CO[2]~s^-1))),xlab="",legend=F,
         panel.first=adderrorbars(x=fits.trt$Tair.mean,y=fits.trt$Rcanopy_umol.mean,SE=fits.trt$Rcanopy_umol.standard.error,direction="updown"))
  lines(x=xvals,y=predA1,col="black",lwd=2)
  lines(x=xvals,y=predE1,col="red",lwd=2)
  magaxis(side=c(1,2,3,4),labels=c(1,1,0,1))
  legend("bottomright","b",bty="n",inset=0.02,cex=1.5)
  
  #- plot respiration per unit leaf area
  plotBy(R_la.mean~Tair.mean|T_treatment,data=fits.trt,type="p",pch=15,ylim=c(0,2.5),xlim=xlims,cex=1.6,cex.lab=1.6,legend=F,xaxt="n",yaxt="n",
         ylab=expression(atop(R["canopy, area"],
                              (mu*mol~CO[2]~m^-2~s^-1))),xlab=expression(T[air]~(degree*C)),
         panel.first=adderrorbars(x=fits.trt$Tair.mean,y=fits.trt$R_la.mean,SE=fits.trt$R_la.standard.error,direction="updown"))
  lines(x=xvals,y=predA2,col="black",lwd=2)
  lines(x=xvals,y=predE2,col="red",lwd=2)
  legend("bottomright","c",bty="n",inset=0.02,cex=1.5)
  magaxis(side=c(1,2,3,4),labels=c(1,1,0,1))
  title(xlab=expression(Temperature~(degree*C)),xpd=NA)
  
  if(export==T) dev.copy2pdf(file="output/Figure2.pdf")
}
#----------------------------------------------------------------------------------------------------------------






#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#- Compare linear and exponential fits of the canopy respiration dataset
compareRlinearRexp <- function(fits.trt){
  xlims <- c(15,35)
  
  #- fit exponential
  fits.trt$invT <- 1/(fits.trt$Tair.mean+273.15)
  fits.trt$logRarea <- log(fits.trt$R_la.mean)
  lmRarea <- lm(logRarea~invT+T_treatment,data=fits.trt)
  Rgas <- 8.314 #J per mol per K
  A_amb2 <- unname(exp(coef(lmRarea)[1]))
  A_ele2 <- unname(exp(coef(lmRarea)[1]+coef(lmRarea)[3]))
  Ea2 <- unname(-1*coef(lmRarea)[2]*Rgas)
  Q10 <- exp(10*Ea2/(Rgas*(25+273.15)^2))
  xvals <- seq(18,32.1, length=101)
  predA2 <- A_amb2*exp(-1*Ea2/(Rgas*(xvals+273.15)))
  predE2 <- A_ele2*exp(-1*Ea2/(Rgas*(xvals+273.15)))
  
  #- fit linear
  lmRlinear <- lm(R_la.mean~Tair.mean*T_treatment,data=fits.trt)
  newdat <- data.frame(expand.grid(Tair.mean=xvals,T_treatment=levels(fits.trt$T_treatment)))
  newdat$pred <- predict(lmRlinear,newdat)
  
  # windows(40,30)
  par(mar=c(5,8,1,3))
  plotBy(R_la.mean~Tair.mean|T_treatment,data=fits.trt,type="p",pch=15,ylim=c(0,2),xlim=xlims,cex=1.6,cex.lab=1.6,legend=F,xaxt="n",yaxt="n",
         ylab=expression(atop(R["canopy, area"],
                              (mu*mol~CO[2]~m^-2~s^-1))),xlab=expression(T[air]~(degree*C)),
         panel.first=adderrorbars(x=fits.trt$Tair.mean,y=fits.trt$R_la.mean,SE=fits.trt$R_la.standard.error,direction="updown"))
  magaxis(side=c(1,2,3,4),labels=c(1,1,0,1))
  
  #- overlay linear predictions
  plotBy(pred~Tair.mean|T_treatment,type="l",lty=2,add=T,data=newdat,lwd=2)
  #- overlay exponential
  lines(x=xvals,y=predA2,col="black",lwd=2)
  lines(x=xvals,y=predE2,col="red",lwd=2)
  
  summary(lmRarea)
  summary(lmRlinear)
  
  

}
#----------------------------------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#- get and plot the Rleaf and Branch data near the end of the experiment
plotRleafRbranch <- function(export=T){

  Rbranch <- read.csv("data/WTC_TEMP_CM_GX-RBRANCH_20140513-20140522_L1_v1.csv")
  Rbranch$date <- as.Date(Rbranch$date)
  Rbranch$campaign <- ifelse(Rbranch$date == as.Date("2014-05-13"),1,2)
  R1 <- subset(Rbranch,campaign==1)
  
  #- average across treatments
  R1.ch <- summaryBy(Rleaf+Rbranch~T_treatment+chamber,data=R1,FUN=c(mean),keep.names=T)
  R1.m <- summaryBy(Rleaf+Rbranch~T_treatment,data=R1.ch,FUN=c(mean,standard.error))
  
  R1.m$Rleaf.cil <- R1.m$Rleaf.mean-R1.m$Rleaf.standard.error
  R1.m$Rleaf.ciu <- R1.m$Rleaf.mean+R1.m$Rleaf.standard.error
  R1.m$Rbranch.cil <- R1.m$Rbranch.mean-R1.m$Rbranch.standard.error
  R1.m$Rbranch.ciu <- R1.m$Rbranch.mean+R1.m$Rbranch.standard.error
  
  #-----------------------------------------------------------------------------------------------------------
  # plot acclimation of tissue-specific R
  # windows(8,8)
  par(mfrow=c(1,2),mar=c(5,3,3,1),oma=c(1,4,0,0),cex.lab=2.5,las=1,cex.axis=1.8,cex.lab=2)
  
  #- plot leaves
  barplot2(height=R1.m$Rleaf.mean,plot.ci=T,ci.l=R1.m$Rleaf.cil,ci.u=R1.m$Rleaf.ciu,yaxt="n",ci.width=0.2,col=c("darkgrey","red"),
           names.arg=c("A","W"),ylim=c(0,4))
  legend("topright","a",bty="n",inset=-0.002,cex=1.5)
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  title(main="Leaves",cex.main=1.5)
  title(xlab="Treatment",outer=F,line=3)
  
  #- plot branches
  barplot2(height=R1.m$Rbranch.mean,plot.ci=T,ci.l=R1.m$Rbranch.cil,ci.u=R1.m$Rbranch.ciu,yaxt="n",ci.width=0.2,col=c("darkgrey","red"),
           names.arg=c("A","W"),ylim=c(0,2))
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  legend("topright","b",bty="n",inset=-0.002,cex=1.5)
  title(main="Branch wood",cex.main=1.5)
  title(xlab="Treatment",outer=F,line=3)
  
  
  title(ylab=expression(Respiration~(nmol~g^-1~s^-1)),outer=T,line=0)
  
  #   
  #   #- statistical analyses
  #   lm.leaf <- lme(Rleaf~T_treatment,random=~1|chamber,data=R1)
  #   
  #   
  #   #look at model diagnostics
  #   plot(lm.leaf,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
  #   plot(lm.leaf,Rleaf~fitted(.),abline=c(0,1))         #predicted vs. fitted for each species
  #   qqnorm(lm.leaf, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot
  #   hist(lm.leaf$residuals)
  #   anova(lm.leaf)
  #   lsmeans(lm.leaf,"T_treatment")
  #   
  #   #- statistical analyses
  #   lm.br <- lme(Rbranch~T_treatment,random=~1|chamber,data=R1)
  #   
  #   
  #   #look at model diagnostics
  #   plot(lm.br,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
  #   plot(lm.br,Rbranch~fitted(.),abline=c(0,1))         #predicted vs. fitted for each species
  #   qqnorm(lm.br, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot
  #   hist(lm.br$residuals)
  #   anova(lm.br)
  #   lsmeans(lm.br,"T_treatment")
  #   
  
  if (export==T) dev.copy2pdf(file="output/Figure3.pdf")
  
}
#----------------------------------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------------------------
#' Adds error bars to a plot
#' 
#' @description Yet another function that adds error bars. The user must specify the length of the error bars.
#' @param x The x coordinates of start of the error bar
#' @param y The y coordinates of start of the error bar
#' @param SE The length of the error bar
#' @param direction One of 'up', 'down', 'right', 'left', 'updown' or 'rightleft'.
#' @param barlen The length of the cross-bar at the end of the error bar.
#' @param \ldots Additional parameters passed to \code{\link{arrows}}, such as the colour (\code{col}).
#' #' @details Simple wrapper for \code{\link{arrows}}, where \code{angle=90} and \code{code=3}. The \code{barlen} argument corresponds to \code{length} in \code{arrows}.
#' @examples
#' # A simple example. Also note that we can specify the colour of the error bars, or other parameters
#' # that arrows() recognizes.
#' x <- rnorm(20)
#' y <- x + rnorm(20)
#' se <- runif(20, 0.2,0.4)
#' plot(x,y,pch=21,bg="white",panel.first=adderrorbars(x,y,se,direction="updown", col="darkgrey"))
#' @export
adderrorbars <- function(x,y,SE,direction,barlen=0.04,...){
  
  if(length(direction)>1)stop("direction must be of length one.")
  #if(direction == "updown")
  #  direction <- c("up","down")
  if(direction == "rightleft" | direction == "leftright")direction <- c("left","right")
  
  if("up" %in% direction)
    arrows(x0=x, x1=x, y0=y, y1=y+SE, code=3, angle=90, length=barlen,...)
  if("down" %in% direction) 
    arrows(x0=x, x1=x, y0=y, y1=y-SE, code=3, angle=90, length=barlen,...)
  if("updown" %in% direction) 
    arrows(x0=x, x1=x, y0=y+SE, y1=y-SE, code=3, angle=90, length=barlen,...)
  if("left" %in% direction) 
    arrows(x0=x, x1=x-SE, y0=y, y1=y, code=3, angle=90, length=barlen,...)
  if("right" %in% direction)
    arrows(x0=x, x1=x+SE, y0=y, y1=y, code=3, angle=90, length=barlen,...)  
  
}
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
standard.error <- function(dat,na.rm=F,...){
  if(na.rm==T){
    dat <- subset(dat,is.na(dat)==F)
  }
  std <- sd(dat)
  n <- length(dat)
  se <- std/sqrt(n)
  return(se)
}
#----------------------------------------------------------------------------------------------------------------






#----------------------------------------------------------------------------------------------------------------
#- function to parition the hourly net flux observations into GPP and Ra using an Arrhenious function
partitionHourlyFluxCUE_arr <- function(dat.hr.gf=dat.hr.gf,Ea=57.69,lagdates,leafRtoTotal = 1,leafRreduction=0){
  rvalue = 8.134
  require(data.table)
  
  #-- convert mmolCO2 s-1 to gC hr-1
  dat.hr.gf$FluxCO2_g <- with(dat.hr.gf,FluxCO2*60*60/1000*12.0107)
  dat.hr.gf$period <- ifelse(dat.hr.gf$PAR>2,"Day","Night")
  
  #-- partition day-time net C exchange into GPP and Ra, similar to how it is done in eddy-covariance.
  #-- create a series of dates
  date.vec <- seq.Date(from=min(dat.hr.gf$Date),to=max(dat.hr.gf$Date),by="day")
  
  #-- estimate R-Tref and Tref for each date for each chamber
  #lagDates <- 3 # establish how many prior days to include
  RTdat <- expand.grid(Date=date.vec,chamber=levels(dat.hr.gf$chamber))
  RTdat$Tref <- RTdat$R_Tref <- NA
  
  
  print("Partitioning Net CO2 fluxes into GPP and Ra")
  #- set up progress bar to track that this is working
  pb <- txtProgressBar(min = 0, max = nrow(RTdat), style = 3)
  
  for (i in 1:nrow(RTdat)){
    #- trial a data.table alternative to speed this up. The filter thing was actually slower.
    
    #dat <- dplyr::filter(dat.hr.gf,chamber==RTdat$chamber[i],Date <= RTdat$Date[i], Date >= (RTdat$Date[i]-lagDates),period =="Night")
    #RTdat$Tref[i] <- mean(dat$Tair_al,na.rm=T)
    #RTdat$R_Tref[i] <- mean(dat$FluxCO2_g,na.rm=T)
    
    inds <- which(dat.hr.gf$chamber==RTdat$chamber[i] & dat.hr.gf$Date <= RTdat$Date[i] & dat.hr.gf$Date >= (RTdat$Date[i]-lagdates) & dat.hr.gf$period =="Night" )
    RTdat$Tref[i] <- mean(dat.hr.gf$Tair_al[inds],na.rm=T)
    RTdat$R_Tref[i] <- mean(dat.hr.gf$FluxCO2_g[inds],na.rm=T)
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  
  RTdat$Tref_K <- with(RTdat,Tref+273.15)
  
  #-- merge these reference data into the gap-filled flux dataframe to estimate Ra during the daytime, and hence GPP
  dat.hr.gf3 <- merge(dat.hr.gf,RTdat,by=c("Date","chamber"))
  #dat.hr.gf3$Ra_est <- with(dat.hr.gf3,R_Tref*Q10^((Tair_al-Tref)/10)) # estimate respiration rate. This is a negative number.
  dat.hr.gf3$Ra_est <- with(dat.hr.gf3,R_Tref*exp((Ea*1000/(rvalue*Tref_K))*(1-Tref_K/(Tair_al+273.15)))) # estimate respiration rate. This is a negative number.
  dat.hr.gf3$Ra_est <- ifelse(dat.hr.gf3$period=="Day",
                              dat.hr.gf3$Ra_est-leafRreduction*leafRtoTotal*dat.hr.gf3$Ra_est,
                              dat.hr.gf3$Ra_est) # estimate respiration rate. This is a negative number. If it's day, subtract 30% from the leaf R fraction
  

  dat.hr.gf3$GPP <- ifelse(dat.hr.gf3$period=="Night",0,dat.hr.gf3$FluxCO2_g-dat.hr.gf3$Ra_est)
  dat.hr.gf3$Ra <- ifelse(dat.hr.gf3$period=="Night",dat.hr.gf3$FluxCO2_g,dat.hr.gf3$Ra_est)
  
  return(dat.hr.gf3)
}
#----------------------------------------------------------------------------------------------------------------














#----------------------------------------------------------------------------------------------------------------
#-- Plots the net CO2 flux and its partitioning into GPP and Ra. Defaults to an example week for chamber 7
plotPartitionedFluxes <- function(dat.hr.gf3=dat.hr.gf3,ch_toplot="C07",startDate="2014-3-22",endDate="2014-3-27",export=F){
  
  startDate <- as.Date(startDate)
  endDate <- as.Date(endDate)
  
  #-- plot an example week, with PAR, Tair, and Cfluxes
  # windows(30,15);
  par(cex.lab=1.5,mar=c(0,0,0,0),oma=c(7,7,2,2),cex.axis=1.5,las=1)
  layout(matrix(c(1,2,3), 3, 1, byrow = TRUE), 
         widths=c(1,1,1), heights=c(1,1,3))
  plotBy(PAR~DateTime,xaxt="n",yaxt="n",data=subset(dat.hr.gf3,chamber==ch_toplot & as.Date(DateTime)>=startDate & 
                                              as.Date(DateTime)<=endDate),ylab="PPFD",lty=1,ylim=c(0,2000),type="l",legend=F,col="black")
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  mtext(text="PPFD",side=2,outer=F,line=4,cex=1.5,las=0)
  plotBy(Tair_al~DateTime,xaxt="n",yaxt="n",data=subset(dat.hr.gf3,chamber==ch_toplot & as.Date(DateTime)>=as.Date("2014-3-22") & 
                                                  as.Date(DateTime)<=as.Date("2014-3-27")),ylab="Tair",lty=1,ylim=c(11,34),type="l",legend=F,col="black")
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  mtext(text="Air T",side=2,outer=F,line=4,cex=1.5,las=0)
  
  
  plotBy(GPP~DateTime,xaxt="n",yaxt="n",data=subset(dat.hr.gf3,chamber==ch_toplot & as.Date(DateTime)>=startDate & 
                                              as.Date(DateTime)<=endDate),ylab="CO2 flux (g hr-1)",pch=16,cex=2,ylim=c(-2,10),type="b",legend=F,col="green")
  plotBy(Ra_est~DateTime,xaxt="n",yaxt="n",data=subset(dat.hr.gf3,chamber==ch_toplot & as.Date(DateTime)>=startDate & 
                                                 as.Date(DateTime)<=endDate),pch=16,col="red",type="b",cex=2,add=T,legend=F)
  plotBy(FluxCO2_g~DateTime,xaxt="n",yaxt="n",data=subset(dat.hr.gf3,chamber==ch_toplot & as.Date(DateTime)>=startDate & 
                                                    as.Date(DateTime)<=endDate),type="b",pch=16,col="black",cex=2,add=T,legend=F)
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  
  mtext(text=expression(CO[2]~flux~(gC~hr^-1)),side=2,outer=F,line=4,cex=1.5,las=0)
  legend("topright",c("GPP","Measured Net Flux","Ra"),pch=16,col=c("green","black","red"),cex=2)
  abline(0,0)
  axis.POSIXct(side=1,at=seq.POSIXt(from=as.POSIXct(startDate),to=as.POSIXct(endDate+1),by="day"),format="%F",cex.axis=2,line=1,tick=F)
  axis.POSIXct(side=1,at=seq.POSIXt(from=as.POSIXct(startDate),to=as.POSIXct(endDate+1),by="day"),format="%F",cex.axis=2,labels=F,tick=T)
  
  mtext(text="Date",side=1,outer=F,line=4,cex=1.5)
  
  if(export==T){dev.copy2pdf(file="output/FigureS2.pdf")}
}
#----------------------------------------------------------------------------------------------------------------






#----------------------------------------------------------------------------------------------------------------
#- function to return daily sums of GPP and Ra from hourly data
returnCUE.day <- function(dat=dat.hr.p){
  
  #- calculate leaf-area specific rates of GPP, convert to umol CO2 m-2 s-1
  dat$GPP_la <- with(dat,GPP/leafArea)
  dat$GPP_la_umol <- with(dat,GPP_la/12*1*10^6/60/60)
  dat$PAR_mol <- dat$PAR*60*60*1*10^-6

  
  
  #- create a date variable that moves the window for "night" observations, to sum observations for the night period following a day period
  dat$Date2 <- as.Date(dat$DateTime)
  dat$hour <- hour(dat$DateTime)
  earlynights <- which(dat$hour < 12 & dat$period =="Night")
  dat$Date2[earlynights] <- dat$Date2[earlynights]-1
  
  #- create daily sums 
  dat.day <- dplyr::summarize(group_by(dat,Date2,chamber,T_treatment,period),
                              GPP = sum(GPP,na.rm=T),
                              Ra = sum(Ra,na.rm=T),
                              FluxCO2_g=sum(FluxCO2_g,na.rm=T),
                              Tair_al=mean(Tair_al,na.rm=T),
                              VPDair=max(VPDair,na.rm=T),
                              PAR=sum(PAR_mol,na.rm=T),
                              leafArea=mean(leafArea,na.rm=T))
  dat.day <- as.data.frame(dat.day)
  dat.day <- dat.day[with(dat.day,order(Date2,chamber)),]
  names(dat.day)[1] <- "Date"
  
  Tair_day <- dplyr::summarize(group_by(dat,Date2,chamber,T_treatment),
                               Tair_al=mean(Tair_al,na.rm=T))
  names(Tair_day)[4] <- "Tair_24hrs"                   
  
  dat.day_tair <- as.data.frame(Tair_day)
  dat.day <- as.data.frame(dat.day)
  
  #- merge night and day data
  dat.day2d <- subset(dat.day,period=="Day")[,c("Date","chamber","T_treatment","GPP","Ra","FluxCO2_g","PAR","Tair_al","VPDair","leafArea")]
  names(dat.day2d)[4:8] <- c("GPP","Raday","Cgain","PAR","T_day")
  dat.day2n <- subset(dat.day,period=='Night')[,c("Date","chamber","T_treatment","Ra","FluxCO2_g","Tair_al","VPDair")]
  names(dat.day2n)[4:7] <- c("Ranight","Closs","T_night","VPDair_night")
  names(dat.day2n)[1] <- c("Date")
  
  #- merge data such that the nightly data are one day ahead of the daily data (i.e., does tonight's respiration depend on today's photosynthesis?)
  cue.day1 <- subset(merge(dat.day2d,dat.day2n,by=c("Date","chamber","T_treatment")))
  cue.day <- merge(cue.day1,dat.day_tair,by.x=c("Date","chamber","T_treatment"),by.y=c("Date2","chamber","T_treatment"))
  cue.day$Ra <- with(cue.day,-1*(Raday+Ranight))
  cue.day$RtoA <- with(cue.day,Ra/GPP)
  cue.day$CUE <- with(cue.day,(1-(Ra/GPP)))
  cue.day$GPP_la <- with(cue.day,GPP/leafArea)
  cue.day$Ra_la <- with(cue.day,Ra/(leafArea)) #g C m-2 day-1
  
  # remove a few days from the beginning of the dataset in C09. Problem with flux data.
  toremove <- which(cue.day$Ra_la < 1.5 & cue.day$chamber=="C07" & cue.day$Date < as.Date("2013-10-3"))
  cue.day <- cue.day[-toremove,] 
  
  
  
  
  #- merge in the key
  key <- read.csv("data/WTC_TEMP_CM_TREATKEY_20121212-20140528_L1_v1.csv")
  cue.day$T_treatment <- ifelse(cue.day$T_treatment=="ambient","ambient","elevated")
  cue.day2 <- merge(cue.day,key,by=c("chamber","Date","T_treatment"))
  
  cue.day.trt <- summaryBy(GPP+GPP_la+Ra+Ra_la+RtoA+PAR+T_day+T_night+Tair_24hrs+VPDair+VPDair_night~Date+T_treatment+Water_treatment,data=cue.day2,FUN=c(mean,standard.error))
  
  return(list(cue.day2,cue.day.trt))
}
#----------------------------------------------------------------------------------------------------------------





#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#- plot environment, CUE, GPP and Ra (normalized by leaf area) over time
plotPAR_AirT_CUE_GPP_Ra <- function(cue.day.trt=cue.day.trt,export=F,lwidth=2.5){
  
  # windows(25,20);
  par(mfrow=c(4,1),mar=c(0,2,0,2),oma=c(6,8,0,6),cex.axis=2,las=1)
  palette(c("black","red"))
  
  #plot temperature
  plotBy(Tair_24hrs.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",
         ylim=c(12,45),
         ylab="",cex.lab=1.5,lwd=3,xaxt="n",yaxt="n")
  mtext(expression(atop(T[air],
                        (degree*C))),side=2,las=0,cex=1.2,line=3)
  legend("topleft",c("A","W"),lty=c(1,1),lwd=3,col=c("black","red"),ncol=2,bty="n",seg.len=3,cex=1.2)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "months"),
            labels = F, tcl = 0.5)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "weeks"),
            labels = F, tcl = 0.2)
  magaxis(side=c(2),labels=c(F),frame.plot=T)
  axis(side=2,at=c(20,30,40),tick=F)
  #-overlay par
  par(new=T)
  plot(PAR.mean~Date,type="l",data=subset(cue.day.trt,Water_treatment=="control" & T_treatment=="ambient"),
       ylim=c(-40,70),col="grey",lwd=2,xaxt="n",yaxt="n")
  magaxis(side=4,labels=F)
  axis(side=4,at=c(0,20,40,60),tick=F)
  mtext(expression(atop(PPFD,
                        (mol~d^-1))),side=4,las=0,cex=1.2,line=6.5)
  legend("topright",c("PPFD"),lty=c(1),lwd=3,col=c("grey"),ncol=1,bty="n",seg.len=3,cex=1.2)
  legend("bottomright","a",bty="n",inset=-0.002,cex=1.5)
  
  
  #-- plot Ra per unit leaf area
  plotBy(Ra_la.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",ylim=c(0.25,2.8),
         ylab="Ra (gC d-1)",cex.lab=1.5,lwd=3,xaxt="n",yaxt="n")
  a <- subset(cue.day.trt,T_treatment=="ambient")
  e <- subset(cue.day.trt,T_treatment=="elevated")
  adderrorbars(x=a$Date-0.25,y=a$Ra_la.mean,SE=a$Ra_la.standard.error,direction="updown",barlen=0,col="black")
  adderrorbars(x=e$Date+0.25,y=e$Ra_la.mean,SE=e$Ra_la.standard.error,direction="updown",barlen=0,col="red")
  plotBy(Ra_la.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",ylim=c(0,90),
         cex.lab=1.5,lwd=lwidth,add=T)
  mtext(expression(atop(R[a],
                        (gC~m^-2~d^-1))),side=2,las=0,cex=1.2,line=3)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "months"),
            labels = F, tcl = 0.5,las=2)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "weeks"),
            labels = F, tcl = 0.2)
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  legend("bottomright","b",bty="n",inset=-0.002,cex=1.5)
  
  #-- plot GPP per unit leaf area
  plotBy(GPP_la.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",ylim=c(0.5,7.7),
         ylab="",cex.lab=1.5,lwd=3,xaxt="n",yaxt="n")
  a <- subset(cue.day.trt,T_treatment=="ambient")
  e <- subset(cue.day.trt,T_treatment=="elevated")
  adderrorbars(x=a$Date-.25,y=a$GPP_la.mean,SE=a$GPP_la.standard.error,direction="updown",barlen=0,col="black")
  adderrorbars(x=e$Date+0.25,y=e$GPP_la.mean,SE=e$GPP_la.standard.error,direction="updown",barlen=0,col="red")
  plotBy(GPP_la.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",ylim=c(0,90),
         lwd=lwidth,add=T)
  mtext(expression(atop(GPP,
                        (gC~m^-2~d^-1))),side=2,las=0,cex=1.2,line=3)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "months"),
            labels = F, tcl = 0.5)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "weeks"),
            labels = F, tcl = 0.2)
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  legend("bottomright","c",bty="n",inset=-0.002,cex=1.5)
  

  #-- plot RtoA
  plotBy(RtoA.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",ylim=c(0,1),
         ylab="",cex.lab=1.5,lwd=1,xaxt="n",yaxt="n")
  a <- subset(cue.day.trt,T_treatment=="ambient")
  e <- subset(cue.day.trt,T_treatment=="elevated")
  adderrorbars(x=a$Date-0.25,y=a$RtoA.mean,SE=a$RtoA.standard.error,direction="updown",barlen=0,col="black")
  adderrorbars(x=e$Date+0.25,y=e$RtoA.mean,SE=e$RtoA.standard.error,direction="updown",barlen=0,col="red")
  plotBy(RtoA.mean~Date|T_treatment,data=subset(cue.day.trt,Water_treatment=="control"),cex=1,pch=16,legend=F,type="l",ylim=c(0,90),
         cex.lab=1.5,lwd=lwidth,add=T)
  mtext(expression(R[a]/GPP),side=2,las=0,cex=1.2,line=5)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "months"),
            labels = T, tcl = 0.5,las=2)
  axis.Date(1, at = seq(as.Date("2013/9/1"), as.Date("2014/6/1"), "weeks"),
            labels = F, tcl = 0.2)
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T)
  legend("bottomright","d",bty="n",inset=-0.002,cex=1.5)
  

  if(export==T) dev.copy2pdf(file="output/Figure4.pdf")
}
#--------------------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------------------------
# plot climate drivers of GPP, Ra, and Ra/GPP (Figure 5)

plotGPP_Ra_CUE_metdrivers <- function(cue.day=cue.day,export=T,shading=0.5,parcut=35){
  
  
  palette(alpha(c("black","red"),0.5))
  
  # windows(20,30);
  par(mfrow=c(3,2),mar=c(0,0,0,0),oma=c(7,7,1,3),las=1,cex.lab=1.8,cex.axis=1.2)
  
  #- plot GPP
  smoothplot(PAR, GPP_la, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0.5,9),
             data=cue.day, kgam=5, axes=F)
  title(ylab=expression(GPP~(gC~m^-2~d^-1)),outer=T,adj=0.95)
  magaxis(side=1:4,labels=c(0,1,0,0))
  legend("topright","a",bty="n",inset=-0.01,cex=1.2)
  
  smoothplot(Tair_24hrs, GPP_la, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0.5,9),
             data=subset(cue.day,PAR>parcut), kgam=5,axes=F)
  magaxis(side=1:4,labels=c(0,0,0,0))
  legend("topright","b",bty="n",inset=0.01,cex=1.2)
  title(ylab=expression(GPP~(gC~m^-2~d^-1)),outer=T,adj=0.95)
  magaxis(side=1:4,labels=c(0,0,0,1))
  
  
  #- plot Ra
  smoothplot(PAR, Ra_la, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0.2,4),
             data=cue.day, kgam=5, axes=F)
  title(ylab=expression(R[a]~(gC~m^-2~d^-1)),outer=T,adj=0.5)
  magaxis(side=1:4,labels=c(0,1,0,0))
  legend("topright","c",bty="n",inset=-0.01,cex=1.2)
  
  smoothplot(Tair_24hrs, Ra_la, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0.2,4),
             data=subset(cue.day,PAR>parcut), kgam=5, axes=F)
  magaxis(side=1:4,labels=c(0,0,0,0))
  legend("topright","d",bty="n",inset=0.01,cex=1.2)
  title(ylab=expression(R[a]~(gC~m^-2~d^-1)),outer=T,adj=0.5)
  magaxis(side=1:4,labels=c(0,0,0,1))
  
  #- plot CUE
  
  smoothplot(PAR, RtoA, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0,1),
             data=cue.day, kgam=5, axes=F)
  title(ylab=expression(R[a]/GPP),outer=T,adj=0.15)
  magaxis(side=1:4,labels=c(1,1,0,0))
  legend("topright","e",bty="n",inset=-0.01,cex=1.2)
  
  smoothplot(Tair_24hrs, RtoA, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0,1),
             data=subset(cue.day,PAR>parcut), kgam=5, axes=F)
  magaxis(side=1:4,labels=c(1,0,0,0))
  legend("topright","f",bty="n",inset=0.01,cex=1.2)
  title(ylab=expression(R[a]/GPP),outer=T,adj=0.15)
  magaxis(side=1:4,labels=c(1,0,0,1))
  
  title(xlab=expression(PPFD~(mol~d^-1)),outer=T,adj=0.08)
  
  title(xlab=expression(T[air]~(degree*C*", 24-"*h~mean)),outer=T,adj=0.95)
  
  if(export==T) dev.copy2pdf(file="output/Figure5.pdf")
}
#----------------------------------------------------------------------------------------------------------------




plotNetC_metdrivers <- function(cue.day=cue.day,export=T,shading=0.5,parcut=35){
  
  
  palette(alpha(c("black","red"),0.5))
  
  # windows(30,30);
  par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(7,7,1,3),las=1,cex.lab=1.8,cex.axis=1.2)
  
  #- plot NetC
  smoothplot(PAR, netC, T_treatment,polycolor=NA,linecol=NA,#polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0,100),
             data=cue.day, kgam=3, axes=F)
  title(ylab=expression(Raw~net~C~gain~(gC~d^-1)),outer=T,adj=0.95)
  magaxis(side=1:4,labels=c(0,1,0,0))
  legend("topright","a",bty="n",inset=-0.01,cex=1.2)
  
  smoothplot(Tair_24hrs, netC, T_treatment,polycolor=NA,linecol=NA,#polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0,100),
             data=subset(cue.day,PAR>parcut), kgam=5,axes=F)
  magaxis(side=1:4,labels=c(0,0,0,0))
  legend("topright","b",bty="n",inset=0.01,cex=1.2)
  magaxis(side=1:4,labels=c(0,0,0,1))
  
  #- plot NetC_la
  smoothplot(PAR, netC_la, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0,9),
             data=cue.day, kgam=5, axes=F)
  title(ylab=expression(Net~C~gain~(gC~m^-2~d^-1)),outer=T,adj=0.05)
  magaxis(side=1:4,labels=c(0,1,0,0))
  legend("topright","c",bty="n",inset=-0.01,cex=1.2)
  
  smoothplot(Tair_24hrs, netC_la, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0,9),
             data=subset(cue.day,PAR>parcut), kgam=5,axes=F)
  magaxis(side=1:4,labels=c(0,0,0,0))
  legend("topright","d",bty="n",inset=0.01,cex=1.2)
  magaxis(side=1:4,labels=c(0,0,0,1))
  
  
  
  title(xlab=expression(PPFD~(mol~d^-1)),outer=T,adj=0.15)
  
  title(xlab=expression(T[air]~(degree*C*", 24-"*h~mean)),outer=T,adj=0.9)
  
  if(export==T) dev.copy2pdf(file="output/FigureS7.pdf")
}




plotClosstogainratio_metdrivers <- function(cue.day=cue.day,export=T,shading=0.5,parcut=35){
  
  
  palette(alpha(c("black","red"),0.5))
  
  # windows(50,30);
  par(mfrow=c(1,2),mar=c(0,0,0,0),oma=c(7,7,1,3),las=1,cex.lab=1.8,cex.axis=1.2)
  
  #- plot raw ratio of C loss to C gain
  smoothplot(PAR, Closs_Cgain, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0,1),
             data=cue.day, kgam=5, axes=F)
  title(ylab=expression(C~loss~"/"~C~gain),outer=T,adj=0.5)
  magaxis(side=1:4,labels=c(1,1,0,0))
  legend("topright","a",bty="n",inset=0.01,cex=1.2)
  
  smoothplot(Tair_24hrs, Closs_Cgain, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",
             ylim=c(0,1),
             data=subset(cue.day,PAR>parcut), kgam=5,axes=F)
  magaxis(side=1:4,labels=c(0,0,0,0))
  legend("topright","b",bty="n",inset=0.01,cex=1.2)
  magaxis(side=1:4,labels=c(1,0,0,1))
  
  
  
  
  title(xlab=expression(PPFD~(mol~d^-1)),outer=T,adj=0.15)
  
  title(xlab=expression(T[air]~(degree*C*", 24-"*h~mean)),outer=T,adj=0.9)
  
  if(export==T) dev.copy2pdf(file="output/FigureS8.pdf")
}






#----------------------------------------------------------------------------------------------------------------
# plot temperature drivers of and Ra/GPP

plot_CUE_temp_drought <- function(cue.day=cue.day,export=T,shading=0.7,parcut=35){
  
  # windows(30,20)
  par(mfrow=c(1,2),mar=c(1,1,1,1),oma=c(6,7,1,3))
  cue.day$combo <- factor(paste(cue.day$T_treatment,cue.day$Water_treatment))
  
  #- plot ambient temperature
  smoothplot(Tair_24hrs, RtoA, combo,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             linecols=c("blue","red"),pointcols=c("blue","red"),
             random="chamber",
             ylim=c(0,1),xlim=c(10,30),
             data=subset(cue.day,PAR>parcut & T_treatment=="ambient" & Date > as.Date("2014-01-1")), kgam=5, axes=F)
  legend("top","Ambient temperature",bty="n",inset=0.01,cex=2)
  legend(x=10,y=0.8,c("Well-watered","Dry-down"),lty=1,lwd=2,bty="n",col=c("blue","red"),cex=2)
  
  title(ylab=expression(R[a]/GPP),outer=T,adj=0.5,cex.lab=3)
  magaxis(side=1:4,labels=c(1,1,0,0),las=1,cex.axis=2)
  
  #- plot warmed
  smoothplot(Tair_24hrs, RtoA, combo,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             linecols=c("blue","red"),pointcols=c("blue","red"),
             random="chamber",
             ylim=c(0,1),xlim=c(10,30),
             data=subset(cue.day,PAR>parcut & T_treatment=="elevated" & Date > as.Date("2014-01-1")), kgam=5, axes=F)
  magaxis(side=1:4,labels=c(1,0,0,1),las=1,cex.axis=2)
  legend("top","Warmed temperature",bty="n",inset=0.01,cex=2)
  
  title(xlab=expression(T[air]~(degree*C*", 24-"*h~mean)),outer=T,adj=0.5,cex.lab=2)
  
  if(export==T) dev.copy2pdf(file="output/RatoGPP_drought.pdf")
}
#----------------------------------------------------------------------------------------------------------------









#----------------------------------------------------------------------------------------------------------------
#' Function for smoothplots of GAMs. 
fitgam <- function(X,Y,dfr, k=-1, R=NULL){
  dfr$Y <- dfr[,Y]
  dfr$X <- dfr[,X]
  if(!is.null(R)){
    dfr$R <- dfr[,R]
    model <- 2
  } else model <- 1
  dfr <- droplevels(dfr)
  
  
  if(model ==1){
    g <- gam(Y ~ s(X, k=k), data=dfr)
  }
  if(model ==2){
    g <- gamm(Y ~ s(X, k=k), random = list(R=~1), data=dfr)
  }
  
  return(g)
}


#' Plot a generalized additive model
#' @param x Variable for X axis (unquoted)
#' @param y Variable for Y axis (unquoted)
#' @param data Dataframe containing x and y
#' @param kgam the \code{k} parameter for smooth terms in gam.
#' @param random An optional random effect (quoted)
#' @param log Whether to add log axes for x or y (but no transformations are done).
#' @param fitoneline Whether to fit only 
smoothplot <- function(x,y,g=NULL,data,
                       fittype=c("gam","lm"),
                       kgam=4,
                       random=NULL,
                       randommethod=c("lmer","aggregate"),
                       log="",
                       fitoneline=FALSE,
                       pointcols=NULL,
                       linecols=NULL, 
                       xlab=NULL, ylab=NULL,
                       polycolor=alpha("lightgrey",0.7),
                       axes=TRUE,
                       ...){
  
  fittype <- match.arg(fittype)
  randommethod <- match.arg(randommethod)
  
  if(!is.null(substitute(g))){
    data$G <- as.factor(eval(substitute(g),data))
  } else {
    fitoneline <- TRUE
    data$G <- 1
  }
  data$X <- eval(substitute(x),data)
  data$Y <- eval(substitute(y),data)
  data <- droplevels(data)
  
  data <- data[!is.na(data$X) & !is.na(data$Y) & !is.na(data$G),]
  
  if(is.null(pointcols))pointcols <- palette()
  if(is.null(linecols))linecols <- palette()
  
  if(is.null(xlab))xlab <- substitute(x)
  if(is.null(ylab))ylab <- substitute(y)
  
  # If randommethod = aggregate, average by group and fit simple gam.
  if(!is.null(random) && randommethod == "aggregate"){
    data$R <- data[,random]
    
    data <- summaryBy(. ~ R, FUN=mean, na.rm=TRUE, keep.names=TRUE, data=data,
                      id=~G)
    R <- NULL
  }
  
  
  if(!fitoneline){
    
    d <- split(data, data$G)
    
    if(fittype == "gam"){
      fits <- lapply(d, function(x)try(fitgam("X","Y",x, k=kgam, R=random)))
      if(!is.null(random))fits <- lapply(fits, "[[", "gam")
    } else {
      fits <- lapply(d, function(x)lm(Y ~ X, data=x))
    }
    hran <- lapply(d, function(x)range(x$X, na.rm=TRUE))
  } else {
    if(fittype == "gam"){
      fits <- list(fitgam("X","Y",data, k=kgam, R=random))
      if(!is.null(random))fits <- lapply(fits, "[[", "gam")
    } else {
      fits <- list(lm(Y ~ X, data=data))
    }
    hran <- list(range(data$X, na.rm=TRUE))
    
  }
  
  with(data, plot(X, Y, xaxt="n",yaxt="n", pch=16, col=pointcols[G],
                  xlab=xlab, ylab=ylab, ...))
  
  if(axes){
    if(log=="xy")magaxis(side=1:2, unlog=1:2)
    if(log=="x"){
      magaxis(side=1, unlog=1)
      axis(2)
      box()
    }
    if(log=="y"){
      magaxis(side=2, unlog=2)
      axis(1)
      box()
    }
    if(log==""){
      axis(1)
      axis(2)
      box()
    }
  }
  
  for(i in 1:length(fits)){
    
    if(fittype == "gam"){
      nd <- data.frame(X=seq(hran[[i]][1], hran[[i]][2], length=101))
      if(!inherits(fits[[i]], "try-error")){
        p <- predict(fits[[i]],nd,se.fit=TRUE)
        addpoly(nd$X, p$fit-2*p$se.fit, p$fit+2*p$se.fit, col=polycolor[i])
        lines(nd$X, p$fit, col=linecols[i], lwd=2)
      }
    }
    if(fittype == "lm"){
      pval <- summary(fits[[i]])$coefficients[2,4]
      LTY <- if(pval < 0.05)1 else 5
      predline(fits[[i]], col=linecols[i], lwd=2, lty=LTY)
    }
  }
  
  return(invisible(fits))
}


addpoly <- function(x,y1,y2,col=alpha("lightgrey",0.7),...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x,rev(x)), c(y1, rev(y2)), col=col, border=NA,...)
}

predline <- function(fit, from=NULL, to=NULL, col=alpha("lightgrey",0.7), ...){
  
  if(is.null(from))from <- min(fit$model[,2], na.rm=TRUE)
  if(is.null(to))to <- max(fit$model[,2], na.rm=TRUE)
  
  newdat <- data.frame(X = seq(from,to, length=101))
  names(newdat)[1] <- names(coef(fit))[2]
  
  pred <- as.data.frame(predict(fit, newdat, se.fit=TRUE, interval="confidence")$fit)
  
  addpoly(newdat[[1]], pred$lwr, pred$upr, col=col)
  
  #ablinepiece(fit, from=from, to=to, ...)
  lines(pred$fit~newdat[,1])
}

#'@title Add a line to a plot
#'@description As \code{abline}, but with \code{from} and \code{to} arguments. 
#'If a fitted linear regression model is used as asn argument, it uses the min and max values of the data used to fit the model.
#'@param a Intercept (optional)
#'@param b Slope (optional)
#'@param reg A fitted linear regression model (output of \code{\link{lm}}).
#'@param from Draw from this X value
#'@param to Draw to this x value
#'@param \dots Further parameters passed to \code{\link{segments}}
#'@export
ablinepiece <- function(a=NULL,b=NULL,reg=NULL,from=NULL,to=NULL,...){
  
  # Borrowed from abline
  if (!is.null(reg)) a <- reg
  
  if (!is.null(a) && is.list(a)) {
    temp <- as.vector(coefficients(a))
    from <- min(a$model[,2], na.rm=TRUE)
    to <- max(a$model[,2], na.rm=TRUE)
    
    if (length(temp) == 1) {
      a <- 0
      b <- temp
    }
    else {
      a <- temp[1]
      b <- temp[2]
    }
  }
  
  segments(x0=from,x1=to,
           y0=a+from*b,y1=a+to*b,...)
  
}

#----------------------------------------------------------------------------------------------------------------







#----------------------------------------------------------------------------------------------------------------
#- plot heat-map style hexagons for GPP relative to air temperaure and PAR (i.e., Figure 6)
plotGPP_hex <- function(dat=dat.hr.p,export=F,shading=0.7){
  
  #- Calculate GPP per unit leaf area, convert to umol CO2 m-2 s-1
  dat$GPP_la <- with(dat,GPP/leafArea)
  dat$GPP_la_umol <- with(dat,GPP_la/12*1*10^6/60/60)

  #creates a scale of colors
  myColorRamp <- function(colors, values) {
    #v <- (values - min(values))/diff(range(values))
    v <- (values - -1)/21 
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
  }
  a <- subset(dat,PAR>10 & Water_treatment=="control" & T_treatment=="ambient")[,c("Tair_al","PAR","GPP","GPP_la","GPP_la_umol")]
  e <- subset(dat,PAR>10 & Water_treatment=="control" & T_treatment=="elevated")[,c("Tair_al","PAR","GPP","GPP_la","GPP_la_umol")]
  all <- subset(dat,PAR>10 & Water_treatment=="control")[,c("T_treatment","chamber","Date","Tair_al","PAR","GPP","GPP_la","GPP_la_umol")]
  all <- subset(all,GPP_la_umol < 20)
  
  #- create hex bins for a, e and all
  ha <- hexbin(a$Tair_al, a$PAR,xbins=30,IDs=TRUE)
  he <- hexbin(e$Tair_al, e$PAR,xbins=30,IDs=TRUE)
  hall <- hexbin(all$Tair_al, all$PAR,xbins=30,IDs=TRUE)

  #average values of points inside hexbins 
  meanHexBina<-data.frame(mean=hexTapply(ha, a$GPP_la_umol, mean)) 
  meanHexBine<-data.frame(mean=hexTapply(he, e$GPP_la_umol, mean)) 
  meanHexBinall<-data.frame(mean=hexTapply(hall, all$GPP_la_umol, mean)) 

  
  colsa <- myColorRamp(c("blue","green","yellow", "red"), meanHexBina$mean)
  colse <- myColorRamp(c("blue","green","yellow", "red"), meanHexBine$mean)
  colsall <- myColorRamp(c("blue","green","yellow", "red"), meanHexBinall$mean)
  

  #-- plot Ambient
  # windows()
  
  ## setup coordinate system of the plot
  par(cex.lab=2,cex.axis=2)
  Pa <- plot(hall, type="n",legend=FALSE,xlab="",ylab="",main="")
  
  ##add hexagons (in the proper viewport):
  pushHexport(Pa$plot.vp)
  
  #plots hexbins based on colors of third column
  #grid.hexagons(hall, style= "lattice", border = gray(.9), pen = colsall,  minarea = 1, maxarea = 1)
  
  grid.hexagons(ha, style= "lattice", border = gray(.9), pen = colsa,  minarea = 1, maxarea = 1)
  if(export==T) dev.copy2pdf(file="output/Figure6a.pdf")
  
  #-- plot elevated
  # windows()
  
  ## setup coordinate system of the plot
  Pe <- plot(hall, type="n",legend=FALSE,xlab="",ylab="",main="")
  
  ##add hexagons (in the proper viewport):
  pushHexport(Pe$plot.vp)
  
  #plots hexbins based on colors of third column
  grid.hexagons(he, style= "lattice", border = gray(.9), pen = colse,  minarea = 1, maxarea = 1)
  if(export==T) dev.copy2pdf(file="output/Figure6b.pdf")
  
  
  
  #- get five hexagons, add them as a legend
  # windows()
  values <- seq(0,20,length=6)
  colors <- myColorRamp(c("blue","green","yellow", "red"), values)
  xloc <- rep(1,6)
  yloc <- c(1,2,3,4,5,6)
  plot(yloc~xloc,pch=18,cex=5,col=colors,ylim=c(0,10),xaxt="n",yaxt="n",xlab="",ylab="")
  text(xloc+0.1,yloc,labels=values,cex=1.5)
  #dev.copy2pdf(file="/output/Figure6_legend.pdf")
  
  
  
  #-- subset data to a par range, plot Tresponse
  toplot <- subset(all, PAR >1200 & PAR <1500 & GPP_la_umol>0)
  
  #- plot smoothplots
  # windows();
  par(las=1)
  smoothplot(Tair_al, GPP_la_umol, T_treatment,pointcols=c(alpha("black",0.3),alpha("red",0.3)),
             linecol=c("black","red"),polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",cex=1,main="",
             xlim=c(15,45),ylim=c(0,25),xlab="",ylab="",
             data=toplot, kgam=4,axes=F)
  legend("topright",pch=16,legend=c("ambient","warmed"),col=c(alpha("black",0.3),alpha("red",0.3)))
  box();axis(side=1,labels=T);axis(side=2,labels=T);axis(side=4,labels=F)
  if(export==T) dev.copy2pdf(file="output/Figure6c.pdf")
  
}
#----------------------------------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------------------------
plotAnet_met_diurnals <- function(export=T,lsize=2,size=2,printANOVAs=F){
  
  #- make multipanel plot of diurnals 
  diurnal <- read.csv("data/WTC_TEMP_CM_GX-DIURNAL_20130710-20140220_L1_v2.csv")
  diurnal$DateTime <- as.POSIXct(diurnal$DateTime,format="%Y-%m-%d %T",tz="GMT")
  diurnal$Date <- as.Date(diurnal$DateTime)
  diurnal$Hour <- hour(diurnal$DateTime)
  diurnal$timepoint <- as.factor(diurnal$timepoint)
  
  #- do statistical tests on each date for photo
  if(printANOVAs==T){
    d.top <- subset(diurnal,position=="top")
    d.l <- split(d.top,d.top$Date)
    atest <- gtest <- list()
    for(i in 1:length(d.l)){
      dat <- d.l[[i]]
      
      print(paste("***"," Date = ",dat$Date[1],sep=""))
      print("Anet")
      atest[[i]] <- lme(Photo~T_treatment*timepoint,random=list(~1|chamber),data=dat)
      print(anova(atest[[i]]))
      
      print("Cond")
      gtest[[i]] <- lme(Cond~T_treatment*timepoint,random=list(~1|chamber),data=dat)
      print(anova(gtest[[i]]))
    }
  }
  
  #-- get treatment averages for plotting
  trtavgs <- summaryBy(Photo+Cond+Tleaf+VpdL+PARi+DateTime+Hour~ T_treatment+Date+timepoint,
                       data=subset(diurnal,position=="top"),FUN=c(mean,standard.error),na.rm=T)
  trtavgs.list <- split(trtavgs,trtavgs$Date)
  ylims <- c(0,25)
  xlims <- c(5,21)
  #xlims=c(0.5,6.5)
  # windows(40,40);
  par(mfrow=c(3,5),mar=c(0,0,0,0),oma=c(5,12,5,3))
  palette(c("black","red"))
  #- plot photo on each date
  for(i in 1:length(trtavgs.list)){
    toplot <- trtavgs.list[[i]]
    
    plotBy(Photo.mean~Hour.mean|T_treatment,data=toplot,legend=F,pch=16,xaxt="n",yaxt="n",type="b",cex=size,ylim=ylims,xlim=xlims,
           panel.first=adderrorbars(x=toplot$Hour.mean,y=toplot$Photo.mean,SE=toplot$Photo.standard.error,direction="updown"))
    if(i==1) magaxis(side=c(1,2,3,4),labels=c(0,1,0,0),las=1)
    if(i >1 & i <5) magaxis(side=c(1,2,3,4),labels=c(0,0,0,0),las=1)
    if(i==5)magaxis(side=c(1,2,3,4),labels=c(0,0,0,1),las=1)
    if(i==1)legend("topleft",c("A","W"),lty=1,lwd=lsize,bty="n",col=c("black","red"),cex=1.2)
    if(i==1)points(x=c(8,8),y=c(24.1,21.9),pch=16,col=c("black","red"),cex=size)
    legend("bottomleft",letters[i],cex=1.2,bty='n')
    if(i==1) title(ylab=expression(atop(A[net],
                                        ~(mu*mol~m^-2~s^-1))),xpd=NA,cex.lab=1.8)
    title(main=format(toplot$Date[1],"%Y-%b-%d"),xpd=NA,cex.lab=1.5,line=1)
    
  }
  
  
  #- plot gs on each date
  for(i in 1:length(trtavgs.list)){
    toplot <- trtavgs.list[[i]]
    
    plotBy(Cond.mean~Hour.mean|T_treatment,data=toplot,legend=F,pch=16,xaxt="n",yaxt="n",type="b",cex=size,ylim=c(0,0.31),xlim=xlims,
           panel.first=adderrorbars(x=toplot$Hour.mean,y=toplot$Cond.mean,SE=toplot$Cond.standard.error,direction="updown"))
    if(i==1) magaxis(side=c(1,2,3,4),labels=c(0,1,0,0),las=1)
    if(i >1 & i <5) magaxis(side=c(1,2,3,4),labels=c(0,0,0,0),las=1)
    if(i==5)magaxis(side=c(1,2,3,4),labels=c(0,0,0,1),las=1)
    legend("bottomleft",letters[i+5],cex=1.2,bty='n')
    if(i==1) title(ylab=expression(atop(g[s],
                                        ~(mol~m^-2~s^-1))),xpd=NA,cex.lab=1.8)
    
  }
  
  #- plot PAR, airT, and VPD on each date
  for(i in 1:length(trtavgs.list)){
    toplot <- trtavgs.list[[i]]
    
    #- plot Tleaf in black and red
    plotBy(Tleaf.mean~Hour.mean,data=subset(toplot,T_treatment=="ambient"),legend=F,pch=16,xaxt="n",yaxt="n",type="l",col=c("black"),
           cex=size,ylim=c(0,45),xlim=xlims,lwd=lsize)
    if(i==1)magaxis(side=c(1,2,3,4),labels=c(1,1,0,0),las=1,col="black")
    
    if(i >1 & i <5) magaxis(side=c(1,2,3,4),labels=c(1,0,0,0),las=1)
    if(i==5)magaxis(side=c(1,2,3,4),labels=c(1,0,0,0),las=1)
    
    if(i==1)legend("topleft",c("T","PPFD","VPD"),lty=1,lwd=lsize,bty="n",col=c("black","blue","red"),cex=1.2)
    if(i==1) title(ylab=expression(Environmental~variables),xpd=NA,cex.lab=1.8,line=7)
    
    #- plot PAR in blue
    par(new=T)
    plotBy(PARi.mean~Hour.mean,data=subset(toplot,T_treatment=="ambient"),legend=F,pch=16,xaxt="n",yaxt="n",type="l",col="blue",cex=size,ylim=c(0,2000),lwd=lsize,xlim=xlims)
    
    if(i==1)axis(2, ylim=c(0,8),lwd=1,line=1.8,col="blue",col.axis="blue",las=1)
    
    #- plot VPD in red
    par(new=T)
    plotBy(VpdL.mean~Hour.mean,data=subset(toplot,T_treatment=="ambient"),lty=1,lwd=lsize,
           legend=F,pch=16,xaxt="n",yaxt="n",type="l",col="red",cex=size,ylim=c(0,7),xlim=xlims)
    #- add VPD of warmed
    #plotBy(VpdL.mean~Hour.mean,data=subset(toplot,T_treatment=="elevated"),lty=2,lwd=lsize,add=T,
    #       legend=F,pch=16,xaxt="n",yaxt="n",type="l",col="red",cex=size,ylim=c(0,7),xlim=xlims)
    
    if(i==3)title(xlab="Hour",cex.lab=1.5,xpd=NA)
    if(i==1)axis(2, ylim=c(0,8),lwd=1,line=4.8,col="red",col.axis="red",las=1)
    legend("bottomleft",letters[i+10],cex=1.2,bty='n')
    
    
    
  }
  
  
  if(export==T)dev.copy2pdf(file="output/Figure7.pdf")
}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------






#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#- compare the direct diurnal measurements with the whole-tree fluxes for a specified focal date.
#-   Take some care here, as the flux data were frequently perterbed by investigators going into chambers on these dates
plotDiurnalvsWTC <- function(fluxdat=dat.hr.p,focaldate=as.Date("2013-12-4"),output=T){
  #------------------------------------
  #- get diurnal data
  diurnal <- read.csv("data/WTC_TEMP_CM_GX-DIURNAL_20130710-20140220_L1_v2.csv")
  diurnal$DateTime <- as.POSIXct(diurnal$DateTime,format="%Y-%m-%d %T",tz="GMT")
  diurnal$Date <- as.Date(diurnal$DateTime)
  diurnal$Hour <- hour(diurnal$DateTime)
  diurnal$timepoint <- as.factor(diurnal$timepoint)
  
  focaldate <- as.Date("2013-12-4")
  
  leaf.m <- summaryBy(Photo+Cond+Tleaf+VpdL+PARi+DateTime+Hour~ T_treatment+Date+timepoint,
                      data=subset(diurnal,position=="top" & Date==focaldate),FUN=c(mean,standard.error),na.rm=T)
  
  #- get teh flux data for that day
  canopyexample <- subset(dat.hr.p,Date==focaldate)
  canopyexample$Hour <- hour(canopyexample$DateTime)
  #- Calculate GPP per unit leaf area, convert to umol CO2 m-2 s-1
  canopyexample$GPP_la <- with(canopyexample,GPP/leafArea)
  canopyexample$GPP_la_umol <- with(canopyexample,GPP_la/12*1*10^6/60/60)
  canopy.m <- summaryBy(GPP_la_umol~ T_treatment+Hour,
                        data=subset(canopyexample,Hour>4 & Hour<22),FUN=c(mean,standard.error),na.rm=T)
  
  # windows();
  par(mar=c(5,7,1,3))
  ylims=c(0,25)
  xlims=c(5,20)
  size=2
  plotBy(Photo.mean~Hour.mean|T_treatment,data=leaf.m,legend=F,pch=16,xaxt="n",yaxt="n",type="b",cex=size,ylim=ylims,xlim=xlims,
         xlab="",ylab="",
         panel.first=adderrorbars(x=leaf.m$Hour.mean,y=leaf.m$Photo.mean,SE=leaf.m$Photo.standard.error,direction="updown"))
  plotBy(GPP_la_umol.mean~Hour|T_treatment,data=canopy.m,legend=F,pch=1,xaxt="n",yaxt="n",type="b",cex=size,ylim=ylims,xlim=xlims,add=T,
         panel.first=adderrorbars(x=canopy.m$Hour,y=canopy.m$GPP_la_umol.mean,SE=canopy.m$GPP_la_umol.standard.error,direction="updown"))
  magaxis(side=c(1,2,4),las=1)
  title(ylab=expression(Photo~(mu*mol~CO[2]~m^-2~s^-1)),
        xlab="Hour",cex.lab=2)
  legend("topright",pch=c(16,16,1,1),col=c("black","red","black","red"),legend=c("Leaf-A","Leaf-W","Canopy-A","Canopy-W"),cex=2)
  #-------------------------------------------------------------------------------------------------------------------
  if(export==T) dev.copy2pdf(file="Output/FigureS9.pdf")
  
}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------







#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#- plot the relationship between GPP, Ra, and CUE for two partitioning methods (Figure S4)
plotCUE_paritioning_method <- function(cue.day1=cue.day1,cue.day2=cue.day2,export=T){
  
  # windows(30,12);
  par(mfrow=c(1,3),mar=c(5,6,1,1),cex.lab=1.5,cex.axis=1.3,las=1)
  
  # GPP
  plot(cue.day2$GPP~cue.day1$GPP,xaxt="n",yaxt="n",xlab="",ylab="",pch="+");abline(0,1)
  magaxis(side=1:4,labels=c(1,1,0,0),frame.plot=T)
  title(xlab=expression(GPP*","~no~light~inhibition),
        ylab=expression(GPP*","~with~light~inhibition))
  lm.gpp <- lm(cue.day2$GPP~cue.day1$GPP);abline(lm.gpp,lty=2)
  legend("topleft",legend=paste("y = ",unname(round(coef(lm.gpp)[2],2)),"x",
                                unname(round(coef(lm.gpp)[1],2)),
                                sep=""),bty="n",cex=1.5)
  legend("bottomright","a",cex=1.5,bty="n",inset=-0.002)
  
  # Ra
  plot(cue.day2$Ra~cue.day1$Ra,xaxt="n",yaxt="n",xlab="",ylab="",pch="+");abline(0,1)
  magaxis(side=1:4,labels=c(1,1,0,0),frame.plot=T)
  title(xlab=expression(Ra*","~no~light~inhibition),
        ylab=expression(Ra*","~with~light~inhibition))
  lm.Ra <- lm(cue.day2$Ra~cue.day1$Ra);abline(lm.Ra,lty=2)
  legend("topleft",legend=paste("y = ",unname(round(coef(lm.Ra)[2],2)),"x","+",
                                unname(round(coef(lm.Ra)[1],2)),
                                sep=""),bty="n",cex=1.5)
  legend("bottomright","b",cex=1.5,bty="n",inset=-0.002)
  
  # CUE
  plot(cue.day2$RtoA~cue.day1$RtoA,xaxt="n",yaxt="n",xlab="",ylab="",pch="+");abline(0,1)
  magaxis(side=1:4,labels=c(1,1,0,0),frame.plot=T)
  title(xlab=expression(R[a]/GPP*","~no~light~inhibition),
        ylab=expression(R[a]/GPP*","~with~light~inhibition))
  y <- cue.day2$RtoA
  x <- cue.day1$RtoA
  lm.CUE <- lm(y~x+I(x^2))
  preds <- predict(lm.CUE,newdata=data.frame(x=seq(0,1,length.out=101)))
  lines(x=seq(0,1,length.out=101),y=preds,lty=2)
  legend("topleft",legend=paste("y = ",unname(round(coef(lm.CUE)[3],3)),expression(x^2),"+",
                                unname(round(coef(lm.CUE)[2],3)),"x","+",
                                unname(round(coef(lm.CUE)[1],3)),
                                sep=""),bty="n",cex=1.5)
  legend("bottomright","c",cex=1.5,bty="n",inset=-0.002)
  if(export==T) dev.copy2pdf(file="Output/FigureS4.pdf")
}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------






#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#- plot Ra/GPP vs. T for the two paritioning methods (with and without light inibition of R)
#- Creates Figure S5
plotCUEvsT_partitioning_method <- function(cue.day1=cue.day,cue.day2=cue.day2,export=T,shading=0.5){
  
  
  palette(c("black","red"))
  
  # windows(20,12);
  par(mfrow=c(1,2),cex.lab=1.5,mar=c(5,5,1,1),las=1)
  
  #- plot method 1; this should be with no light inhibition of leaf R
  smoothplot(Tair_24hrs, RtoA, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             pointcols=c(NA,NA),
             random="chamber",
             ylim=c(0,0.7),
             data=subset(cue.day1,PAR>35), kgam=5, axes=F,xlab="",ylab="")
  magaxis(side=1:4,labels=c(1,1,0,0))
  title(xlab=expression(T[air]~(degree*C)),
        ylab=expression(R[a]/GPP))
  legend("topright","a",bty="n",inset=-0.001,cex=1.2)
  
  #- plot method 2; this should be WITH light inhibition of leaf R
  smoothplot(Tair_24hrs, RtoA, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",pointcols=c(NA,NA),
             ylim=c(0,0.7),
             data=subset(cue.day2,PAR>35), kgam=5, axes=F,xlab="",ylab="")
  magaxis(side=1:4,labels=c(1,1,0,0))
  title(xlab=expression(T[air]~(degree*C)),
        ylab=expression(R[a]/GPP))
  legend("topright","b",bty="n",inset=-0.001,cex=1.2)
  if(export==T) dev.copy2pdf(file="Output/FigureS5.pdf")
}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#- plot the co-variance between air temperature and VPD (Figure S6)
plotVPD_Tair <- function(dat=dat.hr,export=F){
  
  dat$RH <- VPDtoRH(dat.hr$VPDair, dat.hr$Tair_al, Pa = 101)
  toplot <- subset(dat,RH>10 & PAR > 10) #- extract daytime values only, remove a few datapoints with questionable RH data
  toplot$DateTime <- as.factor(toplot$DateTime)
  toplot.hr <- dplyr::summarize(group_by(toplot,DateTime,T_treatment),
                                Tair_al=mean(Tair_al,na.rm=T),
                                VPDair=mean(VPDair,na.rm=T),
                                PAR=mean(PAR,na.rm=T))
  
  
  # windows(20,30);
  par(cex.axis=1.2,cex.lab=1.5,mar=c(5,7,1,1))
  plotBy(VPDair~Tair_al|T_treatment,col=c("black","red"),pch="+",data=toplot.hr,cex=0.5,xlim=c(2,47),ylim=c(0,8),
         xaxt="n",yaxt="n",
         legend=F,ylab=("VPD (kPa)"),xlab=expression(T[air]~(degree*C)))
  Ts <- seq(3,42,length=101)
  RH10 <- RHtoVPD(RH=10, TdegC=Ts, Pa = 101)
  RH20 <- RHtoVPD(RH=20, TdegC=Ts, Pa = 101)
  RH30 <- RHtoVPD(RH=30, TdegC=Ts, Pa = 101)
  RH40 <- RHtoVPD(RH=40, TdegC=Ts, Pa = 101)
  RH50 <- RHtoVPD(RH=50, TdegC=Ts, Pa = 101)
  RH60 <- RHtoVPD(RH=60, TdegC=Ts, Pa = 101)
  RH70 <- RHtoVPD(RH=70, TdegC=Ts, Pa = 101)
  RH80 <- RHtoVPD(RH=80, TdegC=Ts, Pa = 101)
  RH90 <- RHtoVPD(RH=90, TdegC=Ts, Pa = 101)
  lines(RH10~Ts,lty=3)
  lines(RH30~Ts,lty=3)
  lines(RH50~Ts,lty=3)
  lines(RH70~Ts,lty=3)
  lines(RH90~Ts,lty=3)
  xs <- rep(max(Ts),5)+3
  ys <- c(max(RH10),max(RH30),max(RH50),max(RH70),max(RH90))
  labels <- c("10%","30%","50%","70%","90%")
  textxy(X=xs,Y=ys,labs=labels,cex=0.9,offset=0)
  
  legend("topleft",pch=c("+","+",NA),lty=c(NA,NA,3),col=c("black","red","black"),legend=c("Ambient","Warmed","RH isolines"),bty="n")
  magaxis(side=c(1,2,3,4),labels=c(1,1,0,0),xlab="",ylab="",las=1)
  
  if(export==T)dev.copy2pdf(file="Output/FigureS6.pdf")
}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------





#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#- Does Ra correspond to GPP on the prior day? Not really.
Ra_GPP_coupling <- function(cue.day=cue.day,export=export,shading=0.7,kgam=3){
  cue.day <- subset(cue.day,Ranight<0)
  cue.day$Ranight_la <- with(cue.day,-1*Ranight/leafArea)
  cue.day$Cgain_la <- with(cue.day,Cgain/leafArea)

  # windows(20,20)
  palette(c("black","red"))
  par(mar=c(5,7,1,1))
  smoothplot(GPP_la,Ranight_la, T_treatment,polycolor=c(alpha("lightgrey",shading),alpha("lightgrey",shading)),
             random="chamber",xlab="",ylab="",
             xlim=c(0,9),
             ylim=c(0,1.5),
             data=cue.day, kgam=kgam, axes=FALSE)
  title(ylab=expression(Nightly~C~loss~(gC~m^-2~d^-1)),xlab=expression(GPP~(gC~m^-2~d^-1)),outer=F,cex.lab=1.5)

  magaxis(side=1:4,labels=c(1,1,0,0))
  if(export==T) dev.copy2pdf(file="output/FigureS10.pdf")
}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------




















#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#- get the raw tree size data
returnd2h <- function(plotson=0){
  
  
  require(HIEv)
  require(doBy)
  require(plotBy)
  require(zoo)

  # Function for converting all columns in a dataframe to numeric
  numericdfr <- function(dfr){
    for(i in 1:ncol(dfr)){
      options(warn=-1)
      oldna <- sum(is.na(dfr[,i]))
      num <- as.numeric(as.character(dfr[,i]))
      if(sum(is.na(num)) > oldna)
        next
      else
        dfr[,i] <- num
    }
    options(warn=0)
    return(dfr)
  }
  
  
  #------------------------------------------------------------------------------------------------------------
  # get the tree size data
  size <- read.csv("data/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV")
  size$chamber_n <- as.numeric(substr(size$chamber,start=2,stop=3))
  size$DateTime <- as.Date(size$DateTime)
  
  #do some processing to get d2h
  size2 <- base::subset(size,select=c("DateTime","Days_since_transplanting","chamber_n","chamber","T_treatment",
                                      "Water_treatment","Plant_height","Stem_number"))
  size2$diam <- size[,18]
  size2$diam_15 <- size[,16]
  size3 <- subset(size2,chamber_n<=12 & Stem_number==1)
  

  
  # tree 11 is different, because the height of stem 1 is not the actual true height. put the max height of the two stems in as the Plant_height
  tree11 <- subset(size2,chamber_n==11)
  tree11.max <- summaryBy(Plant_height~chamber_n+DateTime+Days_since_transplanting,data=tree11,FUN=max,na.rm=T,keep.names=T)
  size3[which(size3$chamber_n==11),"Plant_height"] <- tree11.max$Plant_height
  
  
  #- convert diameters to cm, heights to m
  size3$diam <- size3$diam/10
  size3$diam_15 <- size3$diam_15/10
  size3$Plant_height <- size3$Plant_height/100
  
  #get rid of NA's
  size4 <- subset(size3,diam>0)
  size4 <- size4[with(size4,order(DateTime,chamber_n)),]
  
  #get just the data with diameter at 15
  size_small <- subset(size3,diam_15>0 & year(DateTime)<2014)
  #------------------------------------------------------------------------------------------------------------
  
  
  # 
  # #------------------------------------------------------------------------------------------------------------
  # # gap fill diameter and height data
  # 
  # # create dataframe for all days
  # alldates <- rep(seq.Date(from=as.Date(range(size4$DateTime)[1]),to=as.Date(range(size4$DateTime)[2]),by="day"),12)
  # chamber_n <- c(rep(1,length(alldates)/12),rep(2,length(alldates)/12),rep(3,length(alldates)/12),rep(4,length(alldates)/12),
  #                rep(5,length(alldates)/12),rep(6,length(alldates)/12),rep(7,length(alldates)/12),rep(8,length(alldates)/12),
  #                rep(9,length(alldates)/12),rep(10,length(alldates)/12),rep(11,length(alldates)/12),rep(12,length(alldates)/12))
  # datedf <- data.frame(chamber_n=chamber_n,DateTime=alldates)                                                                                             
  # 
  # #merge data in with dataframe of all days
  # size5 <-merge(size4,datedf,all=T,by=c("chamber_n","DateTime"))
  # 
  # #break across list, gapfill list
  # size6 <- zoo(size5)
  # size6$Days_since_transplanting <- na.approx(size6$Days_since_transplanting)
  # size6$Plant_height <- na.approx(size6$Plant_height)
  # size6$diam <- na.approx(size6$diam)
  # 
  # 
  # # get it back to a normal dataframe
  # size7 <- numericdfr(fortify.zoo(size6))
  # size7$DateTime <- as.Date(size7$DateTime)
  # size7$d2h <- with(size7,(diam/10)^2*Plant_height)
  # 
  # # put some of the other bits back together
  # size7$chamber <- as.factor(paste0("C",sprintf("%02.0f",size7$chamber_n)))
  # size7$T_treatment <- as.factor(ifelse(size7$chamber_n %% 2 ==1,"ambient","elevated"))
  # size7$Index <- NULL
  # size7$Stem_number <- NULL
  # #plotBy(d2h~DateTime|chamber_n,data=size7)
  # 
  
  if (plotson==1){
    # windows(12,12);
    par(mfrow=c(2,1))
    plotBy(diam~DateTime|chamber_n,data=size4,pch=15,ylim=c(0,100),legend=F,ylab="Diameter (mm)")
    plotBy(diam~DateTime|chamber,type="l",lwd=2,data=size7,add=T,legend=F)
    plotBy(Plant_height~DateTime|chamber_n,data=size4,pch=15,ylim=c(0,1000),legend=F,ylab="Height (cm)")
    plotBy(Plant_height~DateTime|chamber,type="l",lwd=2,data=size7,add=T,legend=F)
  }
  
  return(list(size4,size_small))
}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#- plot tree size over time

plot_tree_size <- function(export=export){
  
  #-------------------
  #- Get the three direct observations of leaf area. They happened on 9 Sept 2013, 10 Feb 2014, and
  #    during the harvest at ~25 May 2014.
  treeMass <- read.csv("data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
  treeMass.sub <- subset(treeMass,as.Date(DateTime) %in% as.Date(c("2013-09-14","2014-02-10","2014-05-27")))
  treeMass.sub$Date <- as.Date(treeMass.sub$DateTime)
  leafArea <- summaryBy(leafArea~Date+chamber+T_treatment,data=subset(treeMass.sub,Water_treatment=="control"),FUN=c(mean),keep.names=T)
  leafArea1 <- summaryBy(leafArea~Date+T_treatment,data=leafArea,FUN=c(mean,standard.error))
  
  #-------------------
  
  #-------------------
  #- Get tree diameters and heights
  size.list <- returnd2h(plotson=0)
  size <- size.list[[1]] # extract just the data after transplanting
  size.m <- summaryBy(diam+Plant_height~DateTime+T_treatment,data=subset(size,Water_treatment=="control"),FUN=c(mean,standard.error))
  
  #- process the diameter at 15cm data, possibly for an inset
  size2 <- subset(size.list[[2]])
  size2.m <- summaryBy(diam_15+Plant_height~T_treatment,data=subset(size2,DateTime==as.Date("2012-12-12")),FUN=c(mean,standard.error))
  
  
  #-------------------
  # windows(20,40);
  par(mfrow=c(3,1),mar=c(0,9,0,3),oma=c(6,0,2,0),las=1,cex.axis=1.5)
  
  #- plot diameter
  xlims <-as.Date(c("2013-3-1","2014-6-1"))
  yearlims <-as.Date(c("2013-1-1","2014-1-1"))
  plotBy(diam.mean~DateTime|T_treatment,data=size.m,pch=16,type="o",ylim=c(0,8),
         xlim=xlims,
         xaxt="n",yaxt="n",xlab="",ylab="",legend=F,
         panel.first=adderrorbars(x=size.m$Date,y=size.m$diam.mean,barlen=0.02,
                                  SE=size.m$diam.standard.error,direction="updown",
                                  col=c("black","red")))
  magaxis(side=c(2,4),labels=c(1,1),frame.plot=T,majorn=3,las=1)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="month"),tcl=0.25,labels=F)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="quarter"),tcl=0.75,labels=F)
  title(ylab=expression(Diameter~(cm)),cex.lab=2)
  abline(v=as.Date("2013-9-13"),lty=2)
  text(x=as.Date("2013-9-13"),y=8.8,labels="Floors sealed",xpd=NA,cex=1.3)
  legend("topleft",pch=16,lty=c(1),col=c("black","red"),legend=c("Ambient","Warmed"),bty="n",cex=1.2)
  legend("bottomright","a",bty="n",inset=0.002,cex=1.2)
  
  #- plot height
  plotBy(Plant_height.mean~DateTime|T_treatment,data=size.m,pch=16,type="o",ylim=c(0,11),
         xlim=xlims,
         xaxt="n",yaxt="n",xlab="",ylab="",legend=F,
         panel.first=adderrorbars(x=size.m$Date,y=size.m$Plant_height.mean,barlen=0.02,
                                  SE=size.m$Plant_height.standard.error,direction="updown",
                                  col=c("black","red")))
  magaxis(side=c(2,4),labels=c(1,1),frame.plot=T,majorn=3,las=1)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="month"),tcl=0.25,labels=F)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="quarter"),tcl=0.75,labels=F)
  title(ylab=expression(Height~(m)),cex.lab=2)
  abline(v=as.Date("2013-9-13"),lty=2)
  legend("bottomright","b",bty="n",inset=0.002,cex=1.2)
  
  #-- plot leaf area over time
  plotBy(leafArea.mean~Date|T_treatment,data=leafArea1,pch=16,type="p",ylim=c(0,25),cex=1.5,
         xlim=xlims,
         xaxt="n",yaxt="n",xlab="",ylab="",legend=F,
         panel.first=adderrorbars(x=leafArea1$Date,y=leafArea1$leafArea.mean,barlen=0.02,
                                  SE=leafArea1$leafArea.standard.error,direction="updown",
                                  col=c("black","red","black","red","black","red")))
  magaxis(side=c(2,4),labels=c(1,1),frame.plot=T,majorn=3,las=1)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="month"),tcl=0.25,labels=F)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="quarter"),tcl=0.75,labels=T,las=2,
            format="%b")
  title(ylab=expression(Total~leaf~area~(m^2)),cex.lab=2)
  abline(v=as.Date("2013-9-13"),lty=2)
  legend("bottomright","c",bty="n",inset=0.002,cex=1.2)
  
  if(export==T) dev.copy2pdf(file="output/treeSize.pdf")
  
}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------





#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#- function to create a table of tree size and leaf area at the beginning and end of the experiment
table_tree_size <- function(){
  
  #-------------------
  #- Get the three direct observations of leaf area. They happened on 9 Sept 2013, 10 Feb 2014, and
  #    during the harvest at ~25 May 2014.
  treeMass <- read.csv("data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
  treeMass.sub <- subset(treeMass,as.Date(DateTime) %in% as.Date(c("2013-09-14","2014-02-10","2014-05-27")))
  treeMass.sub$Date <- as.Date(treeMass.sub$DateTime)
  leafArea <- summaryBy(leafArea~Date+chamber+T_treatment,
                        data=subset(treeMass.sub,Water_treatment=="control" & Date %in% as.Date(c("2013-09-14","2014-05-27"))),
                        FUN=c(mean),keep.names=T)
  leafArea1 <- summaryBy(leafArea~Date+T_treatment,data=leafArea,FUN=c(mean,standard.error))
  
  leafArea$Date <- as.factor(leafArea$Date)
  #-------------------
  
  #-------------------
  #- Get tree diameters and heights
  size.list <- returnd2h(plotson=0)
  size <- size.list[[1]] # extract just the data after transplanting
  size$Date <- as.Date(size$DateTime)
  size.m <- summaryBy(diam+Plant_height~Date+T_treatment,
                      data=subset(size,Water_treatment=="control"& Date %in% as.Date(c("2013-09-05","2014-05-27"))),FUN=c(mean,standard.error))
  
  
  
  
  #-------------------
  # make a table from leafArea1 and size.m
  table1 <- data.frame(Date=leafArea1$Date,Comment=c(rep("Floors sealed",2),rep("Harvest",2)),
                        Treatment=leafArea1$T_treatment)
  table1$Diameter <- paste(sprintf("%.1f",round(size.m$diam.mean,2))," (",sprintf("%.1f",round(size.m$diam.standard.error,1)),")",sep="")
  table1$Height <- paste(sprintf("%.1f",round(size.m$Plant_height.mean,2))," (",sprintf("%.1f",round(size.m$Plant_height.standard.error,1)),")",sep="")
  table1$Canopy <- paste(sprintf("%.1f",round(leafArea1$leafArea.mean,2))," (",sprintf("%.1f",round(leafArea1$leafArea.standard.error,1)),")",sep="")
  
  names(table1) <- c("Date","Comment","Treatment","Diameter (cm)","Height (m)","Total leaf area (m^2)")
  table1 
  
}



