#----------------------------------------------------------------------------------------------------------------
#- This is a collection of many functions that do the actual work of data analysis and plotting. These functions
#    are called by just a few lines of code in "CentralScript.R" to recreate the analyses and figures.
#----------------------------------------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------------------------------------------
#' ## script for converting .Rmd files to .R scripts.

#' #### Kevin Keenan 2014
#' 
#' This function will read a standard R markdown source file and convert it to 
#' an R script to allow the code to be run using the "source" function.
#' 
#' The function is quite simplisting in that it reads a .Rmd file and adds 
#' comments to non-r code sections, while leaving R code without comments
#' so that the interpreter can run the commands.
#' 
#' 
rmd2rscript <- function(infile){
  # read the file
  flIn <- readLines(infile)
  # identify the start of code blocks
  cdStrt <- which(grepl(flIn, pattern = "```{r*", perl = TRUE))
  # identify the end of code blocks
  cdEnd <- sapply(cdStrt, function(x){
    preidx <- which(grepl(flIn[-(1:x)], pattern = "```", perl = TRUE))[1]
    return(preidx + x)
  })
  # define an expansion function
  # strip code block indacators
  flIn[c(cdStrt, cdEnd)] <- ""
  expFun <- function(strt, End){
    strt <- strt+1
    End <- End-1
    return(strt:End)
  }
  idx <- unlist(mapply(FUN = expFun, strt = cdStrt, End = cdEnd, 
                       SIMPLIFY = FALSE))
  # add comments to all lines except code blocks
  comIdx <- 1:length(flIn)
  comIdx <- comIdx[-idx]
  for(i in comIdx){
    flIn[i] <- paste("#' ", flIn[i], sep = "")
  }
  # create an output file
  nm <- strsplit(infile, split = "\\.")[[1]][1]
  flOut <- file(paste(nm, ".R", sep = ""), "w")
  for(i in 1:length(flIn)){
    cat(flIn[i], "\n", file = flOut, sep = "\t")
  }
  close(flOut)
}

#----------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#- function to add regression equation to statistical smoothing 
stat_smooth_func <- function(mapping = NULL, data = NULL,
                             geom = "smooth", position = "identity",
                             ...,
                             method = "auto",
                             formula = y ~ x,
                             se = TRUE,
                             n = 80,
                             span = 0.75,
                             fullrange = FALSE,
                             level = 0.95,
                             method.args = list(),
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE,
                             xpos = NULL,
                             ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}
StatSmoothFunc <- ggproto("StatSmooth", Stat,
                          
                          setup_params = function(data, params) {
                            # Figure out what type of smoothing to do: loess for small datasets,
                            # gam with a cubic regression basis for large data
                            # This is based on the size of the _largest_ group.
                            if (identical(params$method, "auto")) {
                              max_group <- max(table(data$group))
                              
                              if (max_group < 1000) {
                                params$method <- "loess"
                              } else {
                                params$method <- "gam"
                                params$formula <- y ~ s(x, bs = "cs")
                              }
                            }
                            if (identical(params$method, "gam")) {
                              params$method <- mgcv::gam
                            }
                            
                            params
                          },
                          
                          compute_group = function(data, scales, method = "auto", formula = y~x,
                                                   se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                   xseq = NULL, level = 0.95, method.args = list(),
                                                   na.rm = FALSE, xpos=NULL, ypos=NULL) {
                            if (length(unique(data$x)) < 2) {
                              # Not enough data to perform fit
                              return(data.frame())
                            }
                            
                            if (is.null(data$weight)) data$weight <- 1
                            
                            if (is.null(xseq)) {
                              if (is.integer(data$x)) {
                                if (fullrange) {
                                  xseq <- scales$x$dimension()
                                } else {
                                  xseq <- sort(unique(data$x))
                                }
                              } else {
                                if (fullrange) {
                                  range <- scales$x$dimension()
                                } else {
                                  range <- range(data$x, na.rm = TRUE)
                                }
                                xseq <- seq(range[1], range[2], length.out = n)
                              }
                            }
                            # Special case span because it's the most commonly used model argument
                            if (identical(method, "loess")) {
                              method.args$span <- span
                            }
                            
                            if (is.character(method)) method <- match.fun(method)
                            
                            base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                            model <- do.call(method, c(base.args, method.args))
                            
                            m = model
                            eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                             list(a = format(coef(m)[1], digits = 3), 
                                                  b = format(coef(m)[2], digits = 3), 
                                                  r2 = format(summary(m)$r.squared, digits = 3)))
                            func_string = as.character(as.expression(eq))
                            
                            if(is.null(xpos)) xpos = min(data$x)*0.9
                            if(is.null(ypos)) ypos = max(data$y)*1.1
                            data.frame(x=xpos, y=ypos, label=func_string)
                            
                          },
                          
                          required_aes = c("x", "y")
)
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#- function to partition the hourly net flux observations into GPP and Ra using an Arrhenious function
partitionHourlyFluxCUE_arr <- function(data.hr.gf=data.hr.gf,Ea=57.69,lagdates,leafRtoTotal = 1,leafRreduction=0){
  rvalue = 8.134
  require(data.table)
  
  #-- convert mmolCO2 s-1 to gC hr-1
  data.hr.gf$FluxCO2_g <- with(data.hr.gf,FluxCO2*60*60/1000*12.0107)
  data.hr.gf$period <- ifelse(data.hr.gf$PAR>2,"Day","Night")
  
  #-- partition day-time net C exchange into GPP and Ra, similar to how it is done in eddy-covariance.
  #-- create a series of dates
  date.vec <- seq.Date(from=min(data.hr.gf$Date),to=max(data.hr.gf$Date),by="day")
  
  #-- estimate R-Tref and Tref for each date for each chamber
  #lagDates <- 3 # establish how many prior days to include
  RTdat <- expand.grid(Date=date.vec,chamber=levels(data.hr.gf$chamber))
  RTdat$Tref <- RTdat$R_Tref <- NA
  
  
  # print("Partitioning Net CO2 fluxes into GPP and Ra")
  # #- set up progress bar to track that this is working
  # pb <- txtProgressBar(min = 0, max = nrow(RTdat), style = 3)
  
  for (i in 1:nrow(RTdat)){
    #- trial a data.table alternative to speed this up. The filter thing was actually slower.
    
    #dat <- dplyr::filter(data.hr.gf,chamber==RTdat$chamber[i],Date <= RTdat$Date[i], Date >= (RTdat$Date[i]-lagDates),period =="Night")
    #RTdat$Tref[i] <- mean(dat$Tair_al,na.rm=T)
    #RTdat$R_Tref[i] <- mean(dat$FluxCO2_g,na.rm=T)
    
    inds <- which(data.hr.gf$chamber==RTdat$chamber[i] & data.hr.gf$Date <= RTdat$Date[i] & data.hr.gf$Date >= (RTdat$Date[i]-lagdates) & data.hr.gf$period =="Night" )
    RTdat$Tref[i] <- mean(data.hr.gf$Tair_al[inds],na.rm=T)
    RTdat$R_Tref[i] <- mean(data.hr.gf$FluxCO2_g[inds],na.rm=T)
    # setTxtProgressBar(pb, i)
    
  }
  # close(pb)
  
  RTdat$Tref_K <- with(RTdat,Tref+273.15)
  
  #-- merge these reference data into the gap-filled flux dataframe to estimate Ra during the daytime, and hence GPP
  data.hr.gf3 <- merge(data.hr.gf,RTdat,by=c("Date","chamber"))
  #data.hr.gf3$Ra_est <- with(data.hr.gf3,R_Tref*Q10^((Tair_al-Tref)/10)) # estimate respiration rate. This is a negative number.
  data.hr.gf3$Ra_est <- with(data.hr.gf3,R_Tref*exp((Ea*1000/(rvalue*Tref_K))*(1-Tref_K/(Tair_al+273.15)))) # estimate respiration rate. This is a negative number.
  data.hr.gf3$Ra_est <- ifelse(data.hr.gf3$period=="Day",
                              data.hr.gf3$Ra_est-leafRreduction*leafRtoTotal*data.hr.gf3$Ra_est,
                              data.hr.gf3$Ra_est) # estimate respiration rate. This is a negative number. If it's day, subtract 30% from the leaf R fraction
  

  data.hr.gf3$GPP <- ifelse(data.hr.gf3$period=="Night",0,data.hr.gf3$FluxCO2_g-data.hr.gf3$Ra_est)
  data.hr.gf3$Ra <- ifelse(data.hr.gf3$period=="Night",data.hr.gf3$FluxCO2_g,data.hr.gf3$Ra_est)
  
  return(data.hr.gf3)
}
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
#- function to return daily sums of GPP and Ra from hourly data
returnCUE.day <- function(dat=data.hr.p){
  
  # #- calculate leaf-area specific rates of GPP, convert to umol CO2 m-2 s-1
  dat$GPP_la <- with(dat,GPP/leafArea)
  dat$GPP_la_umol <- with(dat,GPP_la/12*1*10^6/60/60)
  dat$PAR_mol <- dat$PAR*60*60*1*10^-6

  #- create a date variable that moves the window for "night" observations, to sum observations for the night period following a day period
  dat$Date2 <- as.Date(dat$DateTime)
  dat$hour <- hour(dat$DateTime)
  earlynights <- which(dat$hour < 12 & dat$period =="Night")
  dat$Date2[earlynights] <- dat$Date2[earlynights]-1
  
  #- Identify and name the drought/watered treatments
  drought.chamb = unique(dat$chamber[ dat$Water_treatment %in% as.factor("drydown")])
  dat$chamber_type = as.factor( ifelse(dat$chamber %in% drought.chamb, "drought", "watered") )
  
  #- create daily sums 
  data.day <- dplyr::summarize(group_by(dat,Date2,chamber,T_treatment,chamber_type,period),
                              GPP = sum(GPP,na.rm=T),
                              Ra = sum(Ra,na.rm=T),
                              FluxCO2_g=sum(FluxCO2_g,na.rm=T),
                              Tair_al=mean(Tair_al,na.rm=T),
                              VPDair=max(VPDair,na.rm=T),
                              PAR=sum(PAR_mol,na.rm=T),
                              leafArea=mean(leafArea,na.rm=T))
  data.day <- as.data.frame(data.day)
  data.day <- data.day[with(data.day,order(Date2,chamber)),]
  names(data.day)[1] <- "Date"
  
  
  Tair_day <- dplyr::summarize(group_by(dat,Date2,chamber,T_treatment,chamber_type),
                               Tair_al=mean(Tair_al,na.rm=T))
  names(Tair_day)[5] <- "Tair_24hrs"                   
  
  data.day_tair <- as.data.frame(Tair_day)
  data.day <- as.data.frame(data.day)
  
  #- merge night and day data
  data.day2d <- subset(data.day,period=="Day")[,c("Date","chamber","T_treatment","GPP","Ra","FluxCO2_g","PAR","Tair_al","VPDair","leafArea","chamber_type")]
  names(data.day2d)[4:8] <- c("GPP","Raday","Cgain","PAR","T_day")
  data.day2n <- subset(data.day,period=='Night')[,c("Date","chamber","T_treatment","Ra","FluxCO2_g","Tair_al","VPDair","chamber_type")]
  names(data.day2n)[4:7] <- c("Ranight","Closs","T_night","VPDair_night")
  names(data.day2n)[1] <- c("Date")
  
  #- merge data such that the nightly data are one day ahead of the daily data (i.e., does tonight's respiration depend on today's photosynthesis?)
  cue.day1 <- subset(merge(data.day2d,data.day2n,by=c("Date","chamber","T_treatment","chamber_type")))
  cue.day <- merge(cue.day1,data.day_tair,by.x=c("Date","chamber","T_treatment","chamber_type"),by.y=c("Date2","chamber","T_treatment","chamber_type"))
  cue.day$Ra <- with(cue.day,-1*(Raday+Ranight))
  cue.day$RtoA <- with(cue.day,Ra/GPP)
  cue.day$CUE <- with(cue.day,(1-(Ra/GPP)))
  cue.day$GPP_la <- with(cue.day,GPP/leafArea)
  cue.day$Ra_la <- with(cue.day,Ra/(leafArea)) #g C m-2 day-1
  
  # remove a few days from the beginning of the dataset in C09. Problem with flux data.
  toremove <- which(cue.day$Ra_la < 1.5 & cue.day$chamber=="C07" & cue.day$Date < as.Date("2013-10-3"))
  cue.day <- cue.day[-toremove,] 
  
  # #- merge in the key
  # key <- read.csv("raw_data/WTC_TEMP_CM_TREATKEY_20121212-20140528_L1_v1.csv")
  # key$chamber_type = as.factor( ifelse(key$chamber %in% drought.chamb, "drought", "watered") )
  # 
  # cue.day$T_treatment <- ifelse(cue.day$T_treatment=="ambient","ambient","elevated")
  # cue.day2 <- merge(cue.day,key,by=c("chamber","Date","T_treatment","chamber_type"))
  
  # cue.day.trt <- summaryBy(GPP+GPP_la+Ra+Ra_la+RtoA+PAR+T_day+T_night+Tair_24hrs+VPDair+VPDair_night~Date+T_treatment+Water_treatment,data=cue.day2,FUN=c(mean,standard.error))
  # cue.day.trt <- summaryBy(GPP+Ra+leafArea~Date+T_treatment+chamber_type,data=cue.day2,FUN=c(mean,standard.error,length))
  # return(list(cue.day2,cue.day.trt))
  cue.day.trt <- summaryBy(GPP+Ra+leafArea~Date+T_treatment,data=cue.day,FUN=c(mean,standard.error))
  return(list(cue.day,cue.day.trt))
}
#----------------------------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#-- Code to return tree mass for every day of the experiment.
#- this code from WTC_Rallometry_script.R in the WTC3 folder. I'm trying to estimate tree mass at the beginning of the experiment
returnTreeMass <- function(plotson=F){
  
  #- get an estimate of d2h for every day of the experiment by interpolation
  size <- returnd2h()
  
  #- get an estimate of volume and mass for wood and bark for each measurement day
  # print("Calculating stem volume")
  vol <- getvol()
  
  #- interpolate volume for every day of the experiment
  # print("Interpolating stem volume")
  vol.all <- gapfillvol(vol)
  
  
  #- get branch mass
  branchMass <- interpolateBranchDW(vol.all=vol.all,plotson=0) 
  stem_branch_mass <- merge(vol.all,branchMass,by=c("chamber","Date"))
  
  
  
  #- now get an estimate of leaf mass for every day of the experiment
  # get the interpolated leaf areas (laDaily)
  
  dat.hr <- data.frame(data.table::fread("raw_data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv"))
  dat.hr$DateTime <- as.POSIXct(dat.hr$DateTime,format="%Y-%m-%d %H:%M:%S",tz="GMT")
  dat.hr$Date <- as.Date(dat.hr$DateTime)
  dat.hr$chamber <- as.factor(dat.hr$chamber)
  laDaily <- summaryBy(leafArea~Date+chamber,data=dat.hr,keep.names=T,na.rm=T)
  
  # get a canopy-weighted SLA estimate
  #setwd("//ad.uws.edu.au/dfshare/HomesHWK$/30035219/My Documents/Work/HFE/WTC3/gx_wtc3/")
  harvest <- read.csv("raw_data/WTC_TEMP_CM_HARVEST-CANOPY_20140526-20140528_L1_v1.csv")
  leafmass <- summaryBy(TotalLeafDM~chamber,data=harvest,FUN=sum,keep.names=F)
  harvest2 <- merge(harvest,leafmass,by="chamber")
  harvest2$weights <- with(harvest2,TotalLeafDM/TotalLeafDM.sum)
  SLA <- plyr::ddply(harvest2,~chamber,summarise, SLA=weighted.mean(SLA,weights))
  
  # merge SLA and leaf area to estimate leaf mass for each day
  leafMass <- merge(laDaily,SLA,by=c("chamber"))
  leafMass$leafMass <- with(leafMass,leafArea/SLA*10000)
  
  treeMass1 <- merge(stem_branch_mass,leafMass,by=c("chamber","Date"))
  treeMass1$boleMass <- with(treeMass1,mass_wood+mass_bark)
  treeMass1$totMass <- rowSums(treeMass1[,c("boleMass","branchMass","leafMass")])
  treeMass <-subset(treeMass1,select=c("chamber","T_treatment","Water_treatment","Date","Measurement","Days_since_transplanting","branchMass","boleMass","leafMass","totMass","SLA"))
  
  # print("Done. Returned treeMass")
  return(treeMass)
}
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#- Read and process the diameter and height dataset, return a dataframe of gapfilled size data.
returnd2h <- function(){
  
  #------------------------------------------------------------------------------------------------------------
  # get the tree size data
  size <- read.csv("raw_data/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV")
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
  
  
  
  #------------------------------------------------------------------------------------------------------------
  # gap fill diameter and height data
  
  # create dataframe for all days
  alldates <- rep(seq.Date(from=as.Date(range(size4$DateTime)[1]),to=as.Date(range(size4$DateTime)[2]),by="day"),12)
  chamber_n <- c(rep(1,length(alldates)/12),rep(2,length(alldates)/12),rep(3,length(alldates)/12),rep(4,length(alldates)/12),
                 rep(5,length(alldates)/12),rep(6,length(alldates)/12),rep(7,length(alldates)/12),rep(8,length(alldates)/12),
                 rep(9,length(alldates)/12),rep(10,length(alldates)/12),rep(11,length(alldates)/12),rep(12,length(alldates)/12))
  datedf <- data.frame(chamber_n=chamber_n,DateTime=alldates)                                                                                             
  
  #merge data in with dataframe of all days
  size5 <-merge(size4,datedf,all=T,by=c("chamber_n","DateTime"))
  
  #break across list, gapfill list
  size6 <- zoo(size5)
  size6$Days_since_transplanting <- na.approx(size6$Days_since_transplanting)
  size6$Plant_height <- na.approx(size6$Plant_height)
  size6$diam <- na.approx(size6$diam)
  
  
  # get it back to a normal dataframe
  size7 <- numericdfr(fortify.zoo(size6))
  size7$Date <- as.Date(size7$DateTime)
  size7$d2h <- with(size7,(diam/10)^2*Plant_height)
  
  # put some of the other bits back together
  size7$chamber <- as.factor(paste0("C",sprintf("%02.0f",size7$chamber_n)))
  size7$T_treatment <- as.factor(ifelse(size7$chamber_n %% 2 ==1,"ambient","elevated"))
  size7$Water_treatment <- "control"
  size7$Index <- NULL
  size7$Stem_number <- NULL
  
  #- establish a treatment key for the Water_treatment variable
  key <- data.frame(chamber=levels(size7$chamber),
                    Water_treatment=c("drydown","control","drydown","drydown","control","drydown",
                                      "control","drydown","control","control","drydown","control"))
  
  size_before <- subset(size7,Date<as.Date("2014-02-4"))
  size_after <- subset(size7,Date>=as.Date("2014-02-4"))
  size_after$Water_treatment <- NULL
  size_after2 <- merge(size_after,key,by="chamber")
  
  #- combined dataframes from before and after the drought began
  size8 <- rbind(size_before,size_after2)
  
  #- clean up dataframe for output
  size_out <- size8[,c("Date","chamber","T_treatment","Water_treatment","diam","Plant_height","d2h")]
  
  
  return(size_out)
}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
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

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
# get an estimate of bole volume for every measurement day
getvol <- function(){

  # browser()
  # get the data. the "interpolated 38" csv was manual editied in excel to linearly interpolate the UPPER stem estimates
  #   on measurement date 38 as the mean of measurements #37 and #39.
  #downloadCSV(filename="Wdata/from HIEv/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV,topath="raw_data/from HIEv/")
  #size <- read.csv("raw_data/from HIEv/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV")
  size <- read.csv("raw_data/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2_interpolated38.csv")
  names(size)[15:47] <- paste("d",substr(names(size)[15:47],start=2,stop=5),sep="")
  size <- size[,1:49]
  size$DateTime <- as.Date(size$DateTime)
  
  #ignore the reference plots
  size2 <- subset(size,T_treatment !="reference")
  size2$T_treatment <- factor(size2$T_treatment)
  size2$Water_treatment <- factor(size2$Water_treatment)
  size2$chamber <- factor(size2$chamber)
  
  #melt into "long" format
  sLong <- reshape::melt(size2,measure.vars=names(size)[15:47],value.name="diam",variable.name="height")
  names(sLong)[(ncol(sLong)-1):ncol(sLong)] <- c("height","diam")
  sLong$height_n <- as.numeric(substr(sLong$height,start=2,stop=4))
  
  sLong <- subset(sLong,height_n>=65)
  
  
  
  # estimate volume
  # split long dataframe into a list for each chamber and each measurement date
  sLong$chamber_date <- paste(sLong$chamber,sLong$Measurement,sep="_")
  sLong.l <- split(sLong,sLong$chamber_date)
  
  #-- add a number for the ground (assume no taper) and a value for the maximum tree height (assume diameter of 0.1mm)
  for (i in 1:length(sLong.l)){
    # add a line to the dataframe for the diameter at floor height
    firstline <- sLong.l[[i]][1,]  
    firstline$height_n[1] <- 0 #. Edited to give Mike an estimate of total tree volume to the ground.
    
    # add a line to the dataframe for the tiny diameter at total plant height
    lastline <- sLong.l[[i]][nrow(sLong.l[[i]]),]
    lastline$height_n[1] <- lastline$Plant_height[1]
    lastline$diam[1] <- 0.1 
    sLong.l[[i]] <- rbind(firstline,sLong.l[[i]],lastline)
    sLong.l[[i]] <- sLong.l[[i]][which(is.na(sLong.l[[i]]$diam)==F),] # get rid of NA's
  }
  treevol(sLong.l[[100]]) # example of volume calculation for a single observation of a single tree
  vols <- do.call(rbind,lapply(sLong.l,FUN=treevol))
  return(vols)
}

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
# estimate tree wood and bark volume and wood and bark mass from a series of measurements of diameter, and the height of that diameter measurement
treevol <- function(dat){
  dat2 <- subset(dat,diam>0)
  
  #-- get the bark depth estimates (along with wood and bark densities)
  density <- getBarkWoodDensity()
  density$barkDepth <- with(density,(diamoverbark-diamunderbark)/2)
  density$log_barkDepth <- log10(density$barkDepth)
  density$log_diamoverbark <- log10(density$diamoverbark)
  lm_bd <- lm(log_barkDepth~log_diamoverbark,data=density) # fit the harvest data
  dens.coefs <- unname(coef(lm_bd))
  dat2$bark_depth <- 10^(dens.coefs[1]+dens.coefs[2]*log10(dat2$diam)) # apply bark depths to estimate inner diameters
  dat2$diam_inner <- with(dat2,diam-bark_depth)  
  
  # estimate density from diameter
  bd.coef <- unname(coef(lm(bark_density~diamoverbark,data=density)))
  wd.coef <- unname(coef(lm(wooddensity~diamoverbark,data=density)))
  
  
  #loop over each stem segment and calculate volume as the fustrum of a cone. Assume the (unmeasured) bottom bit is a cylinder.
  vol <- 0
  vol_inner <- 0
  vol.cum <- 0
  vol.cum.inner <- 0
  mass_wood <- vol_wood <- 0
  mass_bark <- vol_bark <- 0
  for (i in 1:nrow(dat2)){
    # get the vertical segment length
    if(i==1){h1 <- 0}
    if(i > 1){h1 <- dat2[i-1,"height_n"]}
    h2 <- dat2[i,"height_n"]
    h <- h2-h1
    
    # get the radii in cm
    if(i==1){r1 <- dat2[i,"diam"]/20}
    if(i > 1){r1 <- dat2[i-1,"diam"]/20}
    r2 <- dat2[i,"diam"]/20
    
    #outer volume
    vol[i] <- pi*h/3* (r1^2+r1*r2+r2^2) # volume in cm3, fustrum of a cone. See http://jwilson.coe.uga.edu/emt725/Frustum/Frustum.cone.html
    
    
    # get the inner radii in cm
    
    if(i==1){r1.inner <- dat2[i,"diam_inner"]/20}
    if(i > 1){r1.inner <- dat2[i-1,"diam_inner"]/20}
    r2.inner <- dat2[i,"diam_inner"]/20
    
    vol_inner[i] <- pi*h/3* (r1.inner^2+r1.inner*r2.inner+r2.inner^2) # volume in cm3, fustrum of a cone. See http://jwilson.coe.uga.edu/emt725/Frustum/Frustum.cone.html
    
    # get wood volume and mass
    vol_wood[i] <- vol_inner[i]
    mass_wood[i] <- vol_wood[i]*(wd.coef[1]+wd.coef[2]*mean(r1.inner,r2.inner)) #wood mass is volume time density
    
    # get bark volume and mass
    vol_bark[i] <- vol[i] - vol_inner[i]
    mass_bark[i] <- vol_bark[i]*(bd.coef[1]+bd.coef[2]*mean(r1.inner,r2.inner)) #wood mass is volume time density
  }
  
  
  df.out <- dat2[1,c(1,2,3,4,5,6)]
  df.out$diam <- max(dat2$diam)
  df.out$height <- max(dat2$Plant_height)
  df.out$vol <- sum(vol,na.rm=T)
  df.out$vol_wood <- sum(vol_wood,na.rm=T)
  df.out$vol_bark <- sum(vol_bark,na.rm=T)
  df.out$mass_wood <- sum(mass_wood,na.rm=T)
  df.out$mass_bark <- sum(mass_bark,na.rm=T)
  
  return(df.out)
}

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------
#- function to return the measured stemwood and bark densities
getBarkWoodDensity <- function(){
  
  #--- Bark density data
  bark <- read.csv("raw_data/WTC_TEMP_CM_BARKDENSITY_20140528_L1.csv")
  
  
  #--- wood density data
  wood <- read.csv("raw_data/WTC_TEMP_CM_WOODDENSITY_20140528_L1.csv")
  
  wood$chamber_n <- as.numeric(substr(wood$chamber,start=2,stop=3))
  wood$T_treatment <- as.factor(ifelse(wood$chamber_n %% 2 == 1, "ambient","elevated"))
  wood2 <- subset(wood,Layer != "G Tip")
  wood2$Layer <- factor(wood2$Layer)
  wood2$position <- as.factor(ifelse(wood2$Layer == "Base","low",
                                     ifelse(wood2$Layer == "Middle","mid","top")))
  
  #-- merge wood and bark data
  dense <- merge(wood2,bark,by=c("chamber","position","T_treatment"))
  
  dense_out <- dense[,c("chamber","position","T_treatment","wooddensity","bark_density","diamoverbark","diamunderbark")]
  return(dense_out)
}

#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
# gapfill stem volume (for daily data)
gapfillvol <- function(dat){
  
  
  
  dat$chamber_n <- as.numeric(substr(dat$chamber,start=2,stop=3))
  
  # create dataframe for all days
  alldates <- rep(seq.Date(from=as.Date(range(dat$DateTime)[1]),to=as.Date(range(dat$DateTime)[2]),by="day"),12)
  chamber_n <- c(rep(1,length(alldates)/12),rep(2,length(alldates)/12),rep(3,length(alldates)/12),rep(4,length(alldates)/12),
                 rep(5,length(alldates)/12),rep(6,length(alldates)/12),rep(7,length(alldates)/12),rep(8,length(alldates)/12),
                 rep(9,length(alldates)/12),rep(10,length(alldates)/12),rep(11,length(alldates)/12),rep(12,length(alldates)/12))
  datedf <- data.frame(chamber_n=chamber_n,DateTime=alldates)                                                                                             
  
  #merge data in with dataframe of all days, sort it
  dat2 <- merge(dat,datedf,all=T,by=c("chamber_n","DateTime"))
  dat2 <- dat2[with(dat2,order(chamber_n,DateTime)),]
  
  #break across list, gapfill list
  dat3 <- zoo(dat2)
  dat3$Days_since_transplanting <- na.approx(dat3$Days_since_transplanting)
  dat3$vol <- na.approx(dat3$vol)
  dat3$vol_wood <- na.approx(dat3$vol_wood)
  dat3$vol_bark <- na.approx(dat3$vol_bark)
  dat3$mass_wood <- na.approx(dat3$mass_wood)
  dat3$mass_bark <- na.approx(dat3$mass_bark)
  
  
  
  # get it back to a normal dataframe
  dat4 <- numericdfr(fortify.zoo(dat3))
  dat4$DateTime <- as.Date(dat4$DateTime)
  dat4$Date <- dat4$DateTime
  
  # put some of the other bits back together
  dat4$chamber <- as.factor(paste0("C",sprintf("%02.0f",dat4$chamber_n)))
  dat4$T_treatment <- as.factor(ifelse(dat4$chamber_n %% 2 ==1,"ambient","elevated"))
  dat4$Index <- NULL
  dat4$Stem_number <- NULL
  
  
  #- establish a treatment key for the Water_treatment variable
  key <- data.frame(chamber=levels(dat4$chamber),
                    Water_treatment=c("drydown","control","drydown","drydown","control","drydown",
                                      "control","drydown","control","control","drydown","control"))
  
  size_before <- subset(dat4,Date<as.Date("2014-02-4"))
  size_before$Water_treatment <- "control"
  
  size_after <- subset(dat4,Date>=as.Date("2014-02-4"))
  size_after$Water_treatment <- NULL
  size_after2 <- merge(size_after,key,by="chamber")
  
  #- combined dataframes from before and after the drought began
  dat5 <- rbind(size_before,size_after2)
  
  #- clean up for exporting
  dat5_out <- dat5[,c("chamber","Date","T_treatment","Water_treatment","Measurement","Days_since_transplanting","vol","vol_wood","vol_bark","mass_wood","mass_bark")]
  
  return(dat5_out)
  
}

#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Estimate branch biomass through time in WTC3.
interpolateBranchDW <- function(vol.all=vol.all,plotson=0){
  
  bc <- read.csv("raw_data/WTC_TEMP_CM_BRANCHCENSUS_20130910-20140516_L0_v1.csv")
  bc$Date <- as.POSIXct(bc$Date,format="%d/%m/%Y",tz="GMT")
  bc$Date <- as.Date(bc$Date)
  bc$branchid <- as.factor(paste(bc$chamber,bc$branchnumber,sep="-"))
  
  #--------------------------------------------------------------------------------------------------
  # okay, so now can we estimate branch mass allometrically at the harvest? 
  branches <- read.csv("raw_data/WTC_TEMP_CM_GX-RBRANCH_20140513-20140522_L1_v1.csv")
  
  br.allom <- lm(log10(branch_dw)~log10(diam.m),data=branches) # fit log-log branch allometry
  if(plotson==1){
    # windows(12,12)
    x11(12,12)
    plotBy(branch_dw~diam.m|T_treatment,data=branches,log="xy",pch=15,col=c("black","red"),cex.lab=1.5,cex.axis=1.3,
           xlab= "Branch diameter (mm)",ylab="Branch mass (g)");abline(br.allom)
  }
  br.coef <- coef(br.allom) # get the coefficients
  
  #apply allometry to the branch census dataset
  bc$branch_dw <- 10^(br.coef[1]+br.coef[2]*log10(bc$branchdiameter))
  bc <- subset(bc,branch_dw>0)
  
  #get the layer heights
  layers <- read.csv("raw_data/WTC_TEMP_CM_LAYERHEIGHT_20140519_L1.csv")
  cuts <- c(0,mean(layers$toplayer1fromfloor),mean(layers$toplayer2fromfloor))
  layer.df <- data.frame(layer=c("low","mid","top"),layerlimit=cuts)
  
  bc$layer <- ifelse(bc$heightinsertion<=cuts[2],"low",ifelse(bc$heightinsertion<=cuts[3],"mid","top"))
  #sum up estimate of branch mass for each date for each chamber
  bc.sum <- summaryBy(branch_dw~chamber+Date+layer,data=bc,FUN=sum,keep.names=T)
  palette(c("black","blue","green","forestgreen","darkmagenta","chocolate1","cyan","darkgrey","darkgoldenrod1",
            "brown1","darkred","deeppink"))
  #compare branch allometry number to the harvest number
  harvest <- read.csv("raw_data/WTC_TEMP_CM_HARVEST-CANOPY_20140526-20140528_L1_v1.csv")
  harvest.sum <- summaryBy(BranchDM~chamber+layer,data=harvest,FUN=sum,keep.names=T)
  harvest.sum$branch_allom <- subset(bc.sum,Date==as.Date("2014-05-16"))$branch_dw
  
  if(plotson==1){
    X11(12,12)
    plotBy(branch_allom~BranchDM|layer,data=harvest.sum,xlab="Measured branch mass (g)",ylab="Allometric estimate (g)",cex.lab=1.3,xlim=c(0,1000),ylim=c(0,1000),pch=15);abline(0,1)
    textxy(X=harvest.sum$BranchDM,Y=harvest.sum$branch_allom,labs=harvest.sum$chamber,cex=0.9)
  }
  #--------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------
  
  
  bc.sum <- summaryBy(branch_dw~chamber+Date,data=bc,FUN=sum,keep.names=T)
  
  size2 <- merge(bc.sum,vol.all,by=c("Date","chamber"),all=T)
  
  #there is a nice relationship between tree (stem) volume and total branch mass
  if (plotson==1){
    X11(12,12)
    plotBy(branch_dw~vol|chamber,pch=15,type="b",cex=1.2,data=size2,log="xy",xlim=c(300,21000),
           xlab="Stem volume (cm3)",ylab="Branch mass (g)",cex.lab=1.3)
  }
  
  #- fit log-linear regression for each chamber, but all at once
  lm.all <- lm(log10(branch_dw)~log10(vol)*chamber,data=size2)
  
  # fit log- linear regression for each chamber
  dat.l <- split(size2,size2$chamber)
  lms1 <- list()
  for (i in 1:length(dat.l)){
    tofit <- dat.l[[i]]
    lms1[[i]] <- lm(log10(branch_dw)~log10(vol),data=subset(tofit,branch_dw>0))
    if (plotson==1) abline(lms1[[i]],col=i)
  }
  lms1 <- as.data.frame(do.call(rbind,lapply(lms1,coef)))
  lms1$chamber <- levels(size2$chamber)
  names(lms1) <- c("int","slope","chamber")
  
  
  # merge allometric estimates into the large size dataframe
  size3 <- merge(size2,lms1,by="chamber",all=T)
  size3$branch_dw_est <- with(size3,10^(int+slope*log10(vol)))
  #--------------------------------------------------------------------------------------------------
  
  if (plotson==1){
    X11(12,12)
    plotBy(branch_dw_est~Date|chamber,data=size3,ylab="Branch mass (g)",cex.lab=1.3)
  }
  
  
  out <- size3[,c("chamber","Date","branch_dw_est")]
  out <- out[with(out,order(chamber,Date)),]
  names(out)[3] <- "branchMass"
  return(out)
}
#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#- plot metrics of size over time.

plotsize <- function(output=T){
  
  
  
  #------------------------------------------------------------------------------------------------------------
  # #- get initial tree size 
  size <- read.csv("raw_data/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV")
  size$chamber_n <- as.numeric(substr(size$chamber,start=2,stop=3))
  size$DateTime <- as.Date(size$DateTime)
  
  summaryBy(Plant_height~T_treatment,data=subset(size,DateTime==as.Date("2012-12-12")),FUN=c(mean,standard.error))
  #------------------------------------------------------------------------------------------------------------
  
  
  
  
  #------------------------------------------------------------------------------------------------------------
  #- Get the three direct observations of leaf area. They happened on 9 Sept 2013, 10 Feb 2014, and
  #    during the harvest at ~25 May 2014.
  treeMass <- read.csv("raw_data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
  treeMass.sub <- subset(treeMass,as.Date(DateTime) %in% as.Date(c("2013-09-14","2014-02-10","2014-05-27")))
  treeMass.sub$Date <- as.Date(treeMass.sub$DateTime)
  leafArea <- summaryBy(leafArea~Date+chamber+T_treatment+Water_treatment,data=treeMass.sub,FUN=c(mean),keep.names=T)
  leafArea1 <- summaryBy(leafArea~Date+T_treatment+Water_treatment,data=leafArea,FUN=c(mean,standard.error))
  
  #-------------------
  #- get an estimate of volume and mass for wood and bark for each measurement day
  vol <- getvol()
  vol$diam <- vol$diam/10      # convert to cm
  vol$height <- vol$height/100 # convert to m
  vol$vol <- vol$vol/10000     # convert to m3
  
  size.m <- summaryBy(diam+height+vol~DateTime+T_treatment+Water_treatment,
                      data=subset(vol,Days_since_transplanting>0),FUN=c(mean,standard.error))
  
  
  
  # windows(40,40);
  par(mfrow=c(2,2),mar=c(0,5.5,0,3),oma=c(5,2,2,0),las=1,cex.axis=1.5)
  palette(c("#377EB8","#E41A1C"))
  
  #- plot diameter
  xlims <-as.Date(c("2013-3-1","2014-6-1"))
  yearlims <-as.Date(c("2013-1-1","2014-1-1"))
  plotBy(diam.mean~DateTime|T_treatment,data=subset(size.m,Water_treatment=="control"),pch=16,type="o",ylim=c(0,8),
         xlim=xlims,
         xaxt="n",yaxt="n",xlab="",ylab="",legend=F)
  adderrorbars(x=subset(size.m,Water_treatment=="control")$DateTime,
               y=subset(size.m,Water_treatment=="control")$diam.mean,barlen=0.0,
               SE=subset(size.m,Water_treatment=="control")$diam.se,direction="updown",
               col=c("#377EB8","#E41A1C"))       
  plotBy(diam.mean~DateTime|T_treatment,data=subset(size.m,Water_treatment=="drydown"),pch=1,type="o",add=T,legend=F)
  adderrorbars(x=subset(size.m,Water_treatment=="drydown")$DateTime,
               y=subset(size.m,Water_treatment=="drydown")$diam.mean,barlen=0.0,
               SE=subset(size.m,Water_treatment=="drydown")$diam.se,direction="updown",
               col=c("#377EB8","#E41A1C"))    
  magaxis(side=c(2,4),labels=c(1,1),frame.plot=T,majorn=3,las=1)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="month"),tcl=0.25,labels=F)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="quarter"),tcl=0.75,labels=F)
  title(ylab=expression(Diameter~(cm)),cex.lab=2)
  abline(v=as.Date("2013-9-13"),lty=2)
  #text(x=as.Date("2013-9-13"),y=8.8,labels="Flux measurements begin",xpd=NA,cex=1.3)
  legend("topleft",pch=c(16,16,1,1),lty=c(1),col=c("#377EB8","#E41A1C"),seg.len=1.5,
         legend=c("A-Wet","W-Wet","A-Dry","W-Dry"),bty="n",cex=1.2)
  legend("bottomright","a",bty="n",inset=0.002,cex=1.2)
  
  
  #-- plot leaf area over time
  plotBy(leafArea.mean~Date|T_treatment,data=subset(leafArea1,Water_treatment=="control"),pch=16,type="p",ylim=c(0,25),cex=1.5,
         xlim=xlims,
         xaxt="n",yaxt="n",xlab="",ylab="",legend=F)
  adderrorbars(x=subset(leafArea1,Water_treatment=="control")$Date,
               y=subset(leafArea1,Water_treatment=="control")$leafArea.mean,barlen=0.0,
               SE=subset(leafArea1,Water_treatment=="control")$leafArea.se,direction="updown",
               col=c("#377EB8","#E41A1C"))
  plotBy(leafArea.mean~Date|T_treatment,data=subset(leafArea1,Water_treatment=="drydown"),pch=1,type="p",add=T,legend=F,cex=1.5)
  adderrorbars(x=subset(leafArea1,Water_treatment=="drydown")$Date,
               y=subset(leafArea1,Water_treatment=="drydown")$leafArea.mean,barlen=0.0,
               SE=subset(leafArea1,Water_treatment=="drydown")$leafArea.se,direction="updown",
               col=c("#377EB8","#E41A1C"))
  magaxis(side=c(2,4),labels=c(1,1),frame.plot=T,majorn=3,las=1)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="month"),tcl=0.25,labels=F)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="quarter"),tcl=0.75,labels=F)
  title(ylab=expression(Total~leaf~area~(m^2)),cex.lab=2)
  abline(v=as.Date("2013-9-13"),lty=2)
  legend("bottomright","c",bty="n",inset=0.002,cex=1.2)
  
  
  #- plot height
  plotBy(height.mean~DateTime|T_treatment,data=subset(size.m,Water_treatment=="control"),pch=16,type="o",ylim=c(0,11),
         xlim=xlims,
         xaxt="n",yaxt="n",xlab="",ylab="",legend=F)
  adderrorbars(x=subset(size.m,Water_treatment=="control")$DateTime,
               y=subset(size.m,Water_treatment=="control")$height.mean,barlen=0.0,
               SE=subset(size.m,Water_treatment=="control")$height.se,direction="updown",
               col=c("#377EB8","#E41A1C"))
  plotBy(height.mean~DateTime|T_treatment,data=subset(size.m,Water_treatment=="drydown"),pch=1,type="o",add=T,legend=F)
  adderrorbars(x=subset(size.m,Water_treatment=="drydown")$DateTime,
               y=subset(size.m,Water_treatment=="drydown")$height.mean,barlen=0.0,
               SE=subset(size.m,Water_treatment=="drydown")$height.se,direction="updown",
               col=c("#377EB8","#E41A1C"))
  
  magaxis(side=c(2,4),labels=c(1,1),frame.plot=T,majorn=3,las=1)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="month"),tcl=0.25,labels=F)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="quarter"),tcl=0.75,labels=T,las=2,
            format="%b")
  title(ylab=expression(Height~(m)),cex.lab=2)
  abline(v=as.Date("2013-9-13"),lty=2)
  legend("bottomright","b",bty="n",inset=0.002,cex=1.2)
  
  
  #- plot stem volume
  plotBy(vol.mean~DateTime|T_treatment,data=subset(size.m,Water_treatment=="control"),pch=16,type="o",ylim=c(0,2),
         xlim=xlims,
         xaxt="n",yaxt="n",xlab="",ylab="",legend=F)
  adderrorbars(x=subset(size.m,Water_treatment=="control")$DateTime,
               y=subset(size.m,Water_treatment=="control")$vol.mean,barlen=0.0,
               SE=subset(size.m,Water_treatment=="control")$vol.se,direction="updown",
               col=c("#377EB8","#E41A1C"))
  plotBy(vol.mean~DateTime|T_treatment,data=subset(size.m,Water_treatment=="drydown"),pch=1,type="o",add=T,legend=F)
  adderrorbars(x=subset(size.m,Water_treatment=="drydown")$DateTime,
               y=subset(size.m,Water_treatment=="drydown")$vol.mean,barlen=0.0,
               SE=subset(size.m,Water_treatment=="drydown")$vol.se,direction="updown",
               col=c("#377EB8","#E41A1C"))
  
  magaxis(side=c(2,4),labels=c(1,1),frame.plot=T,majorn=3,las=1)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="month"),tcl=0.25,labels=F)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="quarter"),tcl=0.75,labels=T,las=2,
            format="%b")
  title(ylab=expression(Stem~volume~(m^3)),cex.lab=2)
  abline(v=as.Date("2013-9-13"),lty=2)
  legend("bottomright","d",bty="n",inset=0.002,cex=1.2)
  
  
  
  if(output==T) dev.copy2pdf(file="output/treeSize.pdf")
}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#- This is a collection of many functions that do the actual work of data analysis, run the model and figure plotting. 
# These functions are called by just a few lines of code in "CentralScript.R" to recreate the analyses and figures.

# Analysing Plant Storage (TNC) partitioning for Cstorage pool prediction
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
##### This function calculates TNC partitioning to tree organs (without considering organ biomass)
tnc.analysis <- function(carbohydrates,harvest) {
  # carbohydrates = read.csv("raw_data/Duan_carbohydrates.csv")
  carbohydrates = subset(carbohydrates,CO2 == 400 & Temp == "Amb" & Water == "Well watered")
  carbohydrates$tnc = carbohydrates$StarchW + carbohydrates$SolSugW # Unit = mg of tnc per g of dry weight biomass
  carbohydrates$tnc = carbohydrates$tnc / 10 # Unit = % of dry weight biomass
  
  # unit conversion from g of tnc per g of dry weight biomass to gC in tnc per gC in plant biomass
  # 1 g of tnc has 0.4 gC and 1 g of dry weight biomass has 0.48 gC
  carbohydrates$tnc = carbohydrates$tnc * (0.4/0.48) # Unit = gC in tnc / gC in plant biomass
  #-----------------------------------------------------------------------------------------
  ##### Total TNC calculation considering tree organ biomass partitioning
  # harvest = read.csv("raw_data/Duan_harvest.csv")
  harvest = subset(harvest,CO2 == 400 & Temp == "Amb" & Water == "Well watered")
  
  leaf.tnc = subset(carbohydrates,Organ == "Leaf") # Unit = % of dry weight leafmass
  stem.tnc = subset(carbohydrates,Organ == "Stem") # Unit = % of dry weight stemmass
  root.tnc = subset(carbohydrates,Organ == "Root") # Unit = % of dry weight rootmass
  
  tnc = data.frame(harvest$LeafDW,leaf.tnc$tnc,harvest$StemDW,stem.tnc$tnc,harvest$RootDW,root.tnc$tnc)
  names(tnc) <- c("leaf.C","leaf.tnc.C","stem.C","stem.tnc.C","root.C","root.tnc.C") 
  
  tnc$total.leaf.tnc.C = tnc$leaf.tnc.C * tnc$leaf.C / 100 # Unit = gC
  tnc$total.stem.tnc.C = tnc$stem.tnc.C * tnc$stem.C / 100 # Unit = gC
  tnc$total.root.tnc.C = tnc$root.tnc.C * tnc$root.C / 100 # Unit = gC
  
  tnc$leaf_to_all = tnc$total.leaf.tnc.C / (tnc$total.leaf.tnc.C + tnc$total.stem.tnc.C + tnc$total.root.tnc.C) * 100 # Unit = %
  tnc$stem_to_all = tnc$total.stem.tnc.C / (tnc$total.leaf.tnc.C + tnc$total.stem.tnc.C + tnc$total.root.tnc.C) * 100 # Unit = %
  tnc$root_to_all = tnc$total.root.tnc.C / (tnc$total.leaf.tnc.C + tnc$total.stem.tnc.C + tnc$total.root.tnc.C) * 100 # Unit = %
  
  tnc[nrow(tnc)+1, ] = colMeans(tnc, na.rm = TRUE) # R7 = Average of data
  tnc[nrow(tnc)+1, ] = apply(tnc, 2, sd) # R8 = Standard deviation of data
  dimnames(tnc)[[1]] <- c(1:6, "mean", "SD")
  
  return(tnc)
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# This script calcualtes LogLikelihood to find the most accurate model
logLikelihood <- function (data,output,model.comparison) {
  logLi <- matrix(0, nrow=nrow(data), ncol = 1) # Initialising the logLi
  for (i in 1:nrow(data)) {
    if (!is.na(data$Mleaf[i])) {
      logLi[i] = - 0.5*((output$Mleaf[i] - data$Mleaf[i])/data$Mleaf_SD[i])^2 - log(data$Mleaf_SD[i]) - log(2*pi)^0.5
    }
    if (!is.na(data$Mstem[i])) {
      logLi[i] = logLi[i] - 0.5*((output$Mstem[i] - data$Mstem[i])/data$Mstem_SD[i])^2 - log(data$Mstem_SD[i]) - log(2*pi)^0.5
    }
    if (!is.na(data$Mroot[i])) {
      logLi[i] = logLi[i] - 0.5*((output$Mroot[i] - data$Mroot[i])/data$Mroot_SD[i])^2 - log(data$Mroot_SD[i]) - log(2*pi)^0.5
    }
    if (model.comparison==F) {
      if (!is.null(data$Sleaf)) {
        if (!is.na(data$Sleaf[i])) {
          logLi[i] = logLi[i] - 0.5*((output$Sleaf[i] - data$Sleaf[i])/data$Sleaf_SD[i])^2 - log(data$Sleaf_SD[i]) - log(2*pi)^0.5
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
# various temporal scales (e.g. 1,2,...,121 days) to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot) and fluxes
#-------------------------------------------------------------------------------------
CBM.grouping <- function(chainLength, no.param.par.var, vol.group, with.storage, model.comparison, model.optimization) {
  
  # Assign inputs for MCMC
  bunr_in = chainLength * 0.1 # Discard the first 10% iterations for Burn-IN in MCMC (According to Oijen, 2008)
  if (with.storage==T) {
    no.var = 5 # variables to be modelled are: k,Y,af,as,sf
  } else {
    no.var = 4 # variables to be modelled are: k,Y,af,as
  }
  
  # Assign pot volumes and number of parameters per varible in temporal scale
  vol = unique(Cday.data.processed$volume)[order(unique(Cday.data.processed$volume))] # Assign all treatment pot volumes
  
  param.mean = data.frame(matrix(ncol = no.var+1, nrow = length(no.param.par.var)*length(vol.group)))
  if (with.storage==T) {
    names(param.mean) = c("k","Y","af","as","ar","sf")
  } else {
    names(param.mean) = c("Y","af","as","ar","sf")
  }
  aic.bic = data.frame(matrix(ncol = 7, nrow = length(no.param.par.var)*length(vol.group)))
  names(aic.bic) <- c("logLi","aic","bic","time","volume.group","no.param","volume")
  time = data.frame(no.param=rep(no.param.par.var,length(vol.group)),
                    start.time=numeric(length(no.param.par.var)*length(vol.group)),
                    end.time=numeric(length(no.param.par.var)*length(vol.group)),
                    time.taken=numeric(length(no.param.par.var)*length(vol.group)))
  
  q = 0 # Indicates the iteration number
  
  # Start the iteration for different treatment group and number of parameters
  for (v1 in 1:length(vol.group)) {
    v = unlist(vol.group[v1])
    # This script take the subset of processed data for particular treatment group
    source("R/data_processing_CBM.R", local=TRUE)
    
    for (z in 1:length(no.param.par.var)) {
      # Initialize few output data files
      q = q + 1
      time$start.time[q] <- Sys.time()
      param.vary = ceiling(nrow(data)/no.param.par.var[z]) # How many days the parameter set remain unchanged (weekly = 7; monthly = 30; just one parameter = nrow(data))
      no.param = ceiling(nrow(data)/param.vary) # number of parameter set for the whole duration of experiment (121 days)
      
      # This script initializes the parameter setting
      source("R/parameter_setting.R", local=TRUE) # initialize 'sf' prior differently for grouped treatments
      
      # Defining the variance-covariance matrix for proposal generation
      vcov = (0.005*(pMaxima-pMinima))^2
      vcovProposal =  vcov # The higher the coefficient, the higher the deviations in parameter time series
      
      
      # Find the Prior probability density
      prior.dist = vector("list", no.var)
      for (i in 1:no.var) {
        prior.dist[i] = list(log(dnorm(pValues[ , i], (pMinima[ , i] + pMaxima[ , i])/2, (pMaxima[ , i] - pMinima[ , i])/3))) # Prior normal gaussian distribution
      }
      logPrior0 <- sum(unlist(prior.dist))
      
      # Calculating model outputs for the starting point of the chain
      for (j in 1:length(v)) {
        data.set = subset(data,(volume %in% vol[v[j]]))
        Mleaf = Mstem = Mroot = LA = c()
        Mleaf[1] <- data.set$Mleaf[1]
        Mstem[1] <- data.set$Mstem[1]
        Mroot[1] <- data.set$Mroot[1]
        
        if (with.storage==T) {
          output.set = model(data.set$GPP,data.set$Rd,no.param,Mleaf,Mstem,Mroot,pValues$Y,pValues$k,pValues$af,pValues$as,pValues$sf)
        } else {
          output.set = model.without.storage(data.set$GPP,data.set$Rd,no.param,Mleaf,Mstem,Mroot,pValues$Y,pValues$af,pValues$as,pValues$sf)
        }
        
        output.set$volume = as.factor(vol[v[j]])
        if (j == 1) {
          output = output.set
        }
        if (j > 1) {
          output = rbind(output,output.set)
        }
      }
      
      data = data[order(data$volume),]
      logL0 <- logLikelihood(data,output,model.comparison) # Calculate log likelihood of starting point of the chain
      
      if (with.storage==T) {
        pChain[1,] <- c(pValues$k,pValues$Y,pValues$af,pValues$as,pValues$sf,logL0) # Assign the first parameter set with log likelihood
      } else {
        pChain[1,] <- c(pValues$Y,pValues$af,pValues$as,pValues$sf,logL0) # Assign the first parameter set with log likelihood
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
          names(candidatepValues) <- c("k","Y","af","as","sf")
        } else {
          names(candidatepValues) <- c("Y","af","as","sf")
        }
        
        # Reflected back to generate another candidate value
        reflectionFromMin = pmin( unlist(matrix(0,nrow=no.param,ncol=no.var)), unlist(candidatepValues-pMinima) )
        reflectionFromMax = pmax( unlist(list(rep(0, no.param))), unlist(candidatepValues-pMaxima) )
        candidatepValues = candidatepValues - 2 * reflectionFromMin - 2 * reflectionFromMax 
        
        
        # Calculating the prior probability density for the candidate parameter vector
        if (all(candidatepValues>pMinima) && all(candidatepValues<pMaxima)){
          uni.dist = vector("list", no.var)
          for (i in 1:no.var) {
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
            data.set = subset(data,(volume %in% vol[v[j]]))
            Mleaf = Mstem = Mroot = c()
            Mleaf[1] <- data.set$Mleaf[1]
            Mstem[1] <- data.set$Mstem[1]
            Mroot[1] <- data.set$Mroot[1]
            
            
            if (with.storage==T) {
              out.cand.set = model(data.set$GPP,data.set$Rd,no.param,Mleaf,Mstem,Mroot,candidatepValues$Y,
                                   candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf)
            } else {
              out.cand.set = model.without.storage(data.set$GPP,data.set$Rd,no.param,Mleaf,Mstem,Mroot,candidatepValues$Y,
                                                   candidatepValues$af,candidatepValues$as,candidatepValues$sf)
            }
            
            out.cand.set$volume = as.factor(vol[v[j]])
            if (j == 1) {
              out.cand = out.cand.set
            }
            if (j > 1) {
              out.cand = rbind(out.cand,out.cand.set)
            }
          }
          
          data = data[order(data$volume),]
          logL1 <- logLikelihood(data,out.cand,model.comparison) # Calculate log likelihood
          
          # Calculating the logarithm of the Metropolis ratio
          logalpha <- (logPrior1+logL1) - (logPrior0+logL0) 
          
          # Accepting or rejecting the candidate vector
          if ( log(runif(1, min = 0, max =1)) < logalpha && candidatepValues$af[1] + candidatepValues$as[1] <= 1
               && candidatepValues$as[1] >= 0 && candidatepValues$af[1] >= 0) {
            pValues <- candidatepValues
            logPrior0 <- logPrior1
            logL0 <- logL1
          }
        }
        if (with.storage==T) {
          pChain[c,] <- c(pValues$k,pValues$Y,pValues$af,pValues$as,pValues$sf,logL0)
        } else {
          pChain[c,] <- c(pValues$Y,pValues$af,pValues$as,pValues$sf,logL0)
        }
      }
      
      # Discard the first 500 iterations for Burn-IN in MCMC
      pChain <- pChain[(bunr_in+1):nrow(pChain),]
      pChain = as.data.frame(pChain)
      if (with.storage==T) {
        if (no.param.par.var[z]==1) {
          names(pChain) <- c("k1","Y1","af1","as1","sf1","logli")
        } else if (no.param.par.var[z]==2) {
          names(pChain) <- c("k1","k2","Y1","Y2","af1","af2","as1","as2","sf1","sf2","logli")
        } else if (no.param.par.var[z]==3) {
          names(pChain) <- c("k1","k2","k3","Y1","Y2","Y3","af1","af2","af3","as1","as2","as3","sf1","sf2","sf3","logli")
        }
      } else {
        if (no.param.par.var[z]==1) {
          names(pChain) <- c("Y1","af1","as1","sf1","logli")
        } else if (no.param.par.var[z]==2) {
          names(pChain) <- c("Y1","Y2","af1","af2","as1","as2","sf1","sf2","logli")
        } else if (no.param.par.var[z]==3) {
          names(pChain) <- c("Y1","Y2","Y3","af1","af2","af3","as1","as2","as3","sf1","sf2","sf3","logli")
        }
      }
      
      
      # Store the final parameter set values
      param.set = colMeans(pChain[ , 1:(no.param*no.var)])
      param.SD = apply(pChain[ , 1:(no.param*no.var)], 2, sd)
      param.final = data.frame(matrix(ncol = (no.var)*2, nrow = no.param))
      if (with.storage==T) {
        names(param.final) <- c("k","Y","af","as","sf","k_SD","Y_SD","af_SD","as_SD","sf_SD")
        param.final$k = param.set[1:no.param]
        param.final$Y = param.set[(1+no.param):(2*no.param)]
        param.final$af = param.set[(1+2*no.param):(3*no.param)]
        param.final$as = param.set[(1+3*no.param):(4*no.param)]
        param.final$sf = param.set[(1+4*no.param):(5*no.param)]
        
        param.final$k_SD = param.SD[1:no.param]
        param.final$Y_SD = param.SD[(1+no.param):(2*no.param)]
        param.final$af_SD = param.SD[(1+2*no.param):(3*no.param)]
        param.final$as_SD = param.SD[(1+3*no.param):(4*no.param)]
        param.final$sf_SD = param.SD[(1+4*no.param):(5*no.param)]
      } else {
        names(param.final) <- c("Y","af","as","sf","Y_SD","af_SD","as_SD","sf_SD")
        param.final$Y = param.set[1:no.param]
        param.final$af = param.set[(1+no.param):(2*no.param)]
        param.final$as = param.set[(1+2*no.param):(3*no.param)]
        param.final$sf = param.set[(1+3*no.param):(4*no.param)]
        
        param.final$Y_SD = param.SD[1:no.param]
        param.final$af_SD = param.SD[(1+no.param):(2*no.param)]
        param.final$as_SD = param.SD[(1+2*no.param):(3*no.param)]
        param.final$sf_SD = param.SD[(1+3*no.param):(4*no.param)]
      }
      
      # Calculate final output set from the predicted parameter set
      for (j in 1:length(v)) {
        data.set = subset(data,(volume %in% vol[v[j]]))
        Mleaf = Mstem = Mroot = c()
        Mleaf[1] <- data.set$Mleaf[1]
        Mstem[1] <- data.set$Mstem[1]
        Mroot[1] <- data.set$Mroot[1]
        
        if (with.storage==T) {
          output.final.set = model(data.set$GPP,data.set$Rd,no.param,Mleaf,Mstem,Mroot,param.final$Y,
                                   param.final$k,param.final$af,param.final$as,param.final$sf)
        } else {
          output.final.set = model.without.storage(data.set$GPP,data.set$Rd,no.param,Mleaf,Mstem,Mroot,param.final$Y,
                                                   param.final$af,param.final$as,param.final$sf)
        } 
        
        output.final.set$volume = as.factor(vol[v[j]])
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
      
      # Calculate daily parameter values with SD
      Days <- seq(1,nrow(data.set), length.out=nrow(data.set))
      param.daily = param.final[1,]
      
      if (no.param == 1) {
        for (i in 2:length(Days)) {
          param.daily[i,] = param.final[1,]
        }
      }
      if (no.param == 2) {
        for (i in 2:length(Days)) {
          param.daily[i,1:no.var] = param.final[1,1:no.var] + param.final[2,1:no.var] * i
        }
        for (i in (no.var+1):(2*no.var)) {
          param.daily[,i] = ((param.final[1,i]^2 + param.final[2,i]^2)/2)^0.5
        }
      }
      if (no.param == 3) {
        for (i in 2:length(Days)) {
          param.daily[i,1:no.var] = param.final[1,1:no.var] + param.final[2,1:no.var] * i + param.final[3,1:no.var] * i^2
        }
        for (i in (no.var+1):(2*no.var)) {
          param.daily[,i] = ((param.final[1,i]^2 + param.final[2,i]^2 + param.final[3,i]^2)/3)^0.5
        }
      }
      param.daily$ar = 1 - param.daily$af - param.daily$as
      param.daily$ar_SD = with(param.daily, ((af_SD*af_SD + as_SD*as_SD)/2)^0.5)
      param.daily$Date = as.Date(data.set$Date)
      
      
      # Plotting the parameter sets over time
      if (with.storage==T) { 
        melted.param1 = melt(param.daily[,c("k","Y","af","as","ar","sf","Date")], id.vars="Date")
        melted.param2 = melt(param.daily[,c("k_SD","Y_SD","af_SD","as_SD","ar_SD","sf_SD","Date")], id.vars="Date")
      } else {
        melted.param1 = melt(param.daily[,c("Y","af","as","ar","sf","Date")], id.vars="Date")
        melted.param2 = melt(param.daily[,c("Y_SD","af_SD","as_SD","ar_SD","sf_SD","Date")], id.vars="Date")
      }
      melted.param = data.frame(melted.param1$Date, melted.param1$variable, melted.param1$value, melted.param2$value)
      names(melted.param) = c("Date","variable","Parameter","Parameter_SD")
      melted.param$Date = as.Date(melted.param$Date)
      melted.param$volume = list(vol[unlist(vol.group[v1])])
      melted.param$volume.group = as.factor(v1)
      melted.param$no.param = as.factor(no.param.par.var[z])
      
      
      # Plotting the Measured (data) vs Modelled Plant Carbon pools for plotting and comparison
      #----------------------------------------------------------------------------------------------------------------
      # lm.daily = read.csv("processed_data/Cleaf_daily_data.csv") # Unit gC per gC plant
      #----------------------------------------------------------------------------------------------------------------
      
      for (j in 1:length(v)) {
        data.set = subset(data,(volume %in% vol[v[j]]))
        
        #----------------------------------------------------------------------------------------------------------------
        # leafmass.daily = subset(lm.daily,(volume %in% vol[v[j]]))
        # leafmass.daily$Date = as.Date(leafmass.daily$Date)
        # # leafmass.daily.gas = leafmass.daily[leafmass.daily$Date %in% c(unique(as.Date(tnc.data.processed$Date))), ]
        # # browser()
        # data.set$Sleaf = data.set$Sleaf / leafmass.daily$leafmass * 100
        # data.set$Sleaf_SD = ((data.set$Sleaf_SD*data.set$Sleaf_SD + leafmass.daily$leafmass_SE*leafmass.daily$leafmass_SE)/2)^0.5 / leafmass.daily$leafmass * 100
        #----------------------------------------------------------------------------------------------------------------
        
        output.final.set = subset(output.final,(volume %in% vol[v[j]]))
        output.final.set$Date = data.set$Date
        
        if (with.storage==T) { 
          names(output.final.set) = c("Cstorage.modelled","Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled","volume","Date")
          melted.output = melt(output.final.set[,c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled","Date")], id.vars="Date")
          melted.data = melt(data.set[ , c("Mleaf","Mstem","Mroot","Sleaf","Date")], id.vars="Date")
          melted.error = melt(data.set[ , c("Mleaf_SD","Mstem_SD","Mroot_SD","Sleaf_SD","Date")], id.vars="Date")
        } else {
          names(output.final.set) = c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","volume","Date")
          melted.output = melt(output.final.set[,c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Date")], id.vars="Date")
          melted.data = melt(data.set[ , c("Mleaf","Mstem","Mroot","Date")], id.vars="Date")
          melted.error = melt(data.set[ , c("Mleaf_SD","Mstem_SD","Mroot_SD","Date")], id.vars="Date")
        }
        melted.output$Date = as.Date(melted.output$Date)
        melted.output$volume = as.factor(vol[v[j]])
        melted.output$no.param = as.factor(no.param.par.var[z])
        
        if (with.storage==T) { 
          melted.Cstorage = output.final.set[,c("Cstorage.modelled","Date")]
          melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
          melted.Cstorage$volume = as.factor(vol[v[j]])
          melted.Cstorage$no.param = as.factor(no.param.par.var[z])
        }
        
        melted.data$Date = as.Date(melted.data$Date)
        melted.data$volume = as.factor(vol[v[j]])
        
        melted.error$Date = as.Date(melted.error$Date)
        melted.error$volume = as.factor(vol[v[j]])
        melted.error$parameter = melted.data$value
        melted.error$no.param = as.factor(no.param.par.var[z])
        
        if (v1 < 8){
          melted.output$volume.group = as.factor(1)
          if (with.storage==T) { 
            melted.Cstorage$volume.group = as.factor(1)
          }
          melted.error$volume.group = as.factor(1)
        }
        if (v1 == 8){
          melted.output$volume.group = as.factor(2)
          if (with.storage==T) { 
            melted.Cstorage$volume.group = as.factor(2)
          }
          melted.error$volume.group = as.factor(2)
        }
        
        
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
        summary.error = rbind(summary.error,summary.error.set)
        if (z == 1) {
          summary.data = rbind(summary.data,summary.data.set)
        }
      }
      
      # Display the Acceptance rate of the chain
      nAccepted = length(unique(pChain[,1]))
      acceptance = (paste("Volume =",vol[v],", Total Parameter number =",no.param.par.var[z],": ", nAccepted, "out of ", chainLength-bunr_in, "candidates accepted ( = ",
                          round(100*nAccepted/chainLength), "%)"))
      print(acceptance)
      
      
      # # Plotting all parameter whole iterations for Day 1 only to check the convergance
      # png(file = paste("output/Parameter_iterations_day1_",v1,"_vol_",vol[v],"_par_",no.param.par.var[z], ".png", sep = ""))
      # par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
      # plot(pChain[,1],col="red",main="Utilization coefficient at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="k",ylim=c(param.k[1,1],param.k[1,3]))
      # plot(pChain[,1+no.param],col="green",main="Alloc frac to Biomass at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="Y",ylim=c(param.Y[1,1],param.Y[1,3]))
      # plot(pChain[,1+2*no.param],col="magenta",main="Alloc frac to foliage at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="af",ylim=c(param.af[1,1],param.af[1,3]))
      # plot(pChain[,1+3*no.param],col="blue",main="Alloc frac to stem at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="as",ylim=c(param.as[1,1],param.as[1,3]))
      # plot(pChain[,1+4*no.param],col="green",main="Foliage turnover at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="sf",ylim=c(param.sf[1,1],param.sf[1,3]))
      # plot(pChain[,1+5*no.param],col="magenta",main="Log-likelihood",cex.lab = 1.5,xlab="Iterations",ylab="Log-likelihood")
      # title(main = paste("First day Parameter iterations for volume group",v1,"with par",no.param.par.var[z]), outer=TRUE, cex = 1.5)
      # dev.off()
      
      # Calcualte LogLi, AIC, BIC, Time to find the most accurate model for best balance between model fit and complexity
      output.final1 = output.final
      if (with.storage==T) { 
        names(output.final1) = c("Cstorage","Mleaf","Mstem","Mroot","Sleaf","volume") # Rename for the logLikelihood function
      } else {
        names(output.final1) = c("Mleaf","Mstem","Mroot","volume")
      }
      
      data = data[with(data, order(volume)), ]
      row.names(data) = c(1:nrow(data))
      aic.bic[q,1] <- logLikelihood(data,output.final1,model.comparison) # Calculate logLikelihood
      
      k1 = 2 # k = 2 for the usual AIC
      npar = no.param*no.var # npar = total number of parameters in the fitted model
      aic.bic[q,2] = -2*aic.bic[q,1] + k1*npar
      
      if (model.comparison==F) {
        n = sum(!is.na(data$Sleaf)) + sum(!is.na(data$Mleaf)) + sum(!is.na(data$Mstem)) + sum(!is.na(data$Mroot))
      } else {
        n = sum(!is.na(data$Mleaf)) + sum(!is.na(data$Mstem)) + sum(!is.na(data$Mroot))
      }
      k2 = log(n) # n being the number of observations for the so-called BIC
      aic.bic[q,3] = -2*aic.bic[q,1] + k2*npar
      
      time$end.time[q] <- Sys.time()
      time$time.taken[q] <- time$end.time[q] - time$start.time[q]
      aic.bic[q,4] = time$time.taken[q]
      aic.bic$volume.group[q] = v1
      aic.bic[q,6] = no.param.par.var[z]
      aic.bic$volume[q] = list(vol[unlist(vol.group[v1])])
    }
  }
  bic = data.frame(aic.bic[,c("bic","volume.group","no.param","volume")])
  melted.aic.bic = melt(aic.bic[,c(1:6)], id.vars=c("no.param","volume.group"))
  
  # if (model.comparison==T | model.optimization==T) {
  if (model.optimization==T) {
    return(bic)
  } else if (model.optimization==F & with.storage==T) {
    result = list(no.param.par.var,summary.param,summary.data,summary.output,summary.error,bic,summary.Cstorage)
    return(result)
  } else {
    result = list(no.param.par.var,summary.param,summary.data,summary.output,summary.error,bic)
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
# various temporal scales to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot)
#-----------------------------------------------------------------------------------------
# Defining the model to iteratively calculate Cstorage, Cleaf, Cstem, Croot, Sleaf, Sstem, Sroot
model <- function (GPP,Rd,no.param,Mleaf,Mstem,Mroot,Y,k,af,as,sf) {
  Cstorage = Sleaf = Sstem = Sroot = c()
  
  # From Duan's experiment for TNC partitioning to tree organs
  # Leaf TNC C / Leaf C =  0.1167851; Stem TNC C / Stem C =  0.03782242; Root TNC C / Root C =  0.01795031
  Sleaf[1] = Mleaf[1] * tnc$leaf.tnc.C[7]/100
  Sstem[1] = Mstem[1] * tnc$stem.tnc.C[7]/100
  Sroot[1] = Mroot[1] * tnc$root.tnc.C[7]/100
  Cstorage[1] <- Sleaf[1] + Sstem[1] + Sroot[1] 
  
  Cleaf <- Croot <- Cstem <- c()
  Cleaf[1] <- Mleaf[1] - Sleaf[1]
  Cstem[1] <- Mstem[1] - Sstem[1]
  Croot[1] <- Mroot[1] - Sroot[1]
  
  if (no.param == 1) {
    k.i = k[1]; Y.i = Y[1]; af.i = af[1]; as.i = as[1]; sf.i = sf[1]
  }
  for (i in 2:length(GPP)) {
    if (no.param == 2) {
      k.i = k[1] + k[2]*i; Y.i = Y[1]+ Y[2]*i; af.i = af[1]+ af[2]*i; as.i = as[1]+ as[2]*i; sf.i = sf[1]+ sf[2]*i
    }
    if (no.param == 3) {
      k.i = k[1] + k[2]*i + k[3]*i*i; Y.i = Y[1]+ Y[2]*i + Y[3]*i*i; af.i = af[1]+ af[2]*i + af[3]*i*i; 
      as.i = as[1]+ as[2]*i + as[3]*i*i; 
      sf.i = sf[1]+ sf[2]*i + sf[3]*i*i
    }
    Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1]) - k.i*Cstorage[i-1]
    Sleaf[i] <- Cstorage[i] * tnc$leaf_to_all[7]/100 # 75% of storage goes to leaf (Duan's experiment)
    Sstem[i] <- Cstorage[i] * tnc$stem_to_all[7]/100 # 16% of storage goes to stem (Duan's experiment)
    Sroot[i] <- Cstorage[i] * tnc$root_to_all[7]/100 # 9% of storage goes to root (Duan's experiment)
    
    Cleaf[i] <- Cleaf[i-1] + k.i*Cstorage[i-1]*af.i*(1-Y.i) - sf.i*Mleaf[i-1]
    Cstem[i] <- Cstem[i-1] + k.i*Cstorage[i-1]*as.i*(1-Y.i)
    Croot[i] <- Croot[i-1] + k.i*Cstorage[i-1]*(1-af.i-as.i)*(1-Y.i)
    
    Mleaf[i] <- Cleaf[i] + Sleaf[i]
    Mstem[i] <- Cstem[i] + Sstem[i]
    Mroot[i] <- Croot[i] + Sroot[i]
  }
  output = data.frame(Cstorage,Mleaf,Mstem,Mroot,Sleaf)
  
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
# various temporal scales to estimate Carbon pools (Cleaf,Cstem,Croot)
#-----------------------------------------------------------------------------------------
# This version does not consider storage pool

# Defining the model to iteratively calculate Cleaf, Cstem, Croot
model.without.storage <- function (GPP,Rd,no.param,Mleaf,Mstem,Mroot,Y,af,as,sf) {
  
  if (no.param == 1) {
    Y.i = Y[1]; af.i = af[1]; as.i = as[1]; sf.i = sf[1]
  }
  
  for (i in 2:length(GPP)) {
    if (no.param == 2) {
      Y.i = Y[1]+ Y[2]*i; af.i = af[1]+ af[2]*i; as.i = as[1]+ as[2]*i; sf.i = sf[1]+ sf[2]*i
    }
    if (no.param == 3) {
      Y.i = Y[1]+ Y[2]*i + Y[3]*i*i; af.i = af[1]+ af[2]*i + af[3]*i*i; 
      as.i = as[1]+ as[2]*i + as[3]*i*i; 
      sf.i = sf[1]+ sf[2]*i + sf[3]*i*i
    }
    
    Mleaf[i] <- Mleaf[i-1] + (GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1])) * af.i*(1-Y.i) - sf.i*Mleaf[i-1]
    Mstem[i] <- Mstem[i-1] + (GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1])) * as.i*(1-Y.i)
    Mroot[i] <- Mroot[i-1] + (GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1])) * (1-af.i-as.i)*(1-Y.i)
  }
  output = data.frame(Mleaf,Mstem,Mroot)
  return(output)
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# ################ Figure 1 #####################
# # Script to read and gerenate the model diagram
# plot.model <- function() { 
#   img <- readPNG("raw_data/Figure_1.png")
#   
#   #get size
#   h<-dim(img)[1]
#   w<-dim(img)[2]
#   
#   #open new file for saving the image in "output" folder
#   png("output/Figure_1_CBM.png", width=w, height=h)
#   par(mar=c(0,0,0,0), xpd=NA, mgp=c(0,0,0), oma=c(0,0,0,0), ann=F)
#   plot.new()
#   plot.window(0:1, 0:1)
#   
#   #fill plot with image
#   usr<-par("usr")    
#   rasterImage(img, usr[1], usr[3], usr[2], usr[4])
#   
#   #close image
#   dev.off()
# }


#----------------------------------------------------------------------------------------------------------------
################ Figure 1 #####################
# Script to read and gerenate the model diagram
plot.model.wtc3 <- function() { 
  img <- readPNG("raw_data/Figure_1.png")
  
  #get size
  h<-dim(img)[1]
  w<-dim(img)[2]
  
  #open new file for saving the image in "output" folder
  png("output/Figure_1_CBM_wtc3.png", width=w, height=h)
  par(mar=c(0,0,0,0), xpd=NA, mgp=c(0,0,0), oma=c(0,0,0,0), ann=F)
  plot.new()
  plot.window(0:1, 0:1)
  
  #fill plot with image
  usr<-par("usr")    
  rasterImage(img, usr[1], usr[3], usr[2], usr[4])
  #close image
  dev.off()
  grid.raster(img)
}
#----------------------------------------------------------------------------------------------------------------
plot.q10 <- function() { 
  img <- readPNG("output/5.Rdark_vs_T.png")
  
  # #get size
  # h<-dim(img)[1]
  # w<-dim(img)[2]
  # 
  # # #open new file for saving the image in "output" folder
  # # png("output/Figure_1_CBM_wtc3.png", width=w, height=h)
  # par(mar=c(0,0,0,0), xpd=NA, mgp=c(0,0,0), oma=c(0,0,0,0), ann=F)
  # plot.new()
  # plot.window(0:1, 0:1)
  # 
  # #fill plot with image
  # usr<-par("usr")    
  # rasterImage(img, usr[1], usr[3], usr[2], usr[4])
  # #close image
  # dev.off()
  grid.raster(img,  width = 4)
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

################ Figure 2 #####################
# Plot Model Measures ("bic") against "models with and without storage pool" 
# - 3 Grouped treatments with various number of parameters (Constant, Linear and Quadratic)
#-------------------------------------------------------------------------------------
plot.with.without.storage <- function(bic.with.storage, bic.without.storage) { 
  keeps = c("bic", "volume.group", "no.param")
  bic.with.storage = bic.with.storage[ , keeps, drop = FALSE]
  bic.without.storage = bic.without.storage[ , keeps, drop = FALSE]
  bic.with.storage.sub = subset(bic.with.storage, no.param %in% 3)
  bic.without.storage.sub = subset(bic.without.storage, no.param %in% 3)
  bic = data.frame(bic.without.storage.sub$bic, bic.with.storage.sub$bic, bic.with.storage.sub$volume.group)
  names(bic)[1:3] <- c("without storage", "with storage", "Treatment")
  bic$Group = ifelse(bic$Treatment==1,"Small",ifelse(bic$Treatment==2,"Large","Free"))
  bic$Group = factor(bic$Group, levels = bic$Group[c(1,2,3)])
  
  keeps = c("without storage", "with storage", "Group")
  bic = bic[ , keeps, drop = FALSE]
  bic.melt <- melt(bic, id.vars = "Group")
  names(bic.melt)[3] <- "BIC"
  
  p1 = ggplot(data = bic.melt, aes(x = Group, y = BIC, group = variable, fill = variable)) +
    # ggplot(data = bic.melt, aes(x=factor(bic.melt$Treatment, levels=unique(as.character(bic.melt$Treatment))), y=BIC, fill = Model_setting)) +
    geom_bar(stat="identity", width = 0.5, position = "dodge") +
    scale_fill_brewer(palette = "Set1", name = "") +
    ylab("BIC") +
    theme_bw() +
    theme(legend.title = element_text(colour="black", size=10)) +
    theme(legend.text = element_text(colour="black", size=10)) +
    theme(legend.position = c(0.75,0.82)) +
    theme(legend.key.height=unit(1,"line")) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=12)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = 12, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  png("output/Figure_2_bic.with.without.storage.png", units="px", width=800, height=600, res=220)
  print (p1)
  dev.off()
  
  # bic = data.frame(bic.without.storage$bic, bic.with.storage$bic, bic.with.storage$volume.group, bic.with.storage$no.param)
  # names(bic)[1:4] <- c("without storage", "with storage", "Treatment", "no.param")
  # bic$Group = ifelse(bic$Treatment==1,"Group 1",ifelse(bic$Treatment==2,"Group 2","Group 3"))
  # bic$Group = as.factor(bic$Group)
  # bic$Parameter.setting = ifelse(bic$no.param==1,"Constant",ifelse(bic$no.param==2,"Linear","Quadratic"))
  # bic$Parameter.setting = as.factor(bic$Parameter.setting)
  
  # keeps = c("without storage", "with storage", "Group", "Parameter.setting")
  # bic = bic[ , keeps, drop = FALSE]
  # bic.melt <- melt(bic, id.vars = c("Group", "Parameter.setting"))
  # names(bic.melt)[3:4] <- c("Model_setting", "BIC")
  #-------------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------------
  # p1 = ggplot(data = bic.melt, aes(x = Parameter.setting, y = BIC, group = Model_setting, fill = Model_setting)) +
  #   # ggplot(data = bic.melt, aes(x=factor(bic.melt$Treatment, levels=unique(as.character(bic.melt$Treatment))), y=BIC, fill = Model_setting)) +
  #   geom_bar(stat="identity", width = 0.5, position = "dodge") +
  #   facet_grid(. ~ Group) +
  #   scale_fill_brewer(palette = "Set1", name = "Model setting") +
  #   xlab("Parameter setting") +
  #   ylab("BIC") +
  #   # ggtitle("BIC for various model settings") +
  #   theme_bw() +
  #   # theme(plot.title = element_text(size = 12)) +
  #   theme(legend.title = element_text(colour="black", size=10)) +
  #   theme(legend.text = element_text(colour="black", size=10)) +
  #   theme(legend.position = c(0.85,0.8)) +
  #   theme(legend.key.height=unit(1,"line")) +
  #   theme(legend.key = element_blank()) +
  #   theme(text = element_text(size=12)) +
  #   theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  #   theme(axis.title.y = element_text(size = 12, vjust=0.3)) +
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # 
  # png("output/Figure_2_bic.with.without.storage.png", units="px", width=2000, height=1600, res=220)
  # print (p1)
  # dev.off()
}
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

################ Figure 3 #####################
# Plot Model Measures ("bic") against "Treatment groupings" and "number of parameters"
#-------------------------------------------------------------------------------------
plot.parameter.settings <- function(bic.group1, bic.group2, bic.group3, bic.group4) { 
  bic.group = bic.group1[rep(rownames(bic.group1), length(bic.group1$volume[[1]])), ]
  bic.group = bic.group[with(bic.group, order(no.param)), ]
  bic.group$volume = unlist(bic.group1$volume)
  
  bic.group2.order = bic.group2[rep(rownames(bic.group2[1:3,]), length(bic.group2$volume[[1]])), ]
  bic.group2.order = bic.group2.order[with(bic.group2.order, order(no.param)), ]
  bic.group2.order = rbind(bic.group2.order, bic.group2[bic.group2$volume.group==2, ])
  bic.group2.order$volume = unlist(bic.group2$volume)
  bic.group2.order$volume.group = 2
  bic.group = rbind(bic.group, bic.group2.order)
  
  bic.group3.order = bic.group3[rep(rownames(bic.group3[1:6,]), length(bic.group3$volume[[1]])), ]
  bic.group3.order = bic.group3.order[with(bic.group3.order, order(volume.group,no.param)), ]
  bic.group3.order = rbind(bic.group3.order, bic.group3[bic.group3$volume.group==3, ])
  bic.group3.order$volume = unlist(bic.group3$volume)
  bic.group3.order$volume.group = 3
  bic.group = rbind(bic.group, bic.group3.order)
  
  bic.group4$volume.group = 4
  bic.group = rbind(bic.group, bic.group4)
  bic.group$volume = unlist(bic.group$volume)
  bic.group$volume = as.factor(bic.group$volume)
  #-------------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------------
  pd <- position_dodge(0.3)
  p2 = ggplot(data = bic.group, aes(x = factor(bic.group$volume, levels=unique(as.character(bic.group$volume)) ), y = bic, group = interaction(volume.group,no.param), colour=factor(volume.group), shape=factor(no.param))) +
    geom_line(position=pd, data = bic.group, aes(x = factor(bic.group$volume, levels=unique(as.character(bic.group$volume)) ), y = bic, group = interaction(volume.group,no.param), colour=factor(volume.group), linetype=factor(no.param))) +
    geom_point(position=pd, size=2) +
    xlab("Treatments") +
    labs(colour="Grouping options", shape="Parameters", linetype="Parameters") +
    scale_y_continuous(name="BIC",limits = c(0, round_any(max(bic.group$bic), 50, f = ceiling)), breaks=seq(0,round_any(max(bic.group$bic), 50, f = ceiling),250)) +
    scale_x_discrete(name="Treatment", breaks=c("5", "10", "15", "20", "25", "35", "1000"),
                     labels=c("5L", "10L", "15L", "20L", "25L", "35L", "1000L")) +
    scale_color_manual(breaks = c("1", "2", "3", "4"), values=c("red", "green", "blue", "black")) +
    theme_bw() +
    theme(legend.title = element_text(colour="black", size=10)) +
    theme(legend.text = element_text(colour="black", size=10)) +
    theme(legend.key.height=unit(1,"line")) +
    theme(legend.key.width=unit(2,"line")) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=12)) +
    theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
    theme(axis.title.y = element_text(size = 12, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  png("output/Figure_3_bic_treat_group.png", units="px", width=1500, height=1200, res=220)
  print (p2)
  dev.off()
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

################ Figure 4 #####################
# Plot modelled parameters with 3 Grouped treatments and quadratic parameter setting
#-------------------------------------------------------------------------------------
plot.Modelled.parameters <- function(result,with.storage) {
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
  no.param.par.var = result[[1]]
  summary.param = result[[2]]
  summary.data = result[[3]]
  summary.output = result[[4]]
  summary.error = result[[5]]
  
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
plot.Modelled.biomass <- function(result,with.storage) { 
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


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
################ Figure 6A #####################
# Function to gerenate plot for Daily net C assimilation per unit leaf area (Cday) 
plot.Cday <- function(Cday.set, iteration) { 
  ggplot(data = Cday.set, aes(x = Date, y = carbon_day,  group = volume, colour=factor(volume))) +
    geom_point(size=0.01) +
    geom_line(data = Cday.set, aes(x = Date, y = carbon_day,  group = volume, colour=factor(volume))) +
    ylab(expression(C[day]~"(g C "*d^"-1"*" "*leafarea^"-1"*")")) +
    scale_colour_manual(name="", breaks=c("5", "1000"), labels=c("5L", "FS"), values=cbPalette[1:2]) +
    annotate("text", x = min(Cday.set$Date), y = max(Cday.set$carbon_day)*0.98, size = font.size-7, label = paste(title[1])) +
    theme_bw() +
    theme(legend.position = c(0.85,0.85)) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0, 0.45), units="line")) +
    theme(text = element_text(size=font.size)) +
    theme(legend.key.height=unit(0.65,"line")) +
    # theme(legend.key.width=unit(2,"line")) +
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(vjust=-2))
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
################ Figure 6B #####################
# Function to gerenate plot for Day respiration rate (Rd)
plot.Rd <- function(Rd.data, iteration) { 
  plot.shift[[2]] = ggplot(data = Rd.data, aes(x = Date, y = Rd_daily,  group = volume, colour=factor(volume))) +
    geom_point(size=0.01) +
    geom_line(data = Rd.data, aes(x = Date, y = Rd_daily,  group = volume, colour=factor(volume))) +
    # xlab("Month") +
    ylab(expression(R[d]~"(g C "*g^"-1"*" plant "*d^"-1"*")")) + 
    # ggtitle("B - Case 2: Rd (5L pot -> Free)") +
    scale_colour_manual(name="", breaks=c("5", "1000"), labels=c("5L", "FS"), values=cbPalette[2:3]) +
    # scale_y_continuous(limits = c(min(summary.param.set$Parameter), max(summary.param.set$Parameter)), breaks=c(.005,.01,.015),labels=c(.005,.01,.015)) +
    annotate("text", x = min(Rd.data$Date), y = max(Rd.data$Rd_daily)*0.98, size = font.size-7, label = paste(title[2])) +
    theme_bw() +
    theme(legend.position = c(0.85,0.85)) +
    # theme(plot.title = element_text(size = font.size)) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0, 0), units="line")) +
    theme(text = element_text(size=font.size)) +
    theme(legend.key.height=unit(0.65,"line")) +
    # theme(legend.key.width=unit(2,"line")) +
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
################ Figure 6C #####################
# Function to gerenate plot biomass allocation fractions (af, as, ar)
plot.allocation.fractions <- function(summary.param.set, iteration) { 
  ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = interaction(volume.group,variable), colour=factor(volume.group))) +
    # p0 = ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    # geom_point(size=0.01) +
    # geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = interaction(volume.group,variable), colour=factor(volume.group), linetype=factor(variable))) +
    # xlab("Month") +
    # ggtitle(paste(title[p],"- Case",iteration,":",as.character(var[p]),"(5L pot -> Free)")) +
    scale_colour_manual(breaks=c("1", "3"), labels=c("5L", "FS"), values=cbPalette[iteration:(iteration+1)]) +
    # scale_linetype_manual(values=c("solid","dashes", "dotted")) +
    scale_linetype_manual(breaks=c("af","as","ar"), labels=c(expression(a[f]),expression(a[w]),expression(a[r])),values=c(19,17,15)) +
    scale_y_continuous(name=expression(Allocations~"(g C "*g^"-1"*" C "*d^"-1"*")"),limits = c(0,0.9), breaks=seq(0,1,0.2)) +
    annotate("text", x = min(summary.param.set$Date), y = 0.86, size = font.size-7, label = paste(title[iteration])) +
    theme_bw() +
    theme(legend.position = c(0.47,0.83),legend.direction = "vertical",legend.box = "horizontal") +
    # theme(plot.title = element_text(size = font.size)) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0, 0.85), units="line")) +
    theme(text = element_text(size=font.size)) +
    theme(legend.key.height=unit(0.75,"line")) +
    theme(legend.key.width=unit(2,"line")) +
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
################ Figure 6D #####################
# Function to gerenate plot for growth respiration rate (Y) 
plot.Y <- function(summary.param.set, iteration) { 
  ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    geom_point(size=0.01) +
    geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    # geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume, colour=factor(volume))) +
    # xlab("Month") +
    # ggtitle(paste(title[p],"- Case",iteration,":",as.character(var[p]),"(5L pot -> Free)")) +
    scale_colour_manual(breaks=c("1", "3"), labels=c("5L", "FS"), values=cbPalette[iteration:(iteration+1)]) +
    scale_y_continuous(limits = c(min(summary.param.set$Parameter)*0.95, max(summary.param.set$Parameter)*1.05)) +
    annotate("text", x = min(summary.param.set$Date), y = max(summary.param.set$Parameter)*1.04, size = font.size-7, label = paste(title[iteration])) +
    ylab(expression(Y~"(g C "*g^"-1"*" C "*d^"-1"*")")) +
    theme_bw() +
    theme(legend.position = c(0.85,0.55)) +
    # theme(plot.title = element_text(size = font.size)) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0, 0.5), units="line")) +
    theme(text = element_text(size=font.size)) +
    theme(legend.key.height=unit(0.75,"line")) +
    # theme(legend.key.width=unit(3,"line")) +
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
################ Figure 6E #####################
# Function to gerenate plot for leaf turnover rate (sf)
plot.sf <- function(summary.param.set, iteration) { 
  ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    geom_point(size=0.01) +
    geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    # geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume, colour=factor(volume))) +
    # xlab("Month") +
    # ggtitle(paste(title[p],"- Case",iteration,":",as.character(var[p]),"(5L pot -> Free)")) +
    scale_colour_manual(breaks=c("1", "3"), labels=c("5L", "FS"), values=cbPalette[iteration:(iteration+1)]) +
    # scale_y_continuous(limits = c(min(summary.param.set$Parameter)*0.95, max(summary.param.set$Parameter)*1.05), breaks=c(.005,.007,.009),labels=c(.005,.007,.009)) +
    scale_y_continuous(limits = c(min(summary.param.set$Parameter)*0.95, max(summary.param.set$Parameter)*1.05)) +
    annotate("text", x = min(summary.param.set$Date), y = max(summary.param.set$Parameter)*1.04, size = font.size-7, label = paste(title[iteration])) +
    ylab(expression(s[f]~"(g C "*g^"-1"*" C "*d^"-1"*")")) +
    theme_bw() +
    theme(legend.position = c(0.85,0.3)) +
    # theme(plot.title = element_text(size = font.size)) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0, 0.15), units="line")) +
    theme(text = element_text(size=font.size)) +
    theme(legend.key.height=unit(0.75,"line")) +
    # theme(legend.key.width=unit(3,"line")) +
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
################ Figure 6F #####################
# Function to gerenate plot for utilization coefficient (K)
plot.k <- function(summary.param.set, iteration) { 
  ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    geom_point(size=0.01) +
    geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    # geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume, colour=factor(volume))) +
    # xlab("Month") +
    # ggtitle(paste(title[p],"- Case",iteration,":",as.character(var[p]),"(5L pot -> Free)")) +
    scale_colour_manual(breaks=c("1", "3"), labels=c("5L", "FS"), values=cbPalette[iteration:(iteration+1)]) +
    scale_y_continuous(limits = c(min(summary.param.set$Parameter)*0.95, max(summary.param.set$Parameter)*1.05)) +
    annotate("text", x = min(summary.param.set$Date), y = max(summary.param.set$Parameter)*1.04, size = font.size-7, label = paste(title[iteration])) +
    ylab(expression(k~"(g C "*g^"-1"*" C "*d^"-1"*")")) +
    theme_bw() +
    theme(legend.position = c(0.85,0.85)) +
    # theme(plot.title = element_text(size = font.size)) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0.25, 0.5), units="line")) +
    theme(text = element_text(size=font.size)) +
    theme(legend.key.height=unit(0.75,"line")) +
    # theme(legend.key.width=unit(3,"line")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# This sript plots leaf biomass pools (Figure 6G) for various test cases 
# with parameter shifted from potted seedling to free seedling
plot.Mleaf <- function(shift.output.Mleaf) { 
  ggplot() +
    geom_line(data = shift.output.Mleaf, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
    geom_point(size=2) +
    ylab(expression("Foliage mass"~"(g C "*plant^"-1"*")")) + 
    annotate("text", x = min(shift.output.Mleaf$Date), y = max(shift.output.Mleaf$value), size = font.size-7, label = paste(title[7])) +
    scale_colour_manual(breaks=c("0","1","2","3","4","5","6"), labels=c("Baseline (5L)","+ Cday of FS",expression(+ R[d]~"of FS"),
                                                                        expression(+ (a[f] + a[w] + a[r])~"of FS"),"+ Y of FS",expression(+ s[f]~"of FS"),"+ k of FS (Complete FS)"), 
                        values=cbPalette) +
    theme_bw() +
    theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size+2)) +
    theme(legend.key.height=unit(1,"line")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# This sript plots stem biomass pools (Figure 6H) for various test cases 
# with parameter shifted from potted seedling to free seedling
plot.Mstem <- function(shift.output.Mstem) { 
  ggplot() +
    geom_line(data = shift.output.Mstem, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
    geom_point(size=2) +
    ylab(expression("Wood mass"~"(g C "*plant^"-1"*")")) + 
    scale_colour_manual(breaks=c("0","1","2","3","4","5","6"), labels=c("Baseline (5L)","+ Cday of FS",expression(+ R[d]~"of FS"),
                                                                        expression(+ (a[f] + a[w] + a[r])~"of FS"),"+ Y of FS",expression(+ s[f]~"of FS"),"+ k of FS (Complete FS)"), 
                        values=cbPalette) +
    annotate("text", x = min(shift.output.Mstem$Date), y = max(shift.output.Mstem$value), size = font.size-7, label = paste(title[8])) +
    theme_bw() +
    theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size+2)) +
    theme(legend.key.height=unit(1,"line")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# This sript plots root biomass pools (Figure 6I) for various test cases 
# with parameter shifted from potted seedling to free seedling
plot.Mroot <- function(shift.output.Mroot) { 
  ggplot() +
    geom_line(data = shift.output.Mroot, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
    geom_point(size=2) +
    ylab(expression("Root mass"~"(g C "*plant^"-1"*")")) + 
    scale_colour_manual(breaks=c("0","1","2","3","4","5","6"), labels=c("Baseline (5L)","+ Cday of FS",expression(+ R[d]~"of FS"),
                                                                        expression(+ (a[f] + a[w] + a[r])~"of FS"),"+ Y of FS",expression(+ s[f]~"of FS"),"+ k of FS (Complete FS)"), 
                        values=cbPalette) +
    annotate("text", x = min(shift.output.Mroot$Date), y = max(shift.output.Mroot$value), size = font.size-7, label = paste(title[9])) +
    theme_bw() +
    theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size+2)) +
    theme(legend.key.height=unit(1,"line")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


