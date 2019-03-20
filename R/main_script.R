#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-- Script to analysis the raw WTC3 experiment data. 
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

###### R script to import and process the raw WTC3 experiment data 
# to model the carbon pools and fluxes using DA
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Setup

#- If you have cloned the github repository as suggested here https://github.com/kashifmahmud/DA_WTC3 ,
#      then this script should simply work. If you have downloaded a zipfile or otherwise copied the 
#      data and code manually, then you will need to set the working directory to wherever you've put the main directory.
#      The following line is an example- you will need to uncomment this and edit it to match your local machine.
# setwd("C:/Repos/DA_WTC3")


#- Load required libraries. There are quite a few (>20), including some non-standard functions that
#    are not on CRAN. This script will check for required libraries and install any that are missing.
source("R/loadLibraries.R")

# #- Download an the data. This downloads the zipfile from figshare
# download.file("https://ndownloader.figshare.com/files/4857112", "data.zip", mode="wb")
# # Extract data to a folder named "data".
# unzip("data.zip",overwrite=F)

#- Updating packages is highly recommended
#update.packages(ask=F)

#- load the custom analysis and plotting functions that do all of the actual work
source("R/functions.R")

#- export flag. Set to "T" to create pdfs of figures in "output/", or "F" to suppress output.
#- This flag is passed to many of the plotting functions below.
export=T


#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- plot R/A with a leaf-scale model (Figure 1)
plotCUE_conceptual_fig(toexport=export,Tdew=10,Ca=400,Vcmax=100,Jmax=150,Tleaf=20:42,PPFD=1500)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------






#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Plot respiration vs. temperature curves for leaves and whole-tree chambers (Figure 2)

#- get and process the whole-canopy R vs. T datasets
fits.list <- return_Rcanopy_closed()
fits.mass <- fits.list[[1]]     #- tree-level data
fits.trt <- fits.list[[2]]      #- treatment averages

#- plot R vs. T (Figure 2)
plotRvsT_figure2(fits.mass=fits.mass,fits.trt=fits.trt,export=export)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------







#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- plot Rbranch and Rleaf measured at 15 degrees C (Figure 3)
plotRleafRbranch(export=export)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Read and process the hourly flux dataset. Plot Figs 4-5 and Figs S2 and S6.

#- read in the hourly flux data
#dat.hr <- read.csv("data/WTC_TEMP_CM_WTCFLUX_20130910-20140530_L2_V1.csv")
dat.hr <- read.csv("data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
dat.hr$DateTime <- as.POSIXct(dat.hr$DateTime,format="%Y-%m-%d %H:%M:%S",tz="GMT")
dat.hr$Date <- as.Date(dat.hr$DateTime)

#- partition the net fluxes into GPP and Ra components. Specify the actitation energy and the number
#-  of prior nights to include for the estimate of the basal respiration rate (lagdates)
dat.hr.p <- partitionHourlyFluxCUE_arr(dat.hr.gf=dat.hr,Ea=57.69,lagdates=3)

#- plot an example week of partitioned fluxes (Figure S2)
plotPartitionedFluxes(dat.hr.gf3=dat.hr.p,ch_toplot="C07",startDate="2014-3-22",endDate="2014-3-27",export=export)

#- get daily sums from the partitioned hourly data
cue.list <- returnCUE.day(dat=dat.hr.p) # get daily sums from hourly data
cue.day <- cue.list[[1]]                # extract chamber values on each day
cue.day.trt <- cue.list[[2]]            # extract treatment averages

#- plot met, Ra, GPP, and Ra/GPP data over time (Figure 4)
plotPAR_AirT_CUE_GPP_Ra(cue.day.trt=cue.day.trt,export=export,lwidth=2.75)

#- plot PAR and Temperaure dependence of GPP, Ra, and Ra/GPP (Figure 5)
plotGPP_Ra_CUE_metdrivers(cue.day=cue.day,export=export,shading=0.7)

#- plot VPD and Tair dependence (Figure S6)
plotVPD_Tair(dat=dat.hr,export=export)


#- explore simpler metrics including net C flux, create Figs. S7-S8
#- calculate leaf-area specific rates of GPP, convert to umol CO2 m-2 s-1
cue.day$netC <- with(cue.day,Cgain+Closs) #g C day-1
cue.day$netC_la <- with(cue.day,netC/leafArea)
cue.day$Closs_Cgain <- with(cue.day,(-1*Closs)/Cgain)
plotNetC_metdrivers(cue.day=cue.day,export=export,shading=0.7)
plotClosstogainratio_metdrivers(cue.day=cue.day,export=export,shading=0.7)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------






#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Plot the temperature and PAR dependence of GPP per unit leaf area (Figure 6).
#- Note that this produces four separate graphs (panels a, b, and c, plus the legend).
#   These panels were manually combined to create Figure 6.
plotGPP_hex(dat=dat.hr.p,export=export,shading=0.7)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Plot the 5 diurnal observations of leaf-level photosynthesis and stomatal conductance (Figure 7).
#    Set printANOVAs to "T" to print ANOVAs for net photosynthesis (Anet) and stomatal conductance (Cond) on each date.
plotAnet_met_diurnals(export=export,lsize=2,printANOVAs=F)

#- compare the direct diurnal measurements with the whole-tree fluxes for a specified focal date. (Figure S9)
#-   Take some care here, as the flux data were frequently perterbed by investigators going into chambers on these dates
plotDiurnalvsWTC(focaldate=as.Date("2013-12-4"),output=T)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Does inhibition of respiration in the light affect the results? (Figs S4-5)

#- re-partition CO2 fluxes into Ra and GPP assuming R is reduced by 30% in the light
dat.hr.p2 <- partitionHourlyFluxCUE_arr(dat.hr.gf=dat.hr,Ea=57.69,lagdates=3,leafRtoTotal = 1,leafRreduction=0.3)

#- get daily sums again for light-inibited fluxes
cue.list2 <- returnCUE.day(dat=dat.hr.p2) # get daily sums from hourly data
cue.day2 <- cue.list2[[1]]                # extract chamber values on each day
cue.day.trt2 <- cue.list2[[2]]            # extract treatment averages

#- make plots comparing partitioning methods. (Figures S4-5)
plotCUE_paritioning_method(cue.day,cue.day2,export=export)
plotCUEvsT_partitioning_method(cue.day,cue.day2,export=export,shading=0.5)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Does tonight's repiration depend on yesterday's photosynthesis? (Figure S10)
Ra_GPP_coupling(cue.day=cue.day,export=export,shading=0.7)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Calculate climate metrics for methods

#- difference in average daily air temperature
#-- estimate the treatment effect on average daily air temperature

#baddays <- as.Date(c("2014-03-11","2013-12-30","2014-03-15","2014-03-16","2014-03-17")) #these few bad days increase the SD, but that is okay.
baddays <- c()
Tdat <- cue.day[,c("Date","T_treatment","Tair_24hrs")]
Tdat <- Tdat[with(Tdat,order(Date,T_treatment)),]
Tdat.m <- summaryBy(Tair_24hrs~Date+T_treatment,data=Tdat[!(Tdat$Date %in% baddays),] ,FUN=mean)

Tdat.a <- subset(Tdat.m,T_treatment=="ambient")
Tdat.e <- subset(Tdat.m,T_treatment=="elevated")
Tdat.both <- merge(Tdat.a,Tdat.e,by=c("Date"))
Tdat.both$diff <- with(Tdat.both,Tair_24hrs.mean.y-Tair_24hrs.mean.x)
mean(Tdat.both$diff);sd(Tdat.both$diff)


#- average CO2 concentration during sunlit hours
daytime <- subset(dat.hr.p,PAR>100)
hist(daytime$CO2centralCh,xlim=c(350,600),breaks=101)
mean(daytime$CO2centralCh,na.rm=T)
median(daytime$CO2centralCh,na.rm=T)
summaryBy(CO2centralCh~T_treatment,data=daytime,FUN=c(mean,sd),na.rm=T)




#- get relative humidity
daytime$RH <- VPDtoRH(VPD=daytime$VPDair,TdegC=daytime$Tair_al)
summaryBy(RH~T_treatment,data=daytime,FUN=c(mean,sd))
summaryBy(VPDair~T_treatment,data=daytime,FUN=c(mean,sd))


#- get water potentials
wp <- read.csv("data/WTC_TEMP_CM_WATERPOTENTIAL-PREDAWN-MIDDAY_20130515-20140424_L2.csv")
wp$date <- as.Date(wp$date,format="%d/%m/%Y")
wp.ch <- summaryBy(predawn+midday~T_treatment+chamber,data=subset(wp,Water_treatment=="control"),na.rm=T,keep.names=T)
wp.t <- summaryBy(predawn+midday~T_treatment,data=wp.ch,na.rm=T,FUN=c(mean,standard.error))

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Process tree Diameter, height, and leaf area.
plot_tree_size(export=F)
table_tree_size()
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
