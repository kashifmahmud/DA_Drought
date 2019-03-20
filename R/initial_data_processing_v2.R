#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-- Script to analysis the raw WTC3 experiment data. 
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
###### R script to import and process the raw WTC3 experiment data 
# to model the carbon pools and fluxes using DA

# # Clear the workspace (if needed)
# rm(list=ls())
# 
# #- Load required libraries. There are quite a few (>20), including some non-standard functions that
# #    are not on CRAN. This script will check for required libraries and install any that are missing.
# source("R/load_packages_WTC3.R")
# source("R/loadLibraries_WTC3.R")
# 
# #- load the custom analysis and plotting functions that do all of the actual work
# source("R/functions_WTC3.R")

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

####### Based on the height and dia data, we decided to go with Case 3 that represents drought and watered treatments seperately from the beginning of the experiment 
####### so total number of treatment scenario = 4 (ambient+watered, ambient+droughted, elevated+watered, elevated+droughted) with sample size of n = 3 for each option

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
c1 = 0.48 # unit conversion from gDM to gC: 1 gDM = 0.48 gC

#- Script to read and setup a model for GPP, Ra (aboveground respiration) and leaf area
# read the hourly flux data
data.hr <- read.csv("raw_data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
data.hr$DateTime <- as.POSIXct(data.hr$DateTime,format="%Y-%m-%d %H:%M:%S",tz="GMT")
data.hr$Date <- as.Date(data.hr$DateTime)

#- partition the net fluxes into GPP and Ra components. Specify the actitation energy and the number
#-  of prior nights to include for the estimate of the basal respiration rate (lagdates)
data.hr.p <- partitionHourlyFluxCUE_arr(data.hr.gf=data.hr,Ea=57.69,lagdates=3)

#- get daily sums from the partitioned hourly data
cue.list <- returnCUE.day(dat=data.hr.p) # get daily sums from hourly data
cue.day <- cue.list[[1]]                # extract chamber values on each day
cue.day.trt <- cue.list[[2]]            # extract treatment averages

data=cue.day.trt
names(data)[4:9] = c("GPP","Ra","LA","GPP_SE","Ra_SE","LA_SE")
# ignore the first few days of data to start on 2013-09-17 from where we have H & D measurements
data = subset(data, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-05-27"))

# write csv file with daily inputs of GPP, Ra, LA
write.csv(data, file = "processed_raw_data/data_GPP_Ra_LA.csv", row.names = FALSE)

#- plot GPP, Ra, and LA data over time for various treatments
i = 0
font.size = 12
plot = list()
pd <- position_dodge(5)
plots[[1]] = ggplot(data=data, aes(x=Date, y=GPP, group = interaction(T_treatment,chamber_type), colour=T_treatment, shape=chamber_type)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=GPP-GPP_SE, ymax=GPP+GPP_SE), colour="grey", width=1) +
  geom_line(position=pd, data = data, aes(x = Date, y = GPP, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) +
  ylab(expression(GPP~"(g C "*d^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data$Date)-2, max(data$Date)+2)) +
  labs(colour="Temperature", shape="Chamber type", linetype="Chamber type") +
  scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.25,0.9), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

plots[[2]] = ggplot(data=data, aes(x=Date, y=Ra, group = interaction(T_treatment,chamber_type), colour=T_treatment, shape=chamber_type)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=Ra-Ra_SE, ymax=Ra+Ra_SE), colour="grey", width=1) +
  geom_line(position=pd, data = data, aes(x = Date, y = Ra, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) +
  ylab(expression(R[a]~"(g C "*d^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data$Date)-2, max(data$Date)+2)) +
  labs(colour="Temperature", shape="Chamber type", linetype="Chamber type") +
  scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.25,0.9), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

pdf(file = "output/1.GPP_Ra_over_time.pdf", width=20, height=5)
plots[[1]]
plots[[2]]
dev.off()


# #-----------------------------------------------------------------------------------------
# #- plot GPP, Ra, and LA data over time for various treatments
# i = 0
# font.size = 12
# plot = list()
# meas = as.factor(c("GPP","Ra","LA"))
# error = as.factor(c("GPP_SE","Ra_SE","LA_SE"))
# title = as.character(c("A","B","C"))
# pd <- position_dodge(0.5)
# 
# # data.sub <- summaryBy(GPP+Ra+LA ~ Date+T_treatment, data=data, FUN=c(mean,standard.error,length))
# data.sub <- subset(data, chamber_type %in% as.factor("watered"))
# # names(data.sub)[4:9] = c("GPP","Ra","LA","GPP_SE","Ra_SE","LA_SE")
# melted.data = melt(data.sub[,c("GPP","Ra","LA","Date","T_treatment")], id.vars=c("Date","T_treatment"))
# melted.data.error = melt(data.sub[,c("GPP_SE","Ra_SE","LA_SE","Date","T_treatment")], id.vars=c("Date","T_treatment"))
# melted.data$SE = melted.data.error$value
# 
# for (p in 1:length(meas)) {
#   melted.data.sub = subset(melted.data,variable %in% meas[p])
#   i = i + 1
#   plot[[i]] = ggplot(melted.data.sub, aes(x=Date, y=value, group = T_treatment, colour=T_treatment)) + 
#     geom_ribbon(data = melted.data.sub, aes(ymin=value-SE, ymax=value+SE), linetype=2, alpha=0.1,size=0.1) +
#     geom_line(position=pd,data = melted.data.sub, aes(x = Date, y = value, group = T_treatment, colour=T_treatment)) + 
#     scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
#     labs(colour="Temperature") +
#     scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
#     theme_bw() +
#     theme(legend.title = element_text(colour="black", size=font.size)) +
#     theme(legend.text = element_text(colour="black", size = font.size)) +
#     theme(legend.position = c(0.15,0.75)) +
#     theme(legend.key = element_blank()) +
#     theme(text = element_text(size=font.size)) +
#     theme(axis.title.x = element_blank()) +
#     theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
#   
#   if (p==1) {
#     plot[[i]] = plot[[i]] + ylab(expression(GPP~"(g C "*d^"-1"*")"))
#   } else if (p==2) {
#     plot[[i]] = plot[[i]] + ylab(expression(R[a]~"(g C "*d^"-1"*")"))
#   } else if (p==3) {
#     plot[[i]] = plot[[i]] + ylab(expression("Leaf area"~"("*m^"2"*" "*plant^"-1"*")"))
#   }
#   if (p>1) {
#     plot[[i]] = plot[[i]] + guides(colour=FALSE)
#   }
# }
# 
# # png("output/1_1.GPP_Ra_LA.png", units="px", width=2200, height=1600, res=220)
# pdf(file = "output/1.1.GPP_Ra_LA_over_time.pdf")
# print (do.call(grid.arrange,  plot))
# 
# # data.sub <- summaryBy(GPP+Ra+LA ~ Date+chamber_type, data=data, FUN=c(mean,standard.error))
# # names(data.sub)[3:8] = c("GPP","Ra","LA","GPP_SE","Ra_SE","LA_SE")
# data.sub <- subset(data, T_treatment %in% as.factor("ambient"))
# melted.data = melt(data.sub[,c("GPP","Ra","LA","Date","chamber_type")], id.vars=c("Date","chamber_type"))
# melted.data.error = melt(data.sub[,c("GPP_SE","Ra_SE","LA_SE","Date","chamber_type")], id.vars=c("Date","chamber_type"))
# melted.data$SE = melted.data.error$value
# 
# i = 0
# plot = list()
# for (p in 1:length(meas)) {
#   melted.data.sub = subset(melted.data,variable %in% meas[p])
#   i = i + 1
#   plot[[i]] = ggplot(melted.data.sub, aes(x=Date, y=value, group = chamber_type, colour=chamber_type)) + 
#     geom_ribbon(data = melted.data.sub, aes(ymin=value-SE, ymax=value+SE), linetype=2, alpha=0.1,size=0.1) +
#     geom_line(position=pd,data = melted.data.sub, aes(x = Date, y = value, group = chamber_type, colour=chamber_type)) + 
#     scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
#     labs(colour="Treatment") +
#     scale_color_manual(labels = c("drought", "watered"), values = c("red", "blue")) +
#     theme_bw() +
#     theme(legend.title = element_text(colour="black", size=font.size)) +
#     theme(legend.text = element_text(colour="black", size = font.size)) +
#     theme(legend.position = c(0.15,0.75)) +
#     theme(legend.key = element_blank()) +
#     theme(text = element_text(size=font.size)) +
#     theme(axis.title.x = element_blank()) +
#     theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
#   
#   if (p==1) {
#     plot[[i]] = plot[[i]] + ylab(expression(GPP~"(g C "*d^"-1"*")"))
#   } else if (p==2) {
#     plot[[i]] = plot[[i]] + ylab(expression(R[a]~"(g C "*d^"-1"*")"))
#   } else if (p==3) {
#     plot[[i]] = plot[[i]] + ylab(expression("Leaf area"~"("*m^"2"*" "*plant^"-1"*")"))
#   }
#   if (p>1) {
#     plot[[i]] = plot[[i]] + guides(colour=FALSE)
#   }
# }
# 
# print (do.call(grid.arrange,  plot))
# dev.off()
# #-----------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Script to read and plot stem height and diameter for various treatment cases 
# Read data form WTC3 experiment
height.dia.raw <- read.csv("raw_data/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV") # units: Height = cm, dia = mm
height.dia.raw$DateTime = as.Date(height.dia.raw$DateTime)
flux <- read.csv("raw_data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
chambers = unique(flux$chamber)

height.dia = data.frame(chamber = rep(chambers, each = length(unique(height.dia.raw$DateTime))),
                        Date = rep(unique(height.dia.raw$DateTime),length(chambers)),
                        T_treatment = character(length(chambers) * length(unique(height.dia.raw$DateTime))),
                        W_treatment = character(length(chambers) * length(unique(height.dia.raw$DateTime))),
                        diameter = numeric(length(chambers) * length(unique(height.dia.raw$DateTime))),
                        height = numeric(length(chambers) * length(unique(height.dia.raw$DateTime))), stringsAsFactors=FALSE)
height.dia = subset(height.dia, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-05-27"))

for(i in 1:length(chambers)) {
  height.dia.sub = subset(height.dia.raw, chamber %in% as.factor(chambers[i]))
  # height.dia.sub = subset(height.dia, chamber %in% as.factor(chambers[i]) & Water_treatment %in% as.factor("control"))
  if (i==11) {  # Remove the reference measurements made on Chamber 11
    height.dia.sub.1 = subset(height.dia.sub, Stem_number %in% 1)
    height.dia.sub.1.2 = subset(height.dia.sub, Stem_number %in% 1.2)
    height.dia.sub.1[height.dia.sub.1$DateTime >= as.Date("2013-12-24"),"Plant_height"] = height.dia.sub.1.2[,"Plant_height"]
    height.dia.sub = height.dia.sub.1
  }
  keeps <- c("chamber", "DateTime", "T_treatment", "Water_treatment", "Plant_height", "X15", "X65")
  height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
  
  D.15 <- lm(X15 ~ X65, data=height.dia.sub)
  visreg(D.15, "X65", overlay=TRUE)
  summary(D.15)
  eq.D = function(x){coefficients(D.15)[1] + coefficients(D.15)[2] * x }
  
  height.dia.sub$X15 = eq.D(height.dia.sub$X65)
  height.dia.sub = subset(height.dia.sub, DateTime >= as.Date("2013-09-14") & DateTime <= as.Date("2014-05-27"))
  height.dia.sub = height.dia.sub[!is.na(height.dia.sub$X65),]
  keeps <- c("chamber", "DateTime", "T_treatment", "Water_treatment", "X15", "Plant_height")
  height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
  names(height.dia.sub) <- c("chamber","Date","T_treatment","W_treatment","diameter","height")
  height.dia.sub$T_treatment = as.character(height.dia.sub$T_treatment)
  height.dia.sub$W_treatment = as.character(height.dia.sub$W_treatment)
  
  height.dia[(1+(i-1)*length(unique(height.dia$Date))) : (i*length(unique(height.dia$Date))), 
             c("T_treatment","W_treatment","diameter","height")] = height.dia.sub[,c("T_treatment","W_treatment","diameter","height")]
}
height.dia$T_treatment = as.factor(height.dia$T_treatment)
height.dia$W_treatment = as.factor(height.dia$W_treatment)

# Average the ambient and elevated temperature treatments considering the drought/watered treatment seperated from the start of the experiment
# n=3 for whole period
drought.chamb = unique(height.dia$chamber[ height.dia$W_treatment %in% as.factor("drydown")])
height.dia$chamber_type = as.factor( ifelse(height.dia$chamber %in% drought.chamb, "drought", "watered") )
height.dia.final <- summaryBy(height+diameter ~ Date+T_treatment+chamber_type, data=height.dia, FUN=c(mean,standard.error))
names(height.dia.final)[4:7] = c("height", "diameter", "height_SE", "diameter_SE")

#----------------------------------------------------------------------------------------------------------------
# #- plot H and D data over time for various treatments
# i = 0
# font.size = 12
# plot = list()
# meas = as.factor(c("diameter","height"))
# error = as.factor(c("diameter_SE","height_SE"))
# title = as.character(c("A","B","C"))
# pd <- position_dodge(0)
# 
# # height.dia.sub <- summaryBy(height+diameter ~ Date+T_treatment, data=height.dia.final, FUN=c(mean,standard.error))
# # names(height.dia.sub)[3:6] = c("height","diameter","height_SE","diameter_SE")
# height.dia.sub <- subset(height.dia.final, chamber_type %in% as.factor("watered"))
# melted.height.dia = melt(height.dia.sub[,c("height","diameter","Date","T_treatment")], id.vars=c("Date","T_treatment"))
# melted.height.dia.error = melt(height.dia.sub[,c("height_SE","diameter_SE","Date","T_treatment")], id.vars=c("Date","T_treatment"))
# melted.height.dia$SE = melted.height.dia.error$value
# 
# for (p in 1:length(meas)) {
#   melted.height.dia.sub = subset(melted.height.dia,variable %in% meas[p])
#   i = i + 1
#   plot[[i]] = ggplot(melted.height.dia.sub, aes(x=Date, y=value, group = T_treatment, colour=T_treatment)) + 
#     geom_ribbon(data = melted.height.dia.sub, aes(ymin=value-SE, ymax=value+SE), linetype=2, alpha=0.1,size=0.1) +
#     geom_line(position=pd,data = melted.height.dia.sub, aes(x = Date, y = value, group = T_treatment, colour=T_treatment)) + 
#     ylab(paste(as.character(meas[p]),"(mm)")) +
#     xlab("") +
#     scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
#     labs(colour="Temperature") +
#     scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
#     theme_bw() +
#     theme(legend.title = element_text(colour="black", size=font.size)) +
#     theme(legend.text = element_text(colour="black", size = font.size)) +
#     theme(legend.position = c(0.15,0.75)) +
#     theme(legend.key = element_blank()) +
#     theme(text = element_text(size=font.size)) +
#     theme(axis.title.x = element_blank()) +
#     theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
#   
#   if (p==2) {
#     plot[[i]] = plot[[i]] + ylab(paste(as.character(meas[p]),"(cm)"))
#   }
#   if (p>1) {
#     plot[[i]] = plot[[i]] + guides(colour=FALSE)
#   }
# }
# 
# pdf(file = "output/2.H_D.pdf")
# print (do.call(grid.arrange,  plot))
# 
# # height.dia.sub <- summaryBy(height+diameter ~ Date+chamber_type, data=height.dia.final, FUN=c(mean,standard.error))
# # names(height.dia.sub)[3:6] = c("height","diameter","height_SE","diameter_SE")
# height.dia.sub <- subset(height.dia.final, T_treatment %in% as.factor("ambient"))
# melted.height.dia = melt(height.dia.sub[,c("height","diameter","Date","chamber_type")], id.vars=c("Date","chamber_type"))
# melted.height.dia.error = melt(height.dia.sub[,c("height_SE","diameter_SE","Date","chamber_type")], id.vars=c("Date","chamber_type"))
# melted.height.dia$SE = melted.height.dia.error$value
# 
# i = 0
# plot = list()
# for (p in 1:length(meas)) {
#   melted.height.dia.sub = subset(melted.height.dia,variable %in% meas[p])
#   i = i + 1
#   plot[[i]] = ggplot(melted.height.dia.sub, aes(x=Date, y=value, group = chamber_type, colour=chamber_type)) + 
#     geom_ribbon(data = melted.height.dia.sub, aes(ymin=value-SE, ymax=value+SE), linetype=2, alpha=0.1,size=0.1) +
#     geom_line(position=pd,data = melted.height.dia.sub, aes(x = Date, y = value, group = chamber_type, colour=chamber_type)) + 
#     ylab(paste(as.character(meas[p]),"(mm)")) +
#     xlab("") +
#     scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
#     labs(colour="Treatment") +
#     scale_color_manual(labels = c("drought", "watered"), values = c("red", "blue")) +
#     theme_bw() +
#     theme(legend.title = element_text(colour="black", size=font.size)) +
#     theme(legend.text = element_text(colour="black", size = font.size)) +
#     theme(legend.position = c(0.15,0.75)) +
#     theme(legend.key = element_blank()) +
#     theme(text = element_text(size=font.size)) +
#     theme(axis.title.x = element_blank()) +
#     theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
#   
#   if (p==2) {
#     plot[[i]] = plot[[i]] + ylab(paste(as.character(meas[p]),"(cm)"))
#   }
#   if (p>1) {
#     plot[[i]] = plot[[i]] + guides(colour=FALSE)
#   }
# }
# 
# print (do.call(grid.arrange,  plot))
# dev.off()
#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Script to read and setup a model for the root mass from sink limited pot experiment and WTC3 harvest data to predict 
# initial root mass from Height and Diameter
plot.summary.pot = read.csv("data_pot_experiment/plot_summary.csv")
initial.data.pot <- read.csv("data_pot_experiment/seedling_initial.csv")
keeps <- c("root_mass", "diameter_15", "height")
data.pot = initial.data.pot[ , keeps, drop = FALSE]
names(data.pot) <- c("rootmass","diameter","height")
data.pot$volume = as.integer(1000)

harvest.data.pot <- read.csv("data_pot_experiment/seedling_mass.csv")
harvest.data.pot$rootmass = harvest.data.pot$coarseroot + harvest.data.pot$fineroot
keeps <- c("rootmass", "plot", "pot", "volume")
harvest.data.pot = harvest.data.pot[ , keeps, drop = FALSE]

height.dia.pot <- read.csv("data_pot_experiment/height_diameter.csv")
height.dia.pot$Date <- parse_date_time(height.dia.pot$Date,"d m y")
height.dia.pot$Date = as.Date(height.dia.pot$Date, format = "%d/%m/%Y")
height.dia.pot = subset(height.dia.pot,Date=="2013-05-21")

height.dia.pot <- merge(plot.summary.pot, height.dia.pot, by=c("pot","plot"))
harvest.data.pot <- merge(height.dia.pot, harvest.data.pot, by=c("pot","plot","volume"))
keeps <- c("rootmass", "diameter", "height","volume")
harvest.data.pot = harvest.data.pot[ , keeps, drop = FALSE]
data.pot = rbind(data.pot, harvest.data.pot)
# data.pot[,"rootmass"] = data.pot[,"rootmass"] * c1 # unit conversion: gDM to gC 
data.free.pot = subset(data.pot, volume %in% 1000)
data.free.pot = data.free.pot[,-4]

# Read harvest rootmass data form WTC3 experiment
rootmass.harvest = read.csv("raw_data/WTC_TEMP_CM_HARVEST-ROOTS_20140529-20140606_L1_v1.csv")
# rootmass.harvest$chamber_type = as.factor( ifelse(rootmass.harvest$chamber %in% drought.chamb, "drought", "watered") )
rootmass.harvest$Date = as.Date("2014-05-26")
rootmass.harvest = merge(unique(height.dia[,c("chamber","T_treatment","chamber_type")]), rootmass.harvest, by=c("chamber"))

height.dia.sub = subset(height.dia, Date %in% as.Date("2014-05-27") & T_treatment %in% as.factor("ambient") & chamber_type %in% as.factor("watered"))
keeps <- c("chamber", "height", "diameter")
height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
height.dia.sub <- merge(height.dia.sub, rootmass.harvest[,c("chamber","RootDMtotal","chamber_type")], by=c("chamber"))
keeps <- c("RootDMtotal", "diameter", "height")
height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
names(height.dia.sub) <- c("rootmass","diameter","height")

data.free.pot = rbind(data.free.pot, height.dia.sub) # Merge pot exp data (initial and harvest of free seedlings) and WTC3 harvest data (ambient)

# Fit a linear regression to estimate root mass by stem dia and height (ignoring temperature variation)
# the regression is fitted using the Sink limited pot experiemnt of Court Campany for similar seedlings (Euc.Teri.)
rm1 <- lm(log(rootmass) ~ log(diameter) + log(height), data=data.free.pot)
rm2 <- lm(log(rootmass) ~ log(diameter) * log(height), data=data.free.pot)

# png("output/2.model_comparison.png", units="px", width=2200, height=1600, res=220)
# layout(matrix(c(1,2,3,4),2,2,byrow=TRUE), widths=c(1,1), heights=c(10,10))
# par(oma = c(0, 0, 0, 0))
# visreg(rm1, "diameter")
# visreg(rm1, "height")
# visreg(rm2, "diameter", by="height", overlay=TRUE)
# visreg(rm2, "height", by="diameter", overlay=TRUE)
# dev.off()
# 
# sink("Output/2.model_comparison.txt")
# cat("\n\nRootmass models:\n----------------\n### Linear regression ignoring temperature variation:");  summary(rm1)
# cat("\n### Linear regression considering interaction with temperature:"); summary(rm2)
# cat("### Comparison between both models:\n")
# AIC(rm1, rm2); BIC(rm1, rm2)
# sink()

# Estimate the WTC3 root mass from the fitted linear regression equation
# Go for the regression with interaction (rm1) that suits the pot experiment (Euc.Teri.) data best
# otherwise we could only use diameter to predict the rootmass (as intercept and height aren't significant)
# eq.rm = function(x,y){exp(coefficients(rm1)[1] + coefficients(rm1)[2] * log(x)  + coefficients(rm1)[3] * log(y))}
eq.rm = function(x,y){exp(coefficients(rm2)[1] + coefficients(rm2)[2] * log(x)  + coefficients(rm2)[3] * log(y)
                          + coefficients(rm2)[4] * log(x) * log(y) )}

# estimate the initial rootmass from H and D data
height.dia.initial = subset(height.dia.final, Date %in% as.Date("2013-09-17"))
height.dia.initial$RM = eq.rm(height.dia.initial$diameter, height.dia.initial$height)
# height.dia.initial$RM_SE = (((height.dia.initial$diameter_SE/height.dia.initial$diameter)^2 + (height.dia.initial$height_SE/height.dia.initial$height)^2 )^0.5) * height.dia.initial$RM
height.dia.initial$RM_SE = ( (((coefficients(rm2)[2]*height.dia.initial$diameter_SE/height.dia.initial$diameter)^2 + (coefficients(rm2)[3]*height.dia.initial$height_SE/height.dia.initial$height)^2 +
                                 coefficients(rm2)[4]*((height.dia.initial$diameter_SE/height.dia.initial$diameter)^2 + (height.dia.initial$height_SE/height.dia.initial$height)^2) )^0.5) ) * height.dia.initial$RM


# processing the harvest rootmass
rootmass.harvest.mean <- summaryBy(RootDMtotal ~ Date+T_treatment+chamber_type, data=rootmass.harvest, FUN=c(mean,standard.error))
names(rootmass.harvest.mean)[4:5] = c("RM","RM_SE")

# # initialize the biomass data frame
# Dates = rep(seq(as.Date("2013-09-17"), as.Date("2014-05-27"), by="days"), 4)
# biomass = data.frame(Date = Dates,
#                      T_treatment = rep(unique(rootmass.harvest$T_treatment), each=length(Dates)/2))
# biomass = biomass[with(biomass, order(T_treatment,Date)), ]
# biomass$chamber_type = rep(unique(rootmass.harvest$chamber_type), length(Dates)/2)
# 
# # combine the rootmass data with the biomass
# rootmass = merge(height.dia.initial[,c("Date","T_treatment","chamber_type","RM","RM_SE")], rootmass.harvest, all = TRUE)
# biomass = merge(biomass, rootmass, all = TRUE)
# biomass$RM = biomass$RM * c1 # unit conversion from gDM to gC
# biomass$RM_SE = biomass$RM_SE * c1 # unit conversion from gDM to gC

rootmass = merge(height.dia.initial[,c("Date","T_treatment","chamber_type","RM","RM_SE")], rootmass.harvest.mean, all = TRUE)
rootmass$RM = rootmass$RM * c1 # unit conversion from gDM to gC
rootmass$RM_SE = rootmass$RM_SE * c1 # unit conversion from gDM to gC
# data.with.RM = merge(data, rootmass, all = TRUE)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


# #--------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------
# # Script to estimate leaf mass
# harvest.wtc = read.csv("raw_data/WTC_TEMP_CM_HARVEST-CANOPY_20140526-20140528_L1_v1.csv")
# keeps <- c("chamber", "T_treatment", "SLA", "BranchDM", "StemDM", "TotalLeafDM")
# harvest.wtc = harvest.wtc[ , keeps, drop = FALSE]
# 
# # calculate canopy SLA taking the weighted average of 3 laters of canopy (low, mid, top)
# leafmass <- summaryBy(TotalLeafDM~chamber,data=harvest.wtc,FUN=sum,keep.names=F)
# harvest.wtc <- merge(harvest.wtc,leafmass,by="chamber")
# harvest.wtc$weights <- with(harvest.wtc,TotalLeafDM/TotalLeafDM.sum)
# SLA.weighted <- plyr::ddply(harvest.wtc,~chamber,summarise, SLA.weighted=weighted.mean(SLA,weights))
# harvest.wtc <- merge(harvest.wtc,SLA.weighted,by="chamber")
# 
# harvest.wtc <- dplyr::summarize(group_by(harvest.wtc,chamber,T_treatment),
#                                 SLA = mean(SLA.weighted,na.rm=T),
#                                 BranchDM = sum(BranchDM,na.rm=T),
#                                 StemDM = sum(StemDM,na.rm=T),
#                                 TotalLeafDM=sum(TotalLeafDM,na.rm=T),
#                                 Date = max(melted.data$Date))
# harvest.wtc$LMA = 1 / (harvest.wtc$SLA / (100*100))
# harvest.wtc = data.frame(harvest.wtc)
# harvest.wtc$chamber_type = as.factor( ifelse(harvest.wtc$chamber %in% drought.chamb, "drought", "watered") )
# 
# harvest.wtc.mean <- summaryBy(SLA+BranchDM+StemDM+TotalLeafDM+LMA ~ Date+T_treatment+chamber_type, data=harvest.wtc, FUN=c(mean,standard.error,length))
# harvest.wtc.mean = harvest.wtc.mean[,-c((ncol(harvest.wtc.mean)-3) : ncol(harvest.wtc.mean))]
# names(harvest.wtc.mean) = c("Date", "T_treatment", "chamber_type","SLA", "BM", "SM", "LM", "LMA", "SLA_SE", "BM_SE", "SM_SE", "LM_SE", "LMA_SE", "sample.size")
# harvest.wtc.mean$experiment = as.factor("Harvest")
# 
# # read LMA values
# leaf.trait = read.csv("raw_data/WTC_TEMP_CM_TDLSUNSHADE-TRAIT_20131026-20140421.csv")
# names(leaf.trait)[4] = "ch"; names(leaf.trait)[6] = "T_treatment"
# leaf.trait = merge(data.frame(ch=unique(as.factor(leaf.trait$ch)), chamber= as.factor(harvest.wtc$chamber)), leaf.trait, by="ch", all=TRUE)
# leaf.trait$chamber_type = as.factor( ifelse(leaf.trait$chamber %in% drought.chamb, "drought", "watered") )
# 
# leaf.trait$Date = paste( leaf.trait$Year , leaf.trait$Month , 15, sep = "-" )
# leaf.trait$Date <- as.factor(leaf.trait$Date)
# leaf.trait$Date <- parse_date_time(leaf.trait$Date,"y m d")
# leaf.trait$Date = as.Date(leaf.trait$Date, format = "%Y-%m-%d")
# keeps = c("Date", "chamber", "T_treatment", "chamber_type", "leaf", "lma")
# leaf.trait = leaf.trait[ , keeps, drop = FALSE]
# 
# leaf.trait = leaf.trait[,-2]
# leaf.trait.sub <- summaryBy(lma~Date+T_treatment+chamber_type+leaf, data=leaf.trait, FUN=c(mean,standard.error,length))
# names(leaf.trait.sub)[5:7] = c("LMA","LMA_SE","sample.size")
# 
# plots = list() 
# pd <- position_dodge(3)
# plots[[1]] = ggplot(data = leaf.trait.sub, aes(x=Date, y=LMA, group=leaf, colour=as.factor(leaf))) + 
#   geom_point(position=pd) + 
#   geom_errorbar(position=pd,aes(ymin=LMA-LMA_SE, ymax=LMA+LMA_SE), colour="grey", width=2) +
#   geom_line(position=pd, data = leaf.trait.sub, aes(x=Date, y=LMA, group=leaf, colour=as.factor(leaf))) +
#   labs(colour="Treatments", shape="Leaf type", linetype="Leaf type") +
#   facet_wrap(T_treatment ~ chamber_type, scales="free_x") +
#   ylab(expression("LMA"~"(g "*m^"-2"*")")) + theme_bw() +
#   theme(legend.title = element_text(colour="black", size=font.size)) +
#   theme(legend.text = element_text(colour="black", size = font.size)) +
#   theme(legend.position = c(0.15,0.94)) + theme(legend.key.height=unit(0.8,"line")) +
#   theme(legend.key = element_blank()) +
#   theme(text = element_text(size=font.size)) +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
# 
# leaf.trait.sub1 <- summaryBy(lma~Date+T_treatment+chamber_type, data=leaf.trait, FUN=c(mean,standard.error,length))
# names(leaf.trait.sub1)[4:6] = c("LMA","LMA_SE","sample.size")
# leaf.trait.sub1$experiment = as.factor("Court")
# plots[[2]] = ggplot(data = leaf.trait.sub1, aes(x=Date, y=LMA, group=interaction(T_treatment,chamber_type), colour=as.factor(T_treatment), shape=as.factor(chamber_type))) + 
#   geom_point(position=pd) + 
#   geom_errorbar(position=pd,aes(ymin=LMA-LMA_SE, ymax=LMA+LMA_SE), colour="grey", width=2) +
#   geom_line(position=pd, data = leaf.trait.sub1, aes(x=Date, y=LMA, group=interaction(T_treatment,chamber_type), colour=as.factor(T_treatment), linetype=as.factor(chamber_type))) + 
#   labs(colour="Treatments", shape="Chamber type", linetype="Chamber type") +
#   ylab(expression("LMA"~"(g "*m^"-2"*")")) + theme_bw() +
#   theme(legend.title = element_text(colour="black", size=font.size)) +
#   theme(legend.text = element_text(colour="black", size = font.size)) +
#   theme(legend.position = c(0.5,0.9), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
#   theme(legend.key = element_blank()) +
#   theme(text = element_text(size=font.size)) +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# # png("output/3.LMA_over_time.png", units="px", width=2200, height=1600, res=220)
# pdf(file = "output/3.LMA_Court_over_time.pdf")
# # lay <- rbind(c(1,1),c(2,3))
# # grid.arrange(grobs = plots[c(5,1,2)], layout_matrix = lay)
# plots
# dev.off()
# 
# # Merge all LMA valuse from difernt experiments (Mike/Court/Harvest)
# lma.mike = read.csv("raw_data/WTC_TEMP_CM_GX-ASAT_20130515-20140402_L2.csv")
# keeps <- c("date","chamber","leaf","T_treatment","lma")
# lma.mike = lma.mike[ , keeps, drop = FALSE]
# lma.mike = unique(lma.mike)
# lma.mike$chamber_type = as.factor( ifelse(lma.mike$chamber %in% drought.chamb, "drought", "watered") )
# names(lma.mike)[1] = "Date"
# lma.mike$Date = as.Date(lma.mike$Date)
# lma.mike = na.omit(lma.mike)
# 
# lma.mike.sum <- summaryBy(lma ~ Date+T_treatment+chamber_type, data=lma.mike, FUN=c(mean,standard.error,length))
# names(lma.mike.sum)[4:6] = c("LMA","LMA_SE","sample.size")
# lma.mike.sum$experiment = as.factor("Mike")
# 
# # combine all 3 lma data sets
# lma.combined = rbind(lma.mike.sum, leaf.trait.sub1, harvest.wtc.mean[,c("Date","T_treatment","chamber_type","LMA","LMA_SE","sample.size","experiment")])
# 
# plots = list() 
# pd <- position_dodge(2)
# plots = ggplot(data = lma.combined, aes(x=Date, y=LMA, group=experiment, colour=as.factor(experiment))) + 
#   geom_point(position=pd) + 
#   geom_errorbar(position=pd,aes(ymin=LMA-LMA_SE, ymax=LMA+LMA_SE), colour="grey", width=2) +
#   geom_line(position=pd, data = lma.combined, aes(x=Date, y=LMA, group=experiment, colour=as.factor(experiment))) +
#   labs(colour="Treatments", shape="Leaf type", linetype="Leaf type") +
#   facet_wrap(T_treatment ~ chamber_type, scales="free_x") +
#   ylab(expression("LMA"~"(g "*m^"-2"*")")) + theme_bw() +
#   theme(legend.title = element_text(colour="black", size=font.size)) +
#   theme(legend.text = element_text(colour="black", size = font.size)) +
#   theme(legend.position = c(0.15,0.92)) + theme(legend.key.height=unit(0.8,"line")) +
#   theme(legend.key = element_blank()) +
#   theme(text = element_text(size=font.size)) +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
# 
# 
# # png("output/3.LMA_over_time.png", units="px", width=2200, height=1600, res=220)
# pdf(file = "output/3.LMA_combined_over_time.pdf")
# plots
# dev.off()
# 
# #----------------------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------
# 
# ########### how to estimate the LMA over time????????
# 
# # Looking at both figures, we don’t find any particular trend over time for LMA (whether sun/shade or any treatment effect) 
# # ignoring the last points (which could be considered as an artefact) and that’s why we have decided to take one mean LMA 
# # (without considering any time effect) for individual treatment group. 
# leaf.trait.final <- summaryBy(lma ~ T_treatment+chamber_type, data=leaf.trait, FUN=c(mean,standard.error))
# names(leaf.trait.final)[3:4] = c("LMA","LMA_SE")
# 
# #----------------------------------------------------------------------------------------------------------------
# # calculate daily leaf mass
# data.with.LM = merge(data.with.RM, leaf.trait.final, by=c("T_treatment", "chamber_type"), all=TRUE)
# data.with.LM$LM = (data.with.LM$LA * data.with.LM$LMA) 
# data.with.LM$LM = data.with.LM$LM * c1 # unit conversion from gDM to gC
# # data.with.LM$LM_SE = (((data.with.LM$LA_SE * (3)^0.5)^2 * (data.with.LM$LMA_SE * (3)^0.5)^2)/2  )^0.5
# data.with.LM$LM_SE = (((data.with.LM$LA_SE/data.with.LM$LA)^2 + (data.with.LM$LMA_SE/data.with.LM$LMA)^2 )^0.5) * data.with.LM$LM
# data.with.LM$LM_SE = data.with.LM$LM_SE * c1 # unit conversion from gDM to gC
# 
# # Daily LM plot with SE
# plots = list() 
# plots[[1]] = ggplot(data.with.LMA, aes(x=Date, y=LM, group = interaction(T_treatment,chamber_type), colour=T_treatment, shape=chamber_type)) + 
#   # geom_point(position=pd) +
#   geom_ribbon(data = data.with.LMA, aes(ymin=LM-LM_SE, ymax=LM+LM_SE), linetype=2, alpha=0.1,size=0.1) +
#   # geom_errorbar(position=pd,aes(ymin=LMA-LMA_SE, ymax=LMA+LMA_SE), colour="grey", width=2) +
#   geom_line(data = data.with.LMA, aes(x = Date, y = LM, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) + 
#   ylab(expression("LM"~"(g C "*plant^"-1"*")")) +
#   scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.with.LMA$Date), max(data.with.LMA$Date)+1)) +
#   labs(colour="Temperature", shape="Chamber type", linetype="Chamber type") +
#   scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
#   theme_bw() +
#   theme(legend.title = element_text(colour="black", size=font.size)) +
#   theme(legend.text = element_text(colour="black", size = font.size)) +
#   theme(legend.position = c(0.25,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
#   theme(legend.key = element_blank()) +
#   theme(text = element_text(size=font.size)) +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
# 
# # png("output/3.LM.png", units="px", width=2200, height=1600, res=220)
# # plots
# # dev.off()
# 
# #----------------------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
####### From John's script #########
#- Get an estimate of branch, stem, and leaf mass as well as leaf area for each day of the experiment
treeMass.daily <- returnTreeMass()
# treeMassFlux <- merge(dat.hr.gf,treeMass,by=c("chamber","Date","T_treatment","Water_treatment"))
treeMass.daily = subset(treeMass.daily, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-05-26"))

treeMass.daily$Measurement[treeMass.daily$Date %in% as.Date("2014-05-26")] = 40
treeMass = na.omit(treeMass.daily) # consider only fortnightly direct measurements of H and D
# treeMass = treeMass.daily # consider daily data interpolated from fortnightly direct measurements of H and D

treeMass$woodMass = treeMass$branchMass + treeMass$boleMass
treeMass$chamber_type = as.factor( ifelse(treeMass$chamber %in% drought.chamb, "drought", "watered") )
treeMass.sum = summaryBy(leafMass+woodMass ~ Date+T_treatment+chamber_type, data=treeMass, FUN=c(mean,standard.error))
names(treeMass.sum)[4:7] = c("LM","WM","LM_SE","WM_SE")
treeMass.sum[,c(4:7)] = treeMass.sum[,c(4:7)] * c1 # unit conversion from gDM to gC

data.biomass = merge(rootmass, treeMass.sum, by = c("Date", "T_treatment", "chamber_type"), all=TRUE)

# Daily LM (using mean LMA, harvest LMA), WM and RM plot with SEs
# plots[[2]] = ggplot(data.biomass, aes(x=Date, y=LM_harvest, group = interaction(T_treatment,chamber_type), colour=T_treatment, shape=chamber_type)) + 
#   # geom_point(position=pd) +
#   geom_ribbon(data = data.biomass, aes(ymin=LM_harvest-LM_harvest_SE, ymax=LM_harvest+LM_harvest_SE), linetype=2, alpha=0.1,size=0.1) +
#   # geom_errorbar(position=pd,aes(ymin=LMA-LMA_SE, ymax=LMA+LMA_SE), colour="grey", width=2) +
#   geom_line(data = data.biomass, aes(x = Date, y = LM_harvest, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) + 
#   ylab(expression("Leaf Mass using harset LMA"~"(g C "*plant^"-1"*")")) +
#   scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
#   labs(colour="Temperature", shape="Chamber type", linetype="Chamber type") +
#   scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
#   theme_bw() +
#   theme(legend.title = element_text(colour="black", size=font.size)) +
#   theme(legend.text = element_text(colour="black", size = font.size)) +
#   theme(legend.position = c(0.3,0.88), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
#   theme(legend.key = element_blank()) +
#   theme(text = element_text(size=font.size)) +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plots = list()
pd <- position_dodge(5)
plots[[1]] = ggplot(data.biomass, aes(x=Date, y=WM, group = interaction(T_treatment,chamber_type), colour=T_treatment, shape=chamber_type)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=WM-WM_SE, ymax=WM+WM_SE), colour="grey", width=1) +
  # geom_ribbon(data = data.biomass, aes(ymin=WM-WM_SE, ymax=WM+WM_SE), linetype=2, alpha=0.1,size=0.1) +
  # geom_line(data = data.biomass, aes(x = Date, y = WM, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) + 
  ylab(expression("Wood Mass"~"(g C "*plant^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Temperature", shape="Chamber type", linetype="Chamber type") +
  scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.2,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pd <- position_dodge(5)
plots[[2]] = ggplot(data.biomass, aes(x=Date, y=RM, group = interaction(T_treatment,chamber_type), colour=T_treatment, shape=chamber_type)) + 
  geom_point(position=pd,size=2) +
  # geom_ribbon(position=pd,data = data.biomass, aes(ymin=RM-RM_SE, ymax=RM+RM_SE), linetype=2, alpha=0.1,size=0.1) +
  geom_errorbar(position=pd,aes(ymin=RM-RM_SE, ymax=RM+RM_SE), colour="grey", width=1) +
  # geom_line(data = data.biomass, aes(x = Date, y = RM, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) +
  ylab(expression("Root Mass"~"(g C "*plant^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Temperature", shape="Chamber type", linetype="Chamber type") +
  scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.2,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(file = "output/2.Biomass_litter_tnc.pdf")
print (do.call(grid.arrange,  plots))

plots = list()
plots[[1]] = ggplot(data.biomass, aes(x=Date, y=LM, group = interaction(T_treatment,chamber_type), colour=T_treatment, shape=chamber_type)) + 
  geom_point(position=pd) +
  # geom_ribbon(data = data.biomass, aes(ymin=LM-LM_SE, ymax=LM+LM_SE), linetype=2, alpha=0.1,size=0.1) +
  geom_errorbar(position=pd,aes(ymin=LM-LM_SE, ymax=LM+LM_SE), colour="grey", width=1) +
  # geom_line(data = data.biomass, aes(x = Date, y = LM, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) + 
  ylab(expression("Leaf Mass"~"(g C "*plant^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Temperature", shape="Chamber type", linetype="Chamber type") +
  scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.2,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
# Script to setup a model for the leaf mass
litterfall = read.csv("raw_data/WTC_TEMP_CM_LEAFLITTER_20130913-20140528_L1.csv")
litterfall$startDate = as.Date(litterfall$startDate)
litterfall$collectionDate = as.Date(litterfall$collectionDate)
litterfall$Date <- (litterfall$startDate + ((litterfall$collectionDate - litterfall$startDate) / 2))
litterfall = subset(litterfall, Date >= as.Date("2013-09-14") & Date <= as.Date("2014-05-27"))

# convert to data.table in place
litterfall = setDT(litterfall)
# dcast and do individual sums
litterfall.cast = dcast.data.table(litterfall, chamber ~ Date, value.var = 'litter', fun.aggregate = sum)
# cumsum
litterfall.cum <- litterfall.cast[, as.list(cumsum(unlist(.SD))), by = chamber]

litterfall.cum.melt <- melt(litterfall.cum, id.vars = "chamber")
litterfall.cum.melt = merge(litterfall.cum.melt, unique(treeMass[,c("chamber","T_treatment")]), all=TRUE)
litterfall.cum.melt$chamber_type = as.factor( ifelse(litterfall.cum.melt$chamber %in% drought.chamb, "drought", "watered") )
names(litterfall.cum.melt)[2:3] = c("Date","litter")
litterfall.cum.melt$Date = as.Date(litterfall.cum.melt$Date)
litterfall.cum.melt = summaryBy(litter ~ Date+T_treatment+chamber_type, data=litterfall.cum.melt, FUN=c(mean,standard.error))
names(litterfall.cum.melt)[4:5] = c("litter","litter_SE")


litterfall.initial = data.frame(Date = rep(as.Date("2013-09-17"), 4),
                                T_treatment = rep(unique(data$T_treatment), each=2),
                                chamber_type = rep(unique(data$chamber_type), 2),
                                litter = rep(0,4),
                                litter_SE = rep(0,4))
litterfall.cum.melt = rbind(litterfall.initial, litterfall.cum.melt)
litterfall.cum.melt$litter = litterfall.cum.melt$litter * c1 # unit conversion from gDM to gC
litterfall.cum.melt$litter_SE = litterfall.cum.melt$litter_SE * c1 # unit conversion from gDM to gC

data.biomass = merge(data.biomass, litterfall.cum.melt, all=TRUE)

pd <- position_dodge(5)
plots[[2]] = ggplot(data.biomass[complete.cases(data.biomass$litter),], aes(x=Date, y=litter, group = interaction(T_treatment,chamber_type), colour=T_treatment, shape=chamber_type)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=litter-litter_SE, ymax=litter+litter_SE), colour="grey", width=1) +
  # geom_line(position=pd, data = data.biomass[complete.cases(data.biomass$litter),], aes(x = Date, y = litter, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) + 
  ylab(expression("Leaf Litter"~"(g C "*plant^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Temperature", shape="Chamber type", linetype="Chamber type") +
  scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.25,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

print (do.call(grid.arrange,  plots))
# dev.off()

# data.with.litter$litter_ratio = data.with.litter$litter / data.with.litter$LM * 100
# data.with.litter$litter_ratio_SE = (((data.with.litter$litter_SE/(3)^0.5)^2 / (data.with.litter$LM_SE/(3)^0.5)^2 )/2 )^0.5

# plots[[3]] = ggplot(data.with.litter[complete.cases(data.with.litter$litter),], aes(x=Date, y=litter_ratio, group = interaction(T_treatment,chamber_type), colour=T_treatment, shape=chamber_type)) + 
#   geom_point() +
#   geom_errorbar(aes(ymin=litter_ratio - litter_ratio_SE, ymax=litter_ratio + litter_ratio_SE), colour="grey", width=2) +
#   geom_line(data = data.with.litter[complete.cases(data.with.litter$litter),], aes(x = Date, y = litter_ratio, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) + 
#   ylab("Litterfall ratio (% of LM)") +
#   scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(melted.data.sub$Date), max(melted.data.sub$Date)+1)) +
#   labs(colour="Temperature", shape="Chamber type", linetype="Chamber type") +
#   scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
#   theme_bw() +
#   theme(legend.title = element_text(colour="black", size=font.size)) +
#   theme(legend.text = element_text(colour="black", size = font.size)) +
#   theme(legend.position = c(0.25,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
#   theme(legend.key = element_blank()) +
#   theme(text = element_text(size=font.size)) +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# png("output/3.LM_litterfall.png", units="px", width=2200, height=1600, res=220)
# pdf(file = "output/3.LM_litterfall.pdf")
# print (do.call(grid.arrange,  plots))
# dev.off()

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Script to read the leaf storage data 
# read the leaf tnc data over time
tnc = read.csv("raw_data/WTC_TEMP_CM_LEAFCARB_20130515-20140402_L2.csv") # unit = mg of tnc per g of dry leafmass
tnc = na.omit(tnc)
tnc$Date = as.Date(tnc$Date)
tnc = subset(tnc, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-05-27"))
tnc$chamber_type = as.factor( ifelse(tnc$chamber %in% drought.chamb, "drought", "watered") )
keeps <- c("Date", "chamber", "T_treatment", "chamber_type", "TNC_mgg")
tnc = tnc[ , keeps, drop = FALSE]
tnc.mean <- summaryBy(TNC_mgg ~ Date+T_treatment+chamber_type, data=tnc, FUN=c(mean,standard.error))
names(tnc.mean)[4:5] = c("tnc.conc","tnc.conc_SE")

# read the diurnal leaf tnc data - 2-hourly intervals during a sunny (20 February) and an overcast/rainy day (26 March)
tnc.diurnal = read.csv("raw_data/WTC_TEMP_CM_LEAFCARB-DIURNAL_20140220-20140326_R.csv") # unit = mg of tnc per g of dry leafmass
tnc.diurnal$Sample.Date <- parse_date_time(tnc.diurnal$Sample.Date,"d m y")
tnc.diurnal$Sample.Date = as.Date(tnc.diurnal$Sample.Date, format = "%d/%m/%Y")
keeps <- c("Sample.Date", "Chamber", "TNC")
tnc.diurnal = tnc.diurnal[ , keeps, drop = FALSE]
names(tnc.diurnal) = c("Date","chamber","TNC")
tnc.diurnal = merge(tnc.diurnal, unique(height.dia[,c("chamber","T_treatment","chamber_type")]), by="chamber")
tnc.diurnal.mean <- summaryBy(TNC ~ Date+T_treatment+chamber_type, data=tnc.diurnal, FUN=c(mean,standard.error))
names(tnc.diurnal.mean)[4:5] = c("tnc.conc","tnc.conc_SE")

# read the diurnal leaf tnc data - trees were girdled
tnc.girdled = read.csv("raw_data/WTC_TEMP_CM_PETIOLEGIRDLE-LEAFMASS-AREA-CARB_20140507-20140512_L2.csv") # unit = mg of tnc per g of dry leafmass
tnc.girdled$Date <- parse_date_time(tnc.girdled$Date,"y m d")
tnc.girdled$Date = as.Date(tnc.girdled$Date, format = "%Y/%m/%d")
tnc.girdled = subset(tnc.girdled, treatment %in% as.factor("control"))
tnc.girdled$TNC = tnc.girdled$starch + tnc.girdled$soluble_sugars
tnc.girdled = merge(tnc.girdled, unique(height.dia[,c("chamber","T_treatment","chamber_type")]), by="chamber")
tnc.girdled.mean <- summaryBy(TNC ~ Date+T_treatment+chamber_type, data=tnc.girdled, FUN=c(mean,standard.error))
names(tnc.girdled.mean)[4:5] = c("tnc.conc","tnc.conc_SE")

# bind all available tnc data into one dataframe
tnc.mean = rbind(tnc.mean, tnc.diurnal.mean)
tnc.mean = rbind(tnc.mean, tnc.girdled.mean)

# unit conversion: 1 g tnc has 0.4 gC (12/30) and 1 g plant has 0.48 gC
tnc.mean$tnc.conc = tnc.mean$tnc.conc / 1000 # g of tnc per g of dry leafmass
tnc.mean$tnc.conc_SE = tnc.mean$tnc.conc_SE / 1000 # g of tnc per g of dry leafmass
tnc.mean$tnc.conc = tnc.mean$tnc.conc * 0.4 / c1 # gC in tnc per gC of dry leafmass
tnc.mean$tnc.conc_SE = tnc.mean$tnc.conc_SE * 0.4 / c1 # gC in tnc per gC of dry leafmass

# get the daily interpotalerd LM
treeMass.daily$chamber_type = as.factor( ifelse(treeMass.daily$chamber %in% drought.chamb, "drought", "watered") )
treeMass.daily.sum = summaryBy(leafMass ~ Date+T_treatment+chamber_type, data=treeMass.daily, FUN=c(mean,standard.error))
names(treeMass.daily.sum)[4:5] = c("LM","LM_SE")
treeMass.daily.sum[,c(4:5)] = treeMass.daily.sum[,c(4:5)] * c1 # unit conversion from gDM to gC

tnc.final = merge(tnc.mean, treeMass.daily.sum[,c("Date","T_treatment","chamber_type","LM","LM_SE")], by=c("Date","T_treatment","chamber_type"), all=FALSE)
tnc.final$TNC = tnc.final$tnc.conc * tnc.final$LM # Unit = gC in tnc per gC in plant
tnc.final$TNC_SE = tnc.final$tnc.conc_SE * tnc.final$LM # Unit = gC in tnc per gC in plant

data.biomass = merge(data.biomass, tnc.final[,c("Date", "T_treatment", "chamber_type", "TNC", "TNC_SE")], all=TRUE)

plots = list()
pd <- position_dodge(1)
plots[[1]] = ggplot(data.biomass[complete.cases(data.biomass$TNC),], aes(x=Date, y=TNC, group = interaction(T_treatment,chamber_type), colour=T_treatment, shape=chamber_type)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=TNC-TNC_SE, ymax=TNC+TNC_SE), colour="grey", width=1) +
  # geom_line(position=pd, data = data.biomass[complete.cases(data.biomass$TNC),], aes(x = Date, y = TNC, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) + 
  ylab(expression("TNC"~"(g C "*plant^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Temperature", shape="Chamber type", linetype="Chamber type") +
  scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.25,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

df <- data.frame()
plots[[2]] = ggplot(df)

print (do.call(grid.arrange,  plots))
dev.off()

# write the file with all available inputs and data (GPP, Ra, LA, rootmass, woodmass, foliagemass, litterfall, tnc)
write.csv(data.biomass, file = "processed_raw_data/data_biomass_litter_tnc.csv", row.names = FALSE)

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Calculate daily below ground root respiration rates
q25_drake = 2.25

# Coarse root respiration rates
branch.resp = read.csv("raw_data/WTC_TEMP_CM_GX-RBRANCH_20140513-20140522_L1_v1.csv") # Rbranch: branch respiration (nmol CO2 g-1 s-1) 
branch.resp$date = as.Date(branch.resp$date)
branch.resp = subset(branch.resp, date %in% as.Date("2014-05-13")) # Only consider the pre-girdling data
branch.resp$chamber_type = as.factor( ifelse(branch.resp$chamber %in% drought.chamb, "drought", "watered") )
branch.resp$Treatment <- as.factor(paste(branch.resp$T_treatment, branch.resp$chamber_type))

# Test for any significant difference between the treatment groups
boxplot(branch.resp$Rbranch ~ branch.resp$Treatment, xlab="Treatment", ylab=(expression("Branch wood respiration"~"(nmol CO2 "*g^"-1"*" "*s^"-1"*")")))

summary(aov(Rbranch ~ T_treatment * chamber_type, data = branch.resp)) # YES, there is significant difference only accross the temperature treatments
# summary(aov(Rbranch ~ Treatment, data = branch.resp)) # YES, there is significant difference accross the treatments
# t.test(branch.resp$Rbranch ~ branch.resp$T_treatment) # YES, there is significant difference accross temperatures
# t.test(branch.resp$Rbranch ~ branch.resp$chamber_type) # NO, there is no significant difference accross drought/watered treatments

############ So how to group the treatments????????
rd15.root <- summaryBy(Rbranch ~ T_treatment, data=branch.resp, FUN=c(mean))
names(rd15.root)[ncol(rd15.root)] = c("rd15.coarseroot")
rd15.root$rd15.coarseroot = rd15.root$rd15.coarseroot * (10^-9 * 12) * (3600 * 24) * (1/c1) # unit conversion from nmolCO2 g-1 s-1 to gC gC-1 d-1
# rd15.coarseroot$rd15.coarseroot_SE = rd15.coarseroot$rd15.coarseroot_SE * (10^-9 * 12) * (3600 * 24) * (1/c1) # unit conversion from nmolCO2 g-1 s-1 to gC gC-1 d-1

# Bole and big tap root respiration rates
bole.resp = read.csv("raw_data/WTC_TEMP_CM_WTCFLUX-STEM_20140528_L1_v1.csv") # Rbranch: branch respiration (nmol CO2 g-1 s-1) 
bole.resp$chamber_type = as.factor( ifelse(bole.resp$chamber %in% drought.chamb, "drought", "watered") )
bole.resp$Treatment <- as.factor(paste(bole.resp$T_treatment, bole.resp$chamber_type))

# Test for any significant difference between the treatment groups
boxplot(bole.resp$R_stem_nmol ~ bole.resp$Treatment, xlab="Treatment", ylab=(expression("Bole wood respiration"~"(nmol CO2 "*g^"-1"*" "*s^"-1"*")")))

summary(aov(R_stem_nmol ~ T_treatment * chamber_type, data = bole.resp)) # NO, there is no significant difference accross the treatments
# summary(aov(R_stem_nmol ~ Treatment, data = bole.resp)) # NO, there is no significant difference accross the treatments
# t.test(bole.resp$R_stem_nmol ~ bole.resp$T_treatment) # NO, there is no significant difference accross tepmeratures
# t.test(bole.resp$R_stem_nmol ~ bole.resp$chamber_type) # NO, there is no significant difference accross drought/watered treatments

rd15.root$rd15.boleroot = mean(bole.resp$R_stem_nmol)
rd15.root$rd15.boleroot = rd15.root$rd15.boleroot * (10^-9 * 12) * (3600 * 24) * (1/c1) # unit conversion from nmolCO2 g-1 s-1 to gC gC-1 d-1

# Fine root respiration rates (Constant)
# Fine root respiration rate = 10 nmolCO2 g-1 s-1 (Ref: Drake et al. 2017: GREAT exp data; Mark's Email)
rd25.fineroot = 10 * (10^-9 * 12) * (3600 * 24) * (1/c1) # unit conversion from nmolCO2 g-1 s-1 to gC gC-1 d-1
rd15.root$rd15.fineroot = rd25.fineroot * q25_drake^((15-25)/10)

# Intermediate root respiration rates
rd15.root$rd15.intermediateroot = exp ((log(rd15.root$rd15.coarseroot) + log(rd15.root$rd15.fineroot))/2 ) # unit = gC gC-1 d-1


#----------------------------------------------------------------------------------------------------------------
# import site weather data, take only soil temperatures at 10 cm depth, format date stuff
files <- list.files(path = "raw_data/WTC_TEMP_CM_WTCMET", pattern = ".csv", full.names = TRUE)
temp <- lapply(files, fread, sep=",")
met.data <- rbindlist( temp )

met.data <- met.data[ , c("chamber","DateTime","SoilTemp_Avg.1.","SoilTemp_Avg.2.")]
met.data$SoilTemp <- rowMeans(met.data[,c("SoilTemp_Avg.1.","SoilTemp_Avg.2.")], na.rm=TRUE)
met.data$Date <- as.Date(met.data$DateTime)

# need to turn the datetime into hms
met.data$DateTime <- ymd_hms(met.data$DateTime)
met.data$time <- format(met.data$DateTime, format='%H:%M:%S')

# met.data.na = met.data[is.na(met.data$SoilTemp),] # Check any NA values for soil temperature

# subset by Date range of experiment
met.data <- subset(met.data[, c("chamber","Date","time","SoilTemp")], Date  >= "2013-09-17" & Date  <= "2014-05-26")
met.data$chamber = as.factor(met.data$chamber)
met.data = merge(met.data, unique(height.dia[,c("chamber","T_treatment")]), by="chamber")

# need to calculate Rdark through time using rdarkq10 equation by treatment
met.data <- merge(met.data, rd15.root, by=c("T_treatment"))

met.data[,c("Rd.fineroot","Rd.intermediateroot","Rd.coarseroot","Rd.boleroot")] <- 
  with(met.data, met.data[,c("rd15.fineroot","rd15.intermediateroot","rd15.coarseroot","rd15.boleroot")] * 
         q25_drake^((SoilTemp-15)/10)) # unit (gC per gC root per day)

Rd <- summaryBy(Rd.fineroot+Rd.intermediateroot+Rd.coarseroot+Rd.boleroot ~ Date+T_treatment, data=met.data, FUN=mean, keep.names=TRUE , na.rm=TRUE) # Sum of all same day Rd
# colSums(is.na(Rd)) # Check any NA values for Rd

write.csv(Rd, "processed_raw_data/Rd.csv", row.names=FALSE) # unit: gC per gC plant per day

#----------------------------------------------------------------------------------------------------------------
# Plot various root respiration data
Rd.melt <- melt(Rd, id.vars = c("Date","T_treatment"))
i = 0
font.size = 15
plot = list() 
meas = as.factor(c("Rd.fineroot","Rd.intermediateroot","Rd.coarseroot","Rd.boleroot"))
title = as.character(c("A","B","C","D"))
pd <- position_dodge(0) # move the overlapped errorbars horizontally
for (p in 1:length(meas)) {
  Rd.melt.sub = subset(Rd.melt,variable %in% meas[p])
  
  i = i + 1
  plot[[i]] = ggplot(Rd.melt.sub, aes(x=Date, y=value, group = T_treatment, colour=T_treatment)) + 
    geom_point(position=pd) +
    geom_line(position=pd,data = Rd.melt.sub, aes(x = Date, y = value, group = T_treatment, colour=T_treatment)) + 
    ylab(expression(R[d(fineroot)]~"(g C "*g^"-1"*" C "*d^"-1"*")")) + xlab("") +
    scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(Rd$Date)-2, max(Rd$Date)+2)) +
    labs(colour="Temperature") +
    scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
    theme_bw() +
    theme(legend.title = element_text(colour="black", size=font.size)) +
    theme(legend.text = element_text(colour="black", size=font.size)) +
    theme(legend.position = c(0.85,0.85), legend.box = "horizontal") + theme(legend.key.height=unit(0.9,"line")) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  if (p==2) {
    plot[[i]] = plot[[i]] + ylab(expression(R[d(intermediateroot)]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
    # plot[[i]] = plot[[i]] + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 0.75), units="line"))
  } 
  if (p==3) {
    plot[[i]] = plot[[i]] + ylab(expression(R[d(coarseroot)]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
  } 
  if (p==4) {
    plot[[i]] = plot[[i]] + ylab(expression(R[d(boleroot)]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
  }
}
pdf("output/3.Rd.pdf", width=20, height=5)
plot[1]
plot[2]
plot[3]
plot[4]
dev.off()
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Calculate the ratios of harvest rootmass
rootmass.harvest$RootDMbole = rootmass.harvest$RootDMbole + rootmass.harvest$CoreRoots
rootmass.harvest$FRratio = rootmass.harvest$RootDM0to2 / rootmass.harvest$RootDMtotal
rootmass.harvest$IRratio = rootmass.harvest$RootDM2to10 / rootmass.harvest$RootDMtotal
rootmass.harvest$CRratio = rootmass.harvest$RootDM10up / rootmass.harvest$RootDMtotal
rootmass.harvest$BRratio = rootmass.harvest$RootDMbole / rootmass.harvest$RootDMtotal
rootmass.harvest$Treatment <- as.factor(paste(rootmass.harvest$T_treatment, rootmass.harvest$chamber_type))

# Test for any significant difference between the treatment groups
summary(aov(FRratio ~ Treatment, data = rootmass.harvest)) # NO, there is no significant difference accross the treatments
summary(aov(IRratio ~ Treatment, data = rootmass.harvest)) # NO, there is no significant difference accross the treatments
summary(aov(CRratio ~ Treatment, data = rootmass.harvest)) # NO, there is no significant difference accross the treatments
summary(aov(BRratio ~ Treatment, data = rootmass.harvest)) # NO, there is no significant difference accross the treatments

rootmass.harvest.ratio <- summaryBy(FRratio+IRratio+CRratio+BRratio ~ Date, data=rootmass.harvest, FUN=c(mean,standard.error), keep.names=TRUE , na.rm=TRUE) # Sum of all same day Rd
names(rootmass.harvest.ratio)[2:ncol(rootmass.harvest.ratio)] = c("FRratio","IRratio","CRratio","BRratio","FRratio_SE","IRratio_SE","CRratio_SE","BRratio_SE")
rootmass.harvest.ratio.melt = melt(rootmass.harvest.ratio[,c("Date","FRratio","IRratio","CRratio","BRratio")], id.vars = "Date")
rootmass.harvest.ratio.melt$biomass = as.factor("Belowground")
rootmass.harvest.ratio.melt$SE = melt(rootmass.harvest.ratio[,c("Date","FRratio_SE","IRratio_SE","CRratio_SE","BRratio_SE")], id.vars = "Date")[,3]

#----------------------------------------------------------------------------------------------------------------
# Calculate the ratios of aboveground biomass
treeMass$LMratio = treeMass$leafMass / treeMass$totMass
treeMass$BMratio = treeMass$branchMass / treeMass$totMass
treeMass$SMratio = treeMass$boleMass / treeMass$totMass
treeMass$Treatment <- as.factor(paste(treeMass$T_treatment, treeMass$chamber_type))

# Test for any significant difference between the treatment groups
summary(aov(LMratio ~ Treatment, data = treeMass)) # NO, there is no significant difference accross the treatments
summary(aov(BMratio ~ Treatment, data = treeMass)) # NO, there is no significant difference accross the treatments
summary(aov(SMratio ~ Treatment, data = treeMass)) # NO, there is no significant difference accross the treatments


treeMass.ratio = summaryBy(LMratio+BMratio+SMratio ~ Date, data=treeMass, FUN=c(mean,standard.error))
names(treeMass.ratio)[2:ncol(treeMass.ratio)] = c("LMratio","BMratio","SMratio","LMratio_SE","BMratio_SE","SMratio_SE")
#----------------------------------------------------------------------------------------------------------------
# Plot various biomass ratios including different root types partitioning
treeMass.ratio.melt = melt(treeMass.ratio[,c("Date","LMratio","BMratio","SMratio")], id.vars = "Date")
treeMass.ratio.melt$biomass = as.factor("Aboveground")
treeMass.ratio.melt$SE = melt(treeMass.ratio[,c("Date","LMratio_SE","BMratio_SE","SMratio_SE")], id.vars = "Date")[,3]

treeMass.ratio.melt = rbind(treeMass.ratio.melt, rootmass.harvest.ratio.melt)

font.size = 10
plot = list() 
pd <- position_dodge(2) # move the overlapped errorbars horizontally
plot = ggplot(treeMass.ratio.melt, aes(x=Date, y=value, group = interaction(variable,biomass), colour=variable, shape=biomass)) + 
  geom_point(position=pd, size=2) +
  geom_errorbar(position=pd, aes(ymin=value-SE, ymax=value+SE), colour="grey", width=1) +
  geom_line(position=pd,data = treeMass.ratio.melt, aes(x = Date, y = value, group = variable, colour=variable)) + 
  ylab("Biomass ratio") + xlab("") +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(Rd$Date)-2, max(Rd$Date)+2)) +
  labs(colour="Biomass types", shape="") +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size=font.size)) +
  theme(legend.position = c(0.2,0.87), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

pdf("output/4.Biomass_ratio.pdf")
plot
dev.off()

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#- process leaf-scale R vs. T data to find q10 value for WTC3
rvt <- read.csv("raw_data/WTC_TEMP_CM_GX-RdarkVsT_20140207-20140423_L1.csv")
rvt$Date <- as.Date(rvt$Date)
rvt.sub <- subset(rvt,Date==as.Date("2014-02-07")) #- just pull out the data measured prior to the drought

rvt.45 <- subset(rvt.sub, Tleaf > 18 & Tleaf<=40)
rvt.45$lnRmass <- log(rvt.45$Rmass)
rvt.45$varT <- (rvt.45$Tleaf-25)/10

rvt.45$Tleaf_bin <- cut(rvt.45$Tleaf,breaks=seq(from=18,to=40,length=25))
rvt.45$Tleaf_bin_mid <- sapply(strsplit(gsub("^\\W|\\W$", "", rvt.45$Tleaf_bin), ","), function(x)sum(as.numeric(x))/2) 
rvt.treat.bin <- summaryBy(Rmass+lnRmass+varT~date+T_treatment+Water_treatment+Tleaf_bin_mid,data=rvt.45,keep.names=T,FUN=mean)

#- fit arrhenious. 
lmRleaf <- lm(lnRmass~varT+T_treatment+Water_treatment,data=rvt.treat.bin)
summary(lmRleaf)
coef(lmRleaf) 
# q10_amb_wet <- unname(exp(coef(lmRleaf)[2]))
# q10_ele_wet <- unname(exp(coef(lmRleaf)[2]+coef(lmRleaf)[3]))
# q10_amb_dry <- unname(exp(coef(lmRleaf)[2]+coef(lmRleaf)[4]))
# q10_ele_dry <- unname(exp(coef(lmRleaf)[2]+coef(lmRleaf)[3]+coef(lmRleaf)[4]))

# No significant difference for drought/control treatments, so only consider temperature treatment effect
lmRleaf_amb <- lm(lnRmass~varT,data=subset(rvt.treat.bin,T_treatment %in% as.factor("ambient")))
q10_amb <- unname(exp(coef(lmRleaf_amb)[2]))
rd25_amb = unname(exp(coef(lmRleaf_amb)[1]))
lmRleaf_ele <- lm(lnRmass~varT,data=subset(rvt.treat.bin,T_treatment %in% as.factor("elevated")))
q10_ele <- unname(exp(coef(lmRleaf_ele)[2]))
rd25_ele = unname(exp(coef(lmRleaf_ele)[1]))

# get model predictions of R across all temperatures
xvals_Rleaf <- seq(18,40, length=101)
predRleaf_amb <- exp(log(rd25_amb) + (log(q10_amb) * (xvals_Rleaf-25)/10))
predRleaf_ele <- exp(log(rd25_ele) + (log(q10_ele) * (xvals_Rleaf-25)/10))


#----------------------------------------------------------------------------------------------------------------
#- plot R vs. T
rvt.treat <- summaryBy(Rmass~Date+T_treatment+Tleaf_bin_mid,data=rvt.treat.bin,keep.names=F,FUN=c(mean,standard.error))
par(mfrow=c(1,1),cex.lab=1.5,mar=c(2,7,1,2),oma=c(4,1,0,0),las=1)
xlims=c(15,42)

# #- plot MASS BASED leaf R over T for the first date of high resolution T-response curves
plotBy(Rmass.mean~Tleaf_bin_mid|T_treatment,data=rvt.treat,xaxt="n",yaxt="n",
       ylab=expression(R[leaf]~(mu*mol~CO[2]~m^-2~s^-1)),col=c("black","red"),pch=1,
       xlim=xlims,type="p",lwd=3,cex=1,cex.lab=1.6,legend=F,
       panel.first=adderrorbars(x=rvt.treat$Tleaf_bin_mid,y=rvt.treat$Rmass.mean,SE=rvt.treat$Rmass.standard.error,direction="updown"))
magaxis(side=c(1,2,3,4),labels=c(1,1,0,1))
title(xlab=expression(Temperature~(degree*C)),xpd=NA)
lines(x=xvals_Rleaf,y=predRleaf_amb,col="black",lwd=2)
lines(x=xvals_Rleaf,y=predRleaf_ele,col="red",lwd=2)
legend("topleft",legend=unique(rvt.treat$T_treatment),col=c('black','red'),lty=1,bty="n",cex=0.8,pt.cex=1.2)

dev.copy2pdf(file="output/5.Rdark_vs_T.pdf")

# #- plot raw T-response curves
# plotBy(Rmass~Tleaf|chamber,data=subset(rvt.45,T_treatment %in% as.factor("ambient")),type="p",xlim=xlims,ylim=c(0,30),pch=20,size=0.2,cex=1.6,cex.lab=1.6,xaxt="n",yaxt="n",
#        ylab=expression(atop(R[canopy],(mu*mol~CO[2]~s^-1))),xlab="",legend=F)
#        # panel.first=adderrorbars(x=fits.trt$Tair.mean,y=fits.trt$Rcanopy_umol.mean,SE=fits.trt$Rcanopy_umol.standard.error,direction="updown"))
# lines(x=xvals_Rleaf,y=predRleafA,col="black",lwd=2)
# lines(x=xvals_Rleaf,y=predRleafE,col="red",lwd=2)
# legend("topleft",legend=unique(rvt.treat$T_treatment),col=c('black','red'),lty=1,bty="n",cex=0.8,pt.cex=1.2)
# magaxis(side=c(1,2,3,4),labels=c(1,1,0,1))


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# # Script to read and setup a model for the stem mass based on geometry usnig WTC3 harvest data 
# 
# height.dia.harvest = subset(height.dia, Date %in% max(height.dia$Date))
# harvest.wtc.sub = merge(harvest.wtc[,c("Date","chamber","T_treatment","chamber_type","BranchDM","StemDM")], height.dia.harvest[,c("chamber","height","diameter")], by="chamber")
# harvest.wtc.sub$woodmass = harvest.wtc.sub$BranchDM + harvest.wtc.sub$StemDM
#   
# wood.density = read.csv("raw_data/WTC_TEMP_CM_WOODDENSITY_20140528_L1.csv") # unit = g cm-3
# wood.density.mean <- summaryBy(wooddensity ~ chamber, data=wood.density, FUN=c(mean))
# names(wood.density.mean)[2] = c("wooddensity")
# 
# harvest.wtc.sub = merge(harvest.wtc.sub, wood.density.mean, by="chamber")
# harvest.wtc.sub$diameter = harvest.wtc.sub$diameter / 10 # unit consert to cm from mm
# 
# harvest.wtc.sub$stemvolume = pi / 4 * harvest.wtc.sub$height * (harvest.wtc.sub$diameter)^2
# harvest.wtc.sub$StemDM_calc = harvest.wtc.sub$stemvolume * harvest.wtc.sub$wooddensity
# harvest.wtc.sub$TF_total = harvest.wtc.sub$woodmass / harvest.wtc.sub$StemDM_calc
# 
# branch.data = read.csv("raw_data/WTC_TEMP_CM_BRANCHCENSUS_20130910-20140516_L0_v1.csv") # unit: branchlength = cm, branchdiameter = mm
# branch.data$Date = as.POSIXct(branch.data$Date,format="%d/%m/%Y",tz="GMT")
# branch.data$Date = as.Date(branch.data$Date)
# branch.data.harvest = subset(branch.data, Date %in% max(branch.data$Date))
# keeps <- c("Date", "chamber", "branchdiameter", "branchlength")
# branch.data.harvest = branch.data.harvest[ , keeps, drop = FALSE]
# branch.data.harvest$branchdiameter = branch.data.harvest$branchdiameter / 10 # unit consert to cm from mm
# branch.data.harvest = na.omit(branch.data.harvest)
# 
# branch.data.harvest$branchvolume = pi / 4 * branch.data.harvest$branchlength * (branch.data.harvest$branchdiameter)^2
# branch.data.harvest = merge(branch.data.harvest, harvest.wtc.sub[,c("chamber","wooddensity")], by="chamber")
# branch.data.harvest$BranchDM = branch.data.harvest$branchvolume * branch.data.harvest$wooddensity
# branch.data.harvest.mean <- summaryBy(BranchDM ~ chamber, data=branch.data.harvest, FUN=c(sum))
# names(branch.data.harvest.mean)[2] = "BranchDM_calc"
# 
# harvest.wtc.sub = merge(harvest.wtc.sub, branch.data.harvest.mean, by="chamber")
# harvest.wtc.sub$TF_branch = harvest.wtc.sub$BranchDM / harvest.wtc.sub$BranchDM_calc
# harvest.wtc.sub$TF_stem = harvest.wtc.sub$StemDM / harvest.wtc.sub$StemDM_calc
# 
# TF.mean <- summaryBy(TF_total+TF_stem+TF_branch ~ Date+T_treatment+chamber_type, data=harvest.wtc.sub, FUN=c(mean,standard.error))
# # names(TF.mean)[2] = c("wooddensity")



#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------





