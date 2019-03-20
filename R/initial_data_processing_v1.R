#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-- Script to analysis the raw WTC3 experiment data. 
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
###### R script to import and process the raw WTC3 experiment data 
# to model the carbon pools and fluxes using DA

# Clear the workspace (if needed)
rm(list=ls())

#- Load required libraries. There are quite a few (>20), including some non-standard functions that
#    are not on CRAN. This script will check for required libraries and install any that are missing.
source("R/loadLibraries_WTC3.R")

#- load the custom analysis and plotting functions that do all of the actual work
source("R/functions_WTC3.R")

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#- Script to read and setup a model for GPP, Ra (aboveground respiration) and leaf area
# read the hourly flux data
data.hr <- read.csv("data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
data.hr$DateTime <- as.POSIXct(data.hr$DateTime,format="%Y-%m-%d %H:%M:%S",tz="GMT")
data.hr$Date <- as.Date(data.hr$DateTime)

#- partition the net fluxes into GPP and Ra components. Specify the actitation energy and the number
#-  of prior nights to include for the estimate of the basal respiration rate (lagdates)
data.hr.p <- partitionHourlyFluxCUE_arr(data.hr.gf=data.hr,Ea=57.69,lagdates=3)

#- get daily sums from the partitioned hourly data
cue.list <- returnCUE.day(dat=data.hr.p) # get daily sums from hourly data
cue.day <- cue.list[[1]]                # extract chamber values on each day
cue.day.trt <- cue.list[[2]]            # extract treatment averages

#- plot GPP, Ra, and LA data over time
data=cue.day.trt
# data=subset(cue.day.trt,Water_treatment=="control")
names(data)[4:9] = c("GPP","Ra","LA","GPP_SE","Ra_SE","LA_SE")
# ignore the first few days of data to start on 2013-09-17 from where we have H & D measurements
data = subset(data, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-05-27"))

i = 0
font.size = 12
plot = list()
meas = as.factor(c("GPP","Ra","LA"))
error = as.factor(c("GPP_SE","Ra_SE","LA_SE"))
title = as.character(c("A","B","C"))
pd <- position_dodge(0.5)

data.sub <- summaryBy(GPP+Ra+LA ~ Date+T_treatment, data=data, FUN=c(mean,standard.error))
names(data.sub)[3:8] = c("GPP","Ra","LA","GPP_SE","Ra_SE","LA_SE")
melted.data = melt(data.sub[,c("GPP","Ra","LA","Date","T_treatment")], id.vars=c("Date","T_treatment"))
melted.data.error = melt(data.sub[,c("GPP_SE","Ra_SE","LA_SE","Date","T_treatment")], id.vars=c("Date","T_treatment"))
melted.data$SE = melted.data.error$value

for (p in 1:length(meas)) {
  melted.data.sub = subset(melted.data,variable %in% meas[p])
  i = i + 1
  plot[[i]] = ggplot(melted.data.sub, aes(x=Date, y=value, group = T_treatment, colour=T_treatment)) + 
    geom_ribbon(data = melted.data.sub, aes(ymin=value-SE, ymax=value+SE), linetype=2, alpha=0.1,size=0.1) +
    geom_line(position=pd,data = melted.data.sub, aes(x = Date, y = value, group = T_treatment, colour=T_treatment)) + 
    scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
    labs(colour="Temperature") +
    scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
    theme_bw() +
    theme(legend.title = element_text(colour="black", size=font.size)) +
    theme(legend.text = element_text(colour="black", size = font.size)) +
    theme(legend.position = c(0.15,0.75)) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  if (p==1) {
    plot[[i]] = plot[[i]] + ylab(expression(GPP~"(g C "*d^"-1"*")"))
  } else if (p==2) {
    plot[[i]] = plot[[i]] + ylab(expression(R[a]~"(g C "*d^"-1"*")"))
  } else if (p==3) {
    plot[[i]] = plot[[i]] + ylab(expression("Leaf area"~"("*m^"2"*" "*plant^"-1"*")"))
  }
  if (p>1) {
    plot[[i]] = plot[[i]] + guides(colour=FALSE)
  }
}

# png("output/1_1.GPP_Ra_LA_over_time.png", units="px", width=2200, height=1600, res=220)
pdf(file = "output/1_1.GPP_Ra_LA_over_time.pdf")
print (do.call(grid.arrange,  plot))

data.sub <- summaryBy(GPP+Ra+LA ~ Date+chamber_type, data=data, FUN=c(mean,standard.error))
names(data.sub)[3:8] = c("GPP","Ra","LA","GPP_SE","Ra_SE","LA_SE")
melted.data = melt(data.sub[,c("GPP","Ra","LA","Date","chamber_type")], id.vars=c("Date","chamber_type"))
melted.data.error = melt(data.sub[,c("GPP_SE","Ra_SE","LA_SE","Date","chamber_type")], id.vars=c("Date","chamber_type"))
melted.data$SE = melted.data.error$value

i = 0
plot = list()
for (p in 1:length(meas)) {
  melted.data.sub = subset(melted.data,variable %in% meas[p])
  i = i + 1
  plot[[i]] = ggplot(melted.data.sub, aes(x=Date, y=value, group = chamber_type, colour=chamber_type)) + 
    geom_ribbon(data = melted.data.sub, aes(ymin=value-SE, ymax=value+SE), linetype=2, alpha=0.1,size=0.1) +
    geom_line(position=pd,data = melted.data.sub, aes(x = Date, y = value, group = chamber_type, colour=chamber_type)) + 
    scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
    labs(colour="Treatment") +
    scale_color_manual(labels = c("watered_treat", "drought_treat"), values = c("blue", "red")) +
    theme_bw() +
    theme(legend.title = element_text(colour="black", size=font.size)) +
    theme(legend.text = element_text(colour="black", size = font.size)) +
    theme(legend.position = c(0.15,0.75)) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  if (p==1) {
    plot[[i]] = plot[[i]] + ylab(expression(GPP~"(g C "*d^"-1"*")"))
  } else if (p==2) {
    plot[[i]] = plot[[i]] + ylab(expression(R[a]~"(g C "*d^"-1"*")"))
  } else if (p==3) {
    plot[[i]] = plot[[i]] + ylab(expression("Leaf area"~"("*m^"2"*" "*plant^"-1"*")"))
  }
  if (p>1) {
    plot[[i]] = plot[[i]] + guides(colour=FALSE)
  }
}

print (do.call(grid.arrange,  plot))
dev.off()

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Script to read and plot stem height and diameter for various treatment cases 
# Read data form WTC3 experiment
height.dia.raw <- read.csv("data/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV")
height.dia.raw$DateTime = as.Date(height.dia.raw$DateTime)
flux <- read.csv("data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
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
  # keeps <- c("chamber", "DateTime", "T_treatment", "Water_treatment", "X15", "Plant_height")
  keeps <- c("chamber", "DateTime", "T_treatment", "Water_treatment", "X65", "Plant_height")
  height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
  names(height.dia.sub) <- c("chamber","Date","T_treatment","W_treatment","diameter","height")
  height.dia.sub$T_treatment = as.character(height.dia.sub$T_treatment)
  height.dia.sub$W_treatment = as.character(height.dia.sub$W_treatment)
  
  height.dia[(1+(i-1)*length(unique(height.dia$Date))) : (i*length(unique(height.dia$Date))), 
             c("T_treatment","W_treatment","diameter","height")] = height.dia.sub[,c("T_treatment","W_treatment","diameter","height")]
}
height.dia$T_treatment = as.factor(height.dia$T_treatment)
height.dia$W_treatment = as.factor(height.dia$W_treatment)

# Average the ambient and elevated temperature treatments regardless of drought/watered treatment
# n=6 for whole period
height.dia.idn <- summaryBy(height+diameter ~ Date+T_treatment, data=height.dia, FUN=c(mean,standard.error,length))
height.dia.idn$W_treatment = as.factor("all")
height.dia.idn$case = as.factor("1")

# Average the ambient and elevated temperature treatments considering the drought/watered treatment from February 2014
# 1st haf: n=6, 2nd half: n=3
height.dia.final <- summaryBy(height+diameter ~ Date+T_treatment+W_treatment, data=height.dia, FUN=c(mean,standard.error,length))
height.dia.final$case = as.factor("2")
height.dia.final = rbind(height.dia.idn, height.dia.final)

# Average the ambient and elevated temperature treatments considering the drought/watered treatment seperated from the start of the experiment
# n=3 for whole period
drought.chamb = unique(height.dia$chamber[ height.dia$W_treatment %in% as.factor("drydown")])
height.dia.sub = height.dia
height.dia.sub$chamber_type = as.factor( ifelse(height.dia.sub$chamber %in% drought.chamb, "drought_treat", "watered_treat") )
height.dia.idn <- summaryBy(height+diameter ~ Date+T_treatment+chamber_type, data=height.dia.sub, FUN=c(mean,standard.error,length))
names(height.dia.idn)[3] = "W_treatment"
height.dia.idn$case = as.factor("3")
height.dia.final = rbind(height.dia.final, height.dia.idn)
names(height.dia.final)[3:6] = c("height", "diameter", "height_SE", "diameter_SE")


# Plot both height and diameter for various scenarios
melted.height.dia.final = melt(height.dia.final[,c("diameter","height","Date","T_treatment","W_treatment","case")], id.vars=c("Date","T_treatment","W_treatment","case"))
melted.height.dia.error = melt(height.dia.final[,c("diameter_SE","height_SE","Date","T_treatment","W_treatment","case")], id.vars=c("Date","T_treatment","W_treatment","case"))
# melted.height.dia.final = melt(height.dia.final[,c("diameter","height","Date","T_treatment","case")], id.vars=c("Date","T_treatment","case"))
# melted.height.dia.error = melt(height.dia.final[,c("diameter_SE","height_SE","Date","T_treatment","case")], id.vars=c("Date","T_treatment","case"))
melted.height.dia.final$SE = melted.height.dia.error$value

# add.points = subset(melted.height.dia.final,Date %in% as.Date("2014-02-04") & W_treatment %in% as.factor("control"))
add.points = subset(melted.height.dia.final, Date <= as.Date("2014-02-04") & W_treatment %in% as.factor("control"))
add.points$W_treatment = as.factor("drydown")
melted.height.dia.final = rbind(melted.height.dia.final, add.points)

font.size = 12
plot = list() 
pd <- position_dodge(0) # move the overlapped errorbars horizontally
cbPalette = c("black", "green3", "red", "magenta", "blue")
melted.height.dia.final.trait = subset(melted.height.dia.final,variable %in% as.factor("diameter") & T_treatment %in% as.factor("elevated"))

plot = ggplot(melted.height.dia.final.trait, aes(x=Date, y=value, group = W_treatment, colour=W_treatment)) + 
  geom_point(position=pd) +
  geom_line(position=pd,data = melted.height.dia.final.trait, aes(x = Date, y = value, group = W_treatment, colour=W_treatment)) +
  ylab(paste(as.character(meas[p]),"(mm)")) +
  xlab("Month") +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
  labs(colour="Treatments") +
  scale_colour_manual(breaks=c("all","control","watered_treat","drydown","drought_treat"), labels=c("All Warm chambers (n=6)","Wet (n=6:predrought to n=3:postdrought)","Wet (n=3)",
                                                                                                    "Dry (n=6:predrought to n=3:postdrought)","Dry (n=3)"),values=cbPalette) +
  ggtitle(" Elevated temperature treatments") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size-2)) +
  theme(legend.position = c(0.35,0.88)) + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

dev.copy2pdf(file="output/6.Diameter_elevated_predrought_difference.pdf")

i = 0
font.size = 12
plots = list() 
meas = as.factor(c("diameter","height"))
error = as.factor(c("diameter_SE","height_SE"))
temp = as.factor(c("ambient","elevated"))
pd <- position_dodge(0) # move the overlapped errorbars horizontally
for (p in 1:length(meas)) {
  for (q in 1:length(temp)) {
    melted.height.dia.final.trait = subset(melted.height.dia.final,variable %in% meas[p] & T_treatment %in% temp[q])
    # melted.height.dia.final.trait = subset(melted.height.dia.final,variable %in% meas[p] & T_treatment %in% as.factor("elevated") & case %in% as.factor("1"))
    # melted.height.dia.final.trait = subset(melted.height.dia.final,variable %in% meas[p] & T_treatment %in% as.factor("ambient") & case %in% as.factor("3"))
    
    i = 2*(p-1) + q
    plots[[i]] = ggplot(melted.height.dia.final.trait, aes(x=Date, y=value, group = W_treatment, colour=W_treatment)) + 
      geom_point(position=pd) +
      geom_errorbar(position=pd,aes(ymin=value-SE, ymax=value+SE), colour="grey", width=2) +
      geom_line(position=pd,data = melted.height.dia.final.trait, aes(x = Date, y = value, group = W_treatment, colour=W_treatment)) +
      ylab(paste(as.character(meas[p]),"(mm)")) +
      xlab("Month") +
      scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
      labs(colour="Treatments") +
      # scale_colour_manual(breaks=c("all","control","watered_treat","drydown","drought_treat"), labels=c("All Warm chambers (n=6)","Wet (n=6:predrought to n=3:postdrought)","Wet (n=3)",
      #                                                "Dry (n=6:predrought to n=3:postdrought)","Dry (n=3)"),values=cbPalette) +
      ggtitle(paste(as.character(temp[q]), "temperature")) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.title = element_text(colour="black", size=font.size)) +
      theme(legend.text = element_text(colour="black", size = font.size-2)) +
      theme(legend.position = c(0.35,0.88)) + theme(legend.key.height=unit(0.8,"line")) +
      theme(legend.key = element_blank()) +
      theme(text = element_text(size=font.size)) +
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    if (p==2) {
      plots[[i]] = plots[[i]] + ylab(paste(as.character(meas[p]),"(cm)"))
    }
    if (q==1) {
      plots[[i]] = plots[[i]] + guides(colour=FALSE)
    }
  }
}

for (r in 1:length(meas)) {
  melted.height.dia.final.trait = subset(melted.height.dia.final,variable %in% meas[r])
  
  i = i + 1
  plots[[i]] = ggplot(melted.height.dia.final.trait, aes(x=Date, y=value, group = interaction(T_treatment,W_treatment), colour=T_treatment, shape=W_treatment)) + 
    geom_point(position=pd) +
    # geom_errorbar(position=pd,aes(ymin=value-SE, ymax=value+SE), colour="grey", width=2) +
    geom_line(position=pd,data = melted.height.dia.final.trait, aes(x = Date, y = value, group = interaction(T_treatment,W_treatment), colour=T_treatment, linetype=W_treatment)) +
    ylab(paste(as.character(meas[r]),"(mm)")) +
    xlab("Month") +
    scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
    labs(colour="Temperature", shape="Water", linetype="Water") +
    scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
    theme_bw() +
    theme(legend.title = element_text(colour="black", size=font.size)) +
    theme(legend.text = element_text(colour="black", size = font.size)) +
    theme(legend.position = c(0.9,0.3)) + theme(legend.key.height=unit(0.7,"line")) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  if (r==2) {
    plots[[i]] = plots[[i]] + ylab(paste(as.character(meas[r]),"(cm)"))
  }
}

png("output/2.dia_over_time.png", units="px", width=2200, height=1600, res=220)
lay <- rbind(c(1,1),c(2,3))
grid.arrange(grobs = plots[c(5,1,2)], layout_matrix = lay)
# lay <- rbind(c(1,1),c(1,1),c(2,2),c(3,3))
# grid.arrange(grobs = plots[c(5,1,2)], layout_matrix = lay)
dev.off()
png("output/2.height_over_time.png", units="px", width=2200, height=1600, res=220)
grid.arrange(grobs = plots[c(6,3,4)], layout_matrix = lay)
dev.off()

################# Based on the height and dia data, we decided to go with Case 3 that represents drought and watered treatments seperately from the beginning of the experiment 
################# so total number of treatment scenario = 4 (ambient+watered, ambient+droughted, elevated+watered, elevated+droughted) with sample size of n = 3 for each option

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Script to read and setup a model for the root mass from sink limited pot experiment and WTC3 harvest dat to predict 
# initial root mass from Height and Diameter
plot.summary.pot = read.csv("raw_data/plot_summary.csv")
initial.data.pot <- read.csv("raw_data/seedling_initial.csv")
keeps <- c("root_mass", "diameter_15", "height")
data.pot = initial.data.pot[ , keeps, drop = FALSE]
names(data.pot) <- c("rootmass","diameter","height")
data.pot$volume = as.integer(1000)

harvest.data.pot <- read.csv("raw_data/seedling_mass.csv")
harvest.data.pot$rootmass = harvest.data.pot$coarseroot + harvest.data.pot$fineroot
keeps <- c("rootmass", "plot", "pot", "volume")
harvest.data.pot = harvest.data.pot[ , keeps, drop = FALSE]

height.dia.pot <- read.csv("raw_data/height_diameter.csv")
height.dia.pot$Date <- parse_date_time(height.dia.pot$Date,"d m y")
height.dia.pot$Date = as.Date(height.dia.pot$Date, format = "%d/%m/%Y")
height.dia.pot = subset(height.dia.pot,Date=="2013-05-21")

height.dia.pot <- merge(plot.summary.pot, height.dia.pot, by=c("pot","plot"))
harvest.data.pot <- merge(height.dia.pot, harvest.data.pot, by=c("pot","plot","volume"))
keeps <- c("rootmass", "diameter", "height","volume")
harvest.data.pot = harvest.data.pot[ , keeps, drop = FALSE]
data.pot = rbind(data.pot, harvest.data.pot)
# data.pot[,"rootmass"] = data.pot[,"rootmass"] * 0.48 # unit conversion: gDM to gC 
data.free.pot = subset(data.pot, volume %in% 1000)
data.free.pot = data.free.pot[,-4]

# Read data form WTC3 experiment
# flux <- read.csv("data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
# height.dia <- read.csv("data/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV")
# height.dia$DateTime = as.Date(height.dia$DateTime)
rootmass.harvest = read.csv("data/WTC_TEMP_CM_HARVEST-ROOTS_20140529-20140606_L1_v1.csv")
rootmass.harvest$chamber_type = as.factor( ifelse(rootmass.harvest$chamber %in% drought.chamb, "drought_treat", "watered_treat") )
# rootmass.harvest = subset(rootmass.harvest, chamber_type %in% "watered_treat")


# height.dia.initial <- summaryBy(height+diameter ~ Date+T_treatment+W_treatment, data=height.dia.initial, FUN=c(mean,standard.error))
# names(height.dia.initial)[4:7] = c("height","diameter","height_SE","diameter_SE")
# keeps <- c("chamber", "Plant_height", "X15")
# height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]

# chambers = unique(flux$chamber)
height.dia.sub = subset(height.dia, Date %in% as.Date("2014-05-27") & T_treatment %in% as.factor("ambient"))
keeps <- c("chamber", "height", "diameter")
height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
height.dia.sub = height.dia.sub[complete.cases(height.dia.sub), ]
height.dia.sub <- merge(height.dia.sub, rootmass.harvest[,c("chamber","RootDMtotal","chamber_type")], by=c("chamber"))
height.dia.sub = subset(height.dia.sub, chamber_type %in% as.factor("watered_treat"))
keeps <- c("RootDMtotal", "diameter", "height")
height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
names(height.dia.sub) <- c("rootmass","diameter","height")

data.free.pot = rbind(data.free.pot, height.dia.sub) # Merge pot exp data (initial and harvest of free seedlings) and WTC3 harvest data (ambient)

# Fit a linear regression to estimate root mass by stem dia and height (ignoring temperature variation)
# the regression is fitted using the Sink limited pot experiemnt of Court Campany for similar seedlings (Euc.Teri.)
rm1 <- lm(log(rootmass) ~ log(diameter) + log(height), data=data.free.pot)
rm2 <- lm(log(rootmass) ~ log(diameter) * log(height), data=data.free.pot)

png("output/2.model_comparison.png", units="px", width=2200, height=1600, res=220)
layout(matrix(c(1,2,3,4),2,2,byrow=TRUE), widths=c(1,1), heights=c(10,10))
par(oma = c(0, 0, 0, 0))
visreg(rm1, "diameter")
visreg(rm1, "height")
visreg(rm2, "diameter", by="height", overlay=TRUE)
visreg(rm2, "height", by="diameter", overlay=TRUE)
dev.off()

sink("Output/2.model_comparison.txt")
cat("\n\nRootmass models:\n----------------\n### Linear regression ignoring temperature variation:");  summary(rm1)
cat("\n### Linear regression considering interaction with temperature:"); summary(rm2)
cat("### Comparison between both models:\n")
AIC(rm1, rm2); BIC(rm1, rm2)
sink()

# Estimate the WTC3 root mass from the fitted linear regression equation
# Go for the regression with interaction (rm1) that suits the pot experiment (Euc.Teri.) data best
# otherwise we could only use diameter to predict the rootmass (as intercept and height aren't significant)
# eq.rm = function(x,y){exp(coefficients(rm1)[1] + coefficients(rm1)[2] * log(x)  + coefficients(rm1)[3] * log(y))}
eq.rm = function(x,y){exp(coefficients(rm2)[1] + coefficients(rm2)[2] * log(x)  + coefficients(rm2)[3] * log(y)
                          + coefficients(rm2)[4] * log(x) * log(y) )}
height.dia.initial = subset(height.dia.final, Date %in% as.Date("2013-09-17") & case %in% as.factor("3"))
names(height.dia.initial)[9] = "chamber_type"
height.dia.initial$RM = eq.rm(height.dia.initial$diameter, height.dia.initial$height)
height.dia.initial$RM_SE = (((height.dia.initial$diameter_SE)^2 + (height.dia.initial$height_SE)^2 )/2 )^0.5

rootmass.harvest = merge(unique(height.dia[,c("chamber","T_treatment")]), rootmass.harvest[,c("chamber","chamber_type","RootDMtotal")], by=c("chamber"))
rootmass.harvest$Date = as.Date("2014-05-27")
rootmass.harvest <- summaryBy(RootDMtotal ~ Date+T_treatment+chamber_type, data=rootmass.harvest, FUN=c(mean,standard.error))
names(rootmass.harvest)[4:5] = c("RM","RM_SE")

Dates = rep(seq(as.Date("2013-09-17"), as.Date("2014-05-27"), by="days"), 4)
biomass = data.frame(Date = Dates,
                     T_treatment = rep(unique(rootmass.harvest$T_treatment), each=length(Dates)/2))
biomass = biomass[with(biomass, order(T_treatment,Date)), ]
biomass$chamber_type = rep(unique(rootmass.harvest$chamber_type), length(Dates)/2)

biomass = merge(biomass, height.dia.initial[,c("Date","T_treatment","chamber_type","RM","RM_SE")], all = TRUE)
biomass = merge(biomass, rootmass.harvest, all = TRUE)

# biomass = data.frame(chamber = rep(chambers, each = length(unique(height.dia$Date))),
#                      Date = rep(unique(height.dia$Date),length(chambers)),
#                      chamber_type = character(length(chambers) * length(unique(height.dia$Date))),
#                      T_treatment = character(length(chambers) * length(unique(height.dia$Date))),
#                      W_treatment = character(length(chambers) * length(unique(height.dia$Date))), stringsAsFactors=FALSE)
# biomass = subset(biomass, Date >= as.Date("2013-09-14") & Date <= as.Date("2014-05-27"))
# biomass <- merge(biomass, height.dia, by=c("chamber","Date"))


# for(i in 1:length(chambers)) {
#   height.dia.sub = subset(height.dia, chamber %in% as.factor(chambers[i]))
#   # height.dia.sub = subset(height.dia, chamber %in% as.factor(chambers[i]) & Water_treatment %in% as.factor("control"))
#   if (i==11) {
#     height.dia.sub.1 = subset(height.dia.sub, Stem_number %in% 1)
#     height.dia.sub.1.2 = subset(height.dia.sub, Stem_number %in% 1.2)
#     height.dia.sub.1[height.dia.sub.1$DateTime >= as.Date("2013-12-24"),"Plant_height"] = height.dia.sub.1.2[,"Plant_height"]
#     height.dia.sub = height.dia.sub.1
#   }
#   keeps <- c("chamber", "DateTime", "T_treatment", "Plant_height", "X15", "X65")
#   height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
#   
#   D.15 <- lm(X15 ~ X65, data=height.dia.sub)
#   visreg(D.15, "X65", overlay=TRUE)
#   summary(D.15)
#   eq.D = function(x){coefficients(D.15)[1] + coefficients(D.15)[2] * x }
#   
#   height.dia.sub$X15 = eq.D(height.dia.sub$X65)
#   height.dia.sub = subset(height.dia.sub, DateTime >= as.Date("2013-09-14") & DateTime <= as.Date("2014-05-27"))
#   height.dia.sub = height.dia.sub[!is.na(height.dia.sub$X65),]
#   height.dia.sub$rootmass = eq.rm(height.dia.sub$X15, height.dia.sub$Plant_height)
#   keeps <- c("chamber", "DateTime", "T_treatment", "X15", "Plant_height", "rootmass")
#   height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
#   names(height.dia.sub) <- c("chamber","Date","T_treatment","diameter","height","rootmass")
#   height.dia.sub$T_treatment = as.character(height.dia.sub$T_treatment)
#   
#   biomass[(1+(i-1)*length(unique(biomass$Date))) : (i*length(unique(biomass$Date))), 
#           c("T_treatment","diameter","height","rootmass")] = height.dia.sub[,c("T_treatment","diameter","height","rootmass")]
# }
# biomass$T_treatment = as.factor(biomass$T_treatment)
# 
# temp = as.factor(c("ambient","elevated"))
# i=1; j=1
# for(i in 1:length(unique(biomass$T_treatment))) {
#   biomass.idn = subset(biomass,T_treatment==temp[i])
#   for(j in 1:length(unique(biomass.idn$Date))) {
#     biomass.date = subset(biomass.idn, Date == unique(biomass.idn$Date)[j])
#     biomass.date[nrow(biomass.date)+1, c("diameter","height","rootmass")] = colMeans(biomass.date[c("diameter","height","rootmass")], na.rm = TRUE) # R7 = Average of data
#     biomass.date[nrow(biomass.date)+1, c("diameter","height","rootmass")] = (apply(biomass.date[c("diameter","height","rootmass")], 2, sd, na.rm = TRUE))/(nrow(biomass.date)-1)^0.5 # R8 = Standard error of data
#     # height.dia.idn.date[nrow(height.dia.idn.date)+1, 2:ncol(height.dia.idn.date)] = apply(height.dia.idn.date[2:ncol(height.dia.idn.date)], 2, sd, na.rm = TRUE) # R8 = Standard deviation of data
#     biomass.date[(nrow(biomass.date)-1):nrow(biomass.date), c("Date","T_treatment")] = biomass.date[1, c("Date","T_treatment")]
#     dimnames(biomass.date)[[1]] <- c(1:(nrow(biomass.date)-2), "Mean", "SE")
#     if (i == 1 && j == 1) {
#       biomass.final <- biomass.date[0,(2:ncol(biomass.date))]
#     }
#     biomass.final[j+(i-1)*length(unique(biomass.idn$Date)), ] <- biomass.date["Mean", (2:ncol(biomass.date))]
#     biomass.final$diameter_SE[j+(i-1)*length(unique(biomass.idn$Date))] <- biomass.date["SE", 4]
#     biomass.final$height_SE[j+(i-1)*length(unique(biomass.idn$Date))] <- biomass.date["SE", 5]
#     biomass.final$rootmass_SE[j+(i-1)*length(unique(biomass.idn$Date))] <- biomass.date["SE", 6]
#   }
# }
biomass$RM = biomass$RM * 0.48 # unit conversion from gDM to gC
biomass$RM_SE = biomass$RM_SE * 0.48 # unit conversion from gDM to gC










# melted.biomass.final = melt(biomass.final[,c("diameter","height","rootmass","Date","T_treatment")], id.vars=c("Date","T_treatment"))
# melted.biomass.error = melt(biomass.final[,c("diameter_SE","height_SE","rootmass_SE","Date","T_treatment")], id.vars=c("Date","T_treatment"))
# melted.biomass.final$SE = melted.biomass.error$value
# 
# i = 0
# font.size = 12
# plots = list() 
# meas = as.factor(c("diameter","height","rootmass"))
# error = as.factor(c("diameter_SE","height_SE","rootmass_SE"))
# title = as.character(c("A","B","C"))
# pd <- position_dodge(2) # move the overlapped errorbars horizontally
# for (p in 1:length(meas)) {
#   melted.biomass.final.trait = subset(melted.biomass.final,variable %in% meas[p])
#   
#   i = i + 1
#   plots[[i]] = ggplot(melted.biomass.final.trait, aes(x=Date, y=value, group = T_treatment, colour=T_treatment)) + 
#     geom_point(position=pd) +
#     geom_errorbar(position=pd,aes(ymin=value-SE, ymax=value+SE), colour="grey", width=2) +
#     geom_line(position=pd,data = melted.biomass.final.trait, aes(x = Date, y = value, group = T_treatment, colour=T_treatment)) +
#     ylab(paste(as.character(meas[p]),"(g C)")) + xlab("Month") +
#     scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
#     labs(colour="Temperature") +
#     scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
#     theme_bw() +
#     theme(legend.title = element_text(colour="black", size=font.size)) +
#     theme(legend.text = element_text(colour="black", size = font.size)) +
#     theme(legend.position = c(0.9,0.2)) +
#     theme(legend.key = element_blank()) +
#     theme(text = element_text(size=font.size)) +
#     theme(axis.title.x = element_blank()) +
#     theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
#   
#   if (p==1) {
#     plots[[i]] = plots[[i]] + ylab(paste(as.character(meas[p]),"(mm)"))
#   } else if (p==2) {
#     plots[[i]] = plots[[i]] + ylab(paste(as.character(meas[p]),"(cm)"))
#   }
#   if (p>1) {
#     plots[[i]] = plots[[i]] + guides(colour=FALSE)
#   }
# }
# 
# png("output/2.rootmass_height_dia_over_time.png", units="px", width=2200, height=1600, res=220)
# print (do.call(grid.arrange,  plots))
# dev.off()
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Script to setup a model for the leaf mass
harvest.wtc = read.csv("data/WTC_TEMP_CM_HARVEST-CANOPY_20140526-20140528_L1_v1.csv")
keeps <- c("chamber", "T_treatment", "SLA", "StemDM", "TotalLeafDM")
harvest.wtc = harvest.wtc[ , keeps, drop = FALSE]
harvest.wtc <- dplyr::summarize(group_by(harvest.wtc,chamber,T_treatment),
                                SLA = mean(SLA,na.rm=T),
                                StemDM = sum(StemDM,na.rm=T),
                                TotalLeafDM=sum(TotalLeafDM,na.rm=T),
                                Date = max(melted.data$Date))
harvest.wtc$LMA = 1 / (harvest.wtc$SLA / (100*100))
harvest.wtc = data.frame(harvest.wtc)
harvest.wtc$chamber_type = as.factor( ifelse(harvest.wtc$chamber %in% drought.chamb, "drought_treat", "watered_treat") )

harvest.wtc.mean <- summaryBy(SLA+StemDM+TotalLeafDM+LMA ~ Date+T_treatment+chamber_type, data=harvest.wtc, FUN=c(mean,standard.error))
names(harvest.wtc.mean) = c("Date", "T_treatment", "chamber_type","SLA", "SM", "LM", "LMA", "SLA_SE", "SM_SE", "LM_SE", "LMA_SE")

leaf.trait = read.csv("data/WTC_TEMP_CM_TDLSUNSHADE-TRAIT_20131026-20140421.csv")
names(leaf.trait)[4] = "ch"; names(leaf.trait)[6] = "T_treatment"
leaf.trait = merge(data.frame(ch=unique(as.factor(leaf.trait$ch)), chamber= as.factor(harvest.wtc$chamber)), leaf.trait, by="ch", all=TRUE)
leaf.trait$chamber_type = as.factor( ifelse(leaf.trait$chamber %in% drought.chamb, "drought_treat", "watered_treat") )

leaf.trait$Date = paste( leaf.trait$Year , leaf.trait$Month , 15, sep = "-" )
leaf.trait$Date <- as.factor(leaf.trait$Date)
leaf.trait$Date <- parse_date_time(leaf.trait$Date,"y m d")
leaf.trait$Date = as.Date(leaf.trait$Date, format = "%Y-%m-%d")
keeps = c("Date", "chamber", "T_treatment", "chamber_type", "leaf", "lma")
leaf.trait = leaf.trait[ , keeps, drop = FALSE]

# # leaf.trait <- summaryBy(lma~Year+Month+campaign+chamber+T_treatment+chamber_type, data=leaf.trait, FUN=c(mean,length))
# # leaf.trait.1 <- dplyr::summarize(group_by(leaf.trait,Year,Month,chamber,T_treatment),
# #                              lma=mean(lma,na.rm=T))
# leaf.trait.mean <- summaryBy(lma~Year+Month+campaign+T_treatment+chamber_type, data=leaf.trait, FUN=c(mean,standard.error,length))
# names(leaf.trait.mean)[5:6] = c("LMA","LMA_SE")
# 
# leaf.trait <- summaryBy(lma~Year+Month+campaign+chamber+T_treatment+chamber_type, data=leaf.trait, FUN=c(mean,length))
# 
# keeps <- c("Date", "T_treatment", "LMA", "LMA_SE")
# leaf.trait.mean = leaf.trait.mean[ , keeps, drop = FALSE]
# leaf.trait.mean = rbind(leaf.trait.mean, harvest.wtc.mean[ , keeps, drop = FALSE])
# # leaf.trait.mean = data.frame(leaf.trait.mean)

# plots = list()
# # plots[[1]] = plot[[3]]
# pd <- position_dodge(3)
# plots[[1]] = ggplot(leaf.trait.mean, aes(x=Date, y=LMA, group = T_treatment, colour=T_treatment)) + 
#   geom_point(position=pd) +
#   geom_errorbar(position=pd,aes(ymin=LMA-LMA_SE, ymax=LMA+LMA_SE), colour="grey", width=2) +
#   geom_line(position=pd,data = leaf.trait.mean, aes(x = Date, y = LMA, group = T_treatment, colour=T_treatment)) +
#   ylab(expression("LMA"~"(g "*m^"-2"*")")) +
#   scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(melted.data.sub$Date), max(melted.data.sub$Date)+1)) +
#   labs(colour="Temperature") +
#   scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
#   theme_bw() +
#   theme(legend.title = element_text(colour="black", size=font.size)) +
#   theme(legend.text = element_text(colour="black", size = font.size)) +
#   theme(legend.position = c(0.15,0.8)) +
#   theme(legend.key = element_blank()) +
#   theme(text = element_text(size=font.size)) +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
# 
# png("output/3.LMA_over_time.png", units="px", width=2200, height=1600, res=220)
# print (do.call(grid.arrange,  plots))
# dev.off()
leaf.trait = leaf.trait[,-2]
leaf.trait.sub <- summaryBy(lma~Date+T_treatment+chamber_type+leaf, data=leaf.trait, FUN=c(mean,standard.error,length))
names(leaf.trait.sub)[5:6] = c("LMA","LMA_SE")

plots = list() 
pd <- position_dodge(3)
plots[[1]] = ggplot(data = leaf.trait.sub, aes(x=Date, y=LMA, group=leaf, colour=as.factor(leaf))) + 
  geom_point(position=pd) + 
  geom_errorbar(position=pd,aes(ymin=LMA-LMA_SE, ymax=LMA+LMA_SE), colour="grey", width=2) +
  geom_line(position=pd, data = leaf.trait.sub, aes(x=Date, y=LMA, group=leaf, colour=as.factor(leaf))) +
  labs(colour="Treatments", shape="Leaf type", linetype="Leaf type") +
  facet_wrap(T_treatment ~ chamber_type, scales="free_x") +
  ylab(expression("LMA"~"(g "*m^"-2"*")")) + theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.15,0.94)) + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

leaf.trait.sub1 <- summaryBy(lma~Date+T_treatment+chamber_type, data=leaf.trait, FUN=c(mean,standard.error,length))
names(leaf.trait.sub1)[4:5] = c("LMA","LMA_SE")
plots[[2]] = ggplot(data = leaf.trait.sub1, aes(x=Date, y=LMA, group=interaction(T_treatment,chamber_type), colour=as.factor(T_treatment), shape=as.factor(chamber_type))) + 
  geom_point(position=pd) + 
  geom_errorbar(position=pd,aes(ymin=LMA-LMA_SE, ymax=LMA+LMA_SE), colour="grey", width=2) +
  geom_line(position=pd, data = leaf.trait.sub1, aes(x=Date, y=LMA, group=interaction(T_treatment,chamber_type), colour=as.factor(T_treatment), linetype=as.factor(chamber_type))) + 
  labs(colour="Treatments", shape="Chamber type", linetype="Chamber type") +
  ylab(expression("LMA"~"(g "*m^"-2"*")")) + theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.5,0.9), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# png("output/3.LMA_over_time.png", units="px", width=2200, height=1600, res=220)
pdf(file = "output/3.LMA_over_time.pdf")
# lay <- rbind(c(1,1),c(2,3))
# grid.arrange(grobs = plots[c(5,1,2)], layout_matrix = lay)
plots
dev.off()

#----------------------------------------------------------------------------------------------------------------
# calculate daily leaf mass
LM.daily = data.frame(Date = melted.data.sub$Date, T_treatment = melted.data.sub$T_treatment, LA = melted.data.sub$value, LA_SE = melted.data.sub$SE)
LM.daily = merge(LM.daily, leaf.trait.mean, by = c("Date", "T_treatment"), all = TRUE)
LM.daily = LM.daily[with(LM.daily, order(T_treatment,Date)), ]

# interpolate daily LMA from monthly data
for(i in 1:length(temp)) {
  LM.daily.sub = subset(LM.daily, T_treatment %in% temp[i])
  LM.daily$LMA[(1+(i-1)*nrow(LM.daily.sub)) : (i*nrow(LM.daily.sub))] = na.fill(na.approx(LM.daily.sub$LMA), "extend")
  
  ptsLin   <- approx(LM.daily.sub$Date, LM.daily.sub$LMA, method="linear", n=(max(LM.daily$Date)-min(LM.daily$Date)))
  LM.daily$LMA[(1+(i-1)*nrow(LM.daily.sub)) : (i*nrow(LM.daily.sub))] = c(ptsLin$y[1],ptsLin$y)
  ptsLin   <- approx(LM.daily.sub$Date, LM.daily.sub$LMA_SE, method="linear", n=(max(LM.daily$Date)-min(LM.daily$Date)))
  LM.daily$LMA_SE[(1+(i-1)*nrow(LM.daily.sub)) : (i*nrow(LM.daily.sub))] = c(ptsLin$y[1],ptsLin$y)
}

LM.daily$LM = (LM.daily$LA * LM.daily$LMA) / 1000 
# LM.daily$LM_SE = (((LM.daily$LA_SE * (6)^0.5)^2/6 + (LM.daily$LMA_SE * (6)^0.5)^2/6)  )^0.5 
LM.daily$LM_SE = (((LM.daily$LA_SE * (6)^0.5)^2 + (LM.daily$LMA_SE * (6)^0.5)^2)/2  )^0.5 / 1000

# Daily LMA plot with SE
plots[[2]] = ggplot(LM.daily, aes(x=Date, y=LMA, group = T_treatment, colour=T_treatment)) + 
  # geom_point(position=pd) +
  geom_ribbon(data = LM.daily, aes(ymin=LMA-LMA_SE, ymax=LMA+LMA_SE), linetype=2, alpha=0.1,size=0.1) +
  # geom_errorbar(position=pd,aes(ymin=LMA-LMA_SE, ymax=LMA+LMA_SE), colour="grey", width=2) +
  geom_line(data = LM.daily, aes(x = Date, y = LMA, group = T_treatment, colour=T_treatment)) +
  ylab(expression("LMA"~"(g "*m^"-2"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(melted.data.sub$Date), max(melted.data.sub$Date)+1)) +
  labs(colour="Temperature") +
  scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.15,0.8)) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

png("output/3.LMA_data_vs_prediction.png", units="px", width=2200, height=1600, res=220)
print (do.call(grid.arrange,  plots))
dev.off()

# Daily LM plot with SE
plots[[3]] = ggplot(LM.daily, aes(x=Date, y=LM, group = T_treatment, colour=T_treatment)) + 
  # geom_point(position=pd) +
  geom_ribbon(data = LM.daily, aes(ymin=LM-LM_SE, ymax=LM+LM_SE), linetype=2, alpha=0.1,size=0.1) +
  # geom_errorbar(position=pd,aes(ymin=LMA-LMA_SE, ymax=LMA+LMA_SE), colour="grey", width=2) +
  geom_line(data = LM.daily, aes(x = Date, y = LM, group = T_treatment, colour=T_treatment)) +
  ylab(expression("Leaf mass"~"(kg "*plant^"-2"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(melted.data.sub$Date), max(melted.data.sub$Date)+1)) +
  labs(colour="Temperature") +
  scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.15,0.8)) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("output/3.LA_LMA_LM_estimated_daily.png", units="px", width=2200, height=1600, res=220)
print (do.call(grid.arrange,  plots))
dev.off()
#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
# Script to setup a model for the leaf mass
litterfall = read.csv("data/WTC_TEMP_CM_LEAFLITTER_20130913-20140528_L1.csv")
litterfall$startDate = as.Date(litterfall$startDate)
litterfall$collectionDate = as.Date(litterfall$collectionDate)
litterfall$Date <- (litterfall$startDate + ((litterfall$collectionDate - litterfall$startDate) / 2))
litterfall = subset(litterfall, Date >= as.Date("2013-09-14") & Date <= as.Date("2014-05-27"))

# convert to data.table in place
setDT(litterfall)
# dcast and do individual sums
litterfall.cast = dcast.data.table(litterfall, chamber ~ Date, value.var = 'litter', fun.aggregate = sum)
# cumsum
litterfall.cum <- litterfall.cast[, as.list(cumsum(unlist(.SD))), by = chamber]

litterfall.cum.melt <- melt(litterfall.cum, id.vars = "chamber")

litterfall.cum.melt = merge(litterfall.cum.melt, harvest.wtc[,c("chamber","T_treatment")])
names(litterfall.cum.melt)[2:3] = c("Date","litter")
litterfall.cum.melt$Date = as.Date(litterfall.cum.melt$Date)
litterfall.cum.melt = summaryBy(litter ~ Date+T_treatment, data=litterfall.cum.melt, FUN=c(mean,standard.error))
names(litterfall.cum.melt)[3:4] = c("litter","litter_SE")

pd <- position_dodge(3)
plots[[4]] = ggplot(litterfall.cum.melt, aes(x=Date, y=litter, group = T_treatment, colour=T_treatment)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd,aes(ymin=litter-litter_SE, ymax=litter+litter_SE), colour="grey", width=2) +
  geom_line(position=pd,data = litterfall.cum.melt, aes(x = Date, y = litter, group = T_treatment, colour=T_treatment)) +
  ylab(expression("Leaf Litter"~"(g "*plant^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(melted.data.sub$Date), max(melted.data.sub$Date)+1)) +
  labs(colour="Temperature") +
  scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.15,0.8)) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

png("output/4.litterfall.png", units="px", width=2200, height=1600, res=220)
plots[[4]]
dev.off()

png("output/3.LA_LMA_LM_litter_estimated_daily.png", units="px", width=2200, height=1600, res=220)
print (do.call(grid.arrange,  plots))
dev.off()
#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------





