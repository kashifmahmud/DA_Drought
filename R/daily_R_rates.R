# Calculate daily respiration rates for roots, wood and foliage

# Calculate daily below ground root respiration rates
# Coarse root respiration rates
branch.resp = read.csv("raw_data/WTC_TEMP_CM_GX-RBRANCH_20140513-20140522_L1_v1.csv") # Rbranch: branch respiration (nmol CO2 g-1 s-1) 
branch.resp$date = as.Date(branch.resp$date)
branch.resp = subset(branch.resp, date %in% as.Date("2014-05-13")) # Only consider the pre-girdling data
branch.resp$chamber_type = as.factor( ifelse(branch.resp$chamber %in% drought.chamb, "drought", "watered") )
branch.resp$Treatment <- as.factor(paste(branch.resp$T_treatment, branch.resp$chamber_type))

# # Test for any significant difference between the treatment groups
# boxplot(branch.resp$Rbranch ~ branch.resp$Treatment, xlab="Treatment", ylab=(expression("Branch wood respiration"~"(nmol CO2 "*g^"-1"*" "*s^"-1"*")")))

summary(aov(Rbranch ~ T_treatment * chamber_type, data = branch.resp)) # YES, there is significant difference only accross the temperature treatments
# summary(aov(Rbranch ~ Treatment, data = branch.resp)) # YES, there is significant difference accross the treatments
# t.test(branch.resp$Rbranch ~ branch.resp$T_treatment) # YES, there is significant difference accross temperatures
# t.test(branch.resp$Rbranch ~ branch.resp$chamber_type) # NO, there is no significant difference accross drought/watered treatments

############ So how to group the treatments????????
rd15.root <- summaryBy(Rbranch ~ T_treatment, data=branch.resp, FUN=c(mean))
names(rd15.root)[ncol(rd15.root)] = c("rd15.coarseroot")
rd15.root$rd15.coarseroot = rd15.root$rd15.coarseroot * (10^-9 * 12) * (3600 * 24) # unit conversion from nmolCO2 g-1 s-1 to gC gC-1 d-1
# rd15.root$rd15.coarseroot = rd15.root$rd15.coarseroot * (10^-9 * 12) * (3600 * 24) * (1/c1) # unit conversion from nmolCO2 g-1 s-1 to gC gC-1 d-1
# rd15.coarseroot$rd15.coarseroot_SE = rd15.coarseroot$rd15.coarseroot_SE * (10^-9 * 12) * (3600 * 24) * (1/c1) # unit conversion from nmolCO2 g-1 s-1 to gC gC-1 d-1

# Bole and big tap root respiration rates
bole.resp = read.csv("raw_data/WTC_TEMP_CM_WTCFLUX-STEM_20140528_L1_v1.csv") # Bole root respiration (nmol CO2 g-1 s-1) 
bole.resp$chamber_type = as.factor( ifelse(bole.resp$chamber %in% drought.chamb, "drought", "watered") )
bole.resp$Treatment <- as.factor(paste(bole.resp$T_treatment, bole.resp$chamber_type))

# # Test for any significant difference between the treatment groups
# boxplot(bole.resp$R_stem_nmol ~ bole.resp$Treatment, xlab="Treatment", ylab=(expression("Bole wood respiration"~"(nmol CO2 "*g^"-1"*" "*s^"-1"*")")))

summary(aov(R_stem_nmol ~ T_treatment * chamber_type, data = bole.resp)) # NO, there is no significant difference accross the treatments
# summary(aov(R_stem_nmol ~ Treatment, data = bole.resp)) # NO, there is no significant difference accross the treatments
# t.test(bole.resp$R_stem_nmol ~ bole.resp$T_treatment) # NO, there is no significant difference accross tepmeratures
# t.test(bole.resp$R_stem_nmol ~ bole.resp$chamber_type) # NO, there is no significant difference accross drought/watered treatments

rd15.root$rd15.boleroot = mean(bole.resp$R_stem_nmol)
rd15.root$rd15.boleroot = rd15.root$rd15.boleroot * (10^-9 * 12) * (3600 * 24) # unit conversion from nmolCO2 g-1 s-1 to gC gC-1 d-1
# rd15.root$rd15.boleroot = rd15.root$rd15.boleroot * (10^-9 * 12) * (3600 * 24) * (1/c1) # unit conversion from nmolCO2 g-1 s-1 to gC gC-1 d-1

# Fine root respiration rates (Constant)
# Fine root respiration rate = 10 nmolCO2 g-1 s-1 (Ref: Drake et al. 2017: GREAT exp data; Mark's Email)
rd25.fineroot = 10 * (10^-9 * 12) * (3600 * 24) # unit conversion from nmolCO2 g-1 s-1 to gC gC-1 d-1
# rd25.fineroot = 10 * (10^-9 * 12) * (3600 * 24) * (1/c1) # unit conversion from nmolCO2 g-1 s-1 to gC gC-1 d-1
rd15.root$rd15.fineroot = rd25.fineroot * q10^((15-25)/10)

# Intermediate root respiration rates
rd15.root$rd15.intermediateroot = exp ((log(rd15.root$rd15.coarseroot) + log(rd15.root$rd15.fineroot))/2 ) # unit = gC gC-1 d-1


#----------------------------------------------------------------------------------------------------------------
# import site weather data, take only soil temperatures at 10 cm depth, format date stuff
files <- list.files(path = "raw_data/WTC_TEMP_CM_WTCMET", pattern = ".csv", full.names = TRUE)
temp <- lapply(files, fread, sep=",")
met.data <- rbindlist( temp )

met.data <- met.data[ , c("chamber","DateTime","Tair_al","SoilTemp_Avg.1.","SoilTemp_Avg.2.","PPFD_Avg")]
met.data$SoilTemp <- rowMeans(met.data[,c("SoilTemp_Avg.1.","SoilTemp_Avg.2.")], na.rm=TRUE)
met.data$Date <- as.Date(met.data$DateTime)

# need to turn the datetime into hms
met.data$DateTime <- ymd_hms(met.data$DateTime)
met.data$time <- format(met.data$DateTime, format='%H:%M:%S')

# subset by Date range of experiment
met.data <- subset(met.data[, c("chamber","Date","time","Tair_al","SoilTemp","PPFD_Avg")], Date  >= "2013-09-17" & Date  <= "2014-05-26")
met.data$chamber = as.factor(met.data$chamber)
met.data = merge(met.data, unique(height.dia[,c("chamber","T_treatment")]), by="chamber")

met.data$period <- ifelse(met.data$PPFD_Avg>2,"Day","Night")

# Remove the data with missing air and soil temperatures from met data
# met.data = met.data[complete.cases(met.data$SoilTemp),]
# met.data = met.data[complete.cases(met.data$Tair_al),]

# met.data.na1 = met.data[is.na(met.data$SoilTemp),] # Check any NA values for soil temperature
# met.data.na2 = met.data[is.na(met.data$Tair_al),] # Check any NA values for air temperature
# met.data[Date == as.Date("2013-10-06")]

# need to calculate Rdark through time using rdarkq10 equation by treatment
met.data <- merge(met.data, rd15.root, by=c("T_treatment"))
met.data <- merge(met.data, Tair.final[,c("Date","T_treatment","rd25.foliage")], by=c("Date","T_treatment"))

met.data[,c("Rd.fineroot","Rd.intermediateroot","Rd.coarseroot","Rd.boleroot")] = 
  with(met.data, met.data[,c("rd15.fineroot","rd15.intermediateroot","rd15.coarseroot","rd15.boleroot")] * 
         q10^((SoilTemp-15)/10)) # unit (gC per gC root per day)

# calculate daily stem and branch respiration rates
met.data[,c("Rd.stem","Rd.branch")] = 
  with(met.data, met.data[,c("rd15.boleroot","rd15.coarseroot")] * q10^((Tair_al-15)/10)) # unit (gC per gC wood per day)

# calculate foliage respiration rates in 15-mins interval
met.data[,"Rd.foliage"] = with(met.data, met.data[,"rd25.foliage"] * q10^((Tair_al-25)/10)) # unit (gC per gC foliage per day)
met.data[met.data$period == "Day","Rd.foliage"] = 0.7*met.data[met.data$period == "Day","Rd.foliage"] # 30% reduction in leaf respiration during day time

# Calculate daily mean respiration rates for all tree components by summing all 15-mins data for each day
Rd <- summaryBy(Rd.foliage+Rd.stem+Rd.branch+Rd.fineroot+Rd.intermediateroot+Rd.coarseroot+Rd.boleroot ~ Date+T_treatment, data=met.data, FUN=mean, na.rm=TRUE) # Sum of all same day Rd
# colSums(is.na(Rd)) # Check any NA values for Rd
# Rd.na = Rd[is.na(Rd$Rd.foliage.mean),]

# Fill missing values due to atmospheric data gaps
Rd.sub1 = subset(Rd, T_treatment %in% as.factor("ambient"))
for (i in 3:ncol(Rd.sub1)) {
  Rd.sub1[,i] = na.approx(Rd.sub1[,..i])
}
Rd.sub2 = subset(Rd, T_treatment %in% as.factor("elevated"))
for (i in 3:ncol(Rd.sub2)) {
  Rd.sub2[,i] = na.approx(Rd.sub2[,..i])
}
Rd = rbind(Rd.sub1,Rd.sub2)
# names(Rd)[3:ncol(Rd)] = c("Rd.foliage","Rd.stem","Rd.branch","Rd.fineroot","Rd.intermediateroot","Rd.coarseroot","Rd.boleroot")
# colSums(is.na(Rd.fill)) # Check any NA values for Rd

# Merge respiration rates with daily woodmass, rootmass partitioning, GPP, Ra, LA, mass pool data
# Rd = subset(Rd, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-02-12"))
data.all = merge(data.all, Rd, by=c("Date","T_treatment"), all=TRUE)
# data.all$Treatment <- as.factor(paste(data.all$T_treatment, data.all$chamber_type))
data.all$Treatment <- as.factor(paste(data.all$T_treatment))

# #----------------------------------------------------------------------------------------------------------------
# # write csv file with daily respiration rates for roots, wood and foliage
# write.csv(Rd, "processed_data/Rd.csv", row.names=FALSE) # unit: gC per gC plant per day

#----------------------------------------------------------------------------------------------------------------
# Plot daily respiration rates for roots, wood and foliage
Rd.melt <- melt(Rd, id.vars = c("Date","T_treatment"))
i = 0
font.size = 10
plot = list() 
meas = as.factor(c("Rd.foliage.mean","Rd.stem.mean","Rd.branch.mean","Rd.fineroot.mean","Rd.intermediateroot.mean","Rd.coarseroot.mean","Rd.boleroot.mean"))
# title = as.character(c("A","B","C","D"))
pd <- position_dodge(0) # move the overlapped errorbars horizontally
for (p in 1:length(meas)) {
  Rd.melt.sub = subset(Rd.melt,variable %in% meas[p])
  
  i = i + 1
  plot[[i]] = ggplot(Rd.melt.sub, aes(x=Date, y=value, group = T_treatment, colour=T_treatment)) + 
    geom_point(position=pd) +
    geom_line(position=pd,data = Rd.melt.sub, aes(x = Date, y = value, group = T_treatment, colour=T_treatment)) + 
    ylab(expression(R[foliage]~"(g C "*g^"-1"*" C "*d^"-1"*")")) + xlab("") +
    scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(Rd$Date)-2, max(Rd$Date)+2)) +
    labs(colour="Temperature") +
    scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
    theme_bw() +
    theme(legend.title = element_text(colour="black", size=font.size)) +
    theme(legend.text = element_text(colour="black", size=font.size)) +
    theme(legend.position = c(0.9,0.75), legend.box = "horizontal") + theme(legend.key.height=unit(0.9,"line")) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  if (p==2) {
    plot[[i]] = plot[[i]] + ylab(expression(R[bolewood]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
    # plot[[i]] = plot[[i]] + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 0.75), units="line"))
  }
  if (p==3) {
    plot[[i]] = plot[[i]] + ylab(expression(R[branchwood]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
    # plot[[i]] = plot[[i]] + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 0.75), units="line"))
  }
  if (p==4) {
    plot[[i]] = plot[[i]] + ylab(expression(R[fineroot]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
    # plot[[i]] = plot[[i]] + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 0.75), units="line"))
  }
  if (p==5) {
    plot[[i]] = plot[[i]] + ylab(expression(R[intermediateroot]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
    # plot[[i]] = plot[[i]] + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 0.75), units="line"))
  } 
  if (p==6) {
    plot[[i]] = plot[[i]] + ylab(expression(R[coarseroot]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
  } 
  if (p==7) {
    plot[[i]] = plot[[i]] + ylab(expression(R[boleroot]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
  }
}

png("output/3.Rd.png", units="px", width=2200, height=1500, res=130)
# do.call(grid.arrange,  plot)
ncols = 1
do.call("grid.arrange", c(plot, ncol=ncols))
dev.off()

do.call("grid.arrange", c(plot, ncol=ncols))

