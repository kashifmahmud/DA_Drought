
# Estimate one-point Vcmax from data on A, Ci & Tleaf
# Function call: Vcmax1pt(Photo,Ci,Tleaf,Tcorr)
# Function needs input of photosynthesis, Ci, and Tleaf
# If Tcorr = T (default) it returns a value corrected to 25 deg 
# Otherwise value is at measurement temperature

# Functions
.Rgas <- function()8.314
Tk <- function(x)x+273.15
# Arrhenius
arrh <- function(Tleaf, Ea){
  exp((Ea * (Tk(Tleaf) - 298.15)) / (298.15 * .Rgas() * Tk(Tleaf))) 
}
#Gamma Star
TGammaStar <- function(Tleaf,  Egamma=37830.0, value25=42.75){  
  value25*arrh(Tleaf,Egamma)
}
# Km 
TKm <- function(Tleaf) {
  Kc <- 404.9*arrh(Tleaf,79430)
  Ko <- 278400*arrh(Tleaf,36380)
  Oi <- 205000
  Km <- Kc * (1+Oi/Ko)
  return(Km)
}
# Vcmax
peaked_arrh <- function(Tleaf, Ea, delS) {
  TleafK <- Tk(Tleaf)
  Ed <- 2E5
  func <- exp((Ea*(TleafK - 298.15))/(298.15 * .Rgas() * TleafK)) * 
    (1 + exp((298.15*delS - Ed)/(298.15 * .Rgas()))) / 
    (1 + exp((TleafK*delS - Ed)/(TleafK * .Rgas())))
  return(func)
}

# Main function
Vcmax1pt <- function(Asat,Tcorr) {
  Asat.final = Asat
  Photo = Asat$Photo
  Ci = Asat$Ci
  Tleaf = Asat$Tleaf
  Rday = Asat$Rd_pred
  # RVrat <- 0.015 # assume Rd/Vcmax ratio equals 0.015
  Gstar <- TGammaStar(Tleaf)
  Km <- TKm(Tleaf)
  Vcmax <- (Photo - Rday) / ((Ci-Gstar)/(Ci+Km)) 
  Vcmax25 <- Vcmax / peaked_arrh(Tleaf, Ea = 67338, delS = 631.7)
  if (Tcorr == T) {
    Asat.final$Vcmax25 = Vcmax25
  } else {
    Asat.final$Vcmax = Vcmax
  }
  return(Asat.final)
}

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# Read in Asat spot measurements
Asat <- read.csv("raw_data/GHS30_Eglob-TxCxW_GEasat_20110117-20110321_L1.csv")
# Asat = subset(Asat, Cond >= 0) # Filter Licor data having conductance less than 0.05
Asat = subset(Asat, Cond >= 0 & Photo >=0) # Filter Licor data having conductance less than 0.05

Asat$Date <- as.Date(Asat$Date, format="%Y/%m/%d")
# resp_dark$ID <- paste(resp_dark$plot, resp_dark$pot, sep = "-")

# leafmass <- read.csv("raw data/seedling leaf mass area.csv")

Asat <- (subset(Asat, select = c("Date", "Potnum", "Temp", "CO2", "Water", "Photo", "Cond", "Ci", "VpdL", "Tleaf", "CO2S")))
Asat <- subset(Asat, Temp == "Amb" & CO2 == 400) # only consider the ambient drought treatment

# Asat = summaryBy(Photo+Cond+Ci+VpdL+Tleaf+CO2S ~ Date+Water, data=Asat, FUN=c(mean,standard.error))
# names(Asat) = c("Date", "Water", "Photo", "Cond", "Ci", "VpdL", "Tleaf", "CO2S", "Photo_se", "Cond_se", "Ci_se", "VpdL_se", "Tleaf_se", "CO2S_se")
# keeps = c("Date", "Water", "Photo", "Cond", "Ci", "VpdL", "Tleaf", "CO2S")
# Asat = Asat[ , keeps, drop = FALSE]

# boxplot(Photo~ Water, data=rdark)

#---------------------------------------------------------------------------------------------
# Plot Vcmax25 over time and treatment
Asat.df = summaryBy(Photo ~ Date+Water, data=Asat, FUN=c(mean,standard.error))
names(Asat.df) = c("Date", "Water", "Photo", "Photo_se")

font.size = 12
pd <- position_dodge(2)
p0 = ggplot(Asat.df, aes(x=Date, y=Photo, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=Photo-Photo_se, ymax=Photo+Photo_se), colour="grey", width=3) +
  geom_line(position=pd, data = Asat.df, aes(x = Date, y = Photo, group = Water, colour=Water)) +
  ylab(expression(A[sat] ~ (mu~mol ~ m^{-2} ~ s^{-1}))) +
  # scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Treatment") +
  # scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.75,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

png(file = "output/1.Asat.png")
print (p0)
dev.off()

#-------------------------------------------------------------------------------------
# process met data
tempdata = read.csv("raw_data/GHS30_Eglob-TxCxW_ClimateTemp_20110117-20110321_L1.csv")
tempdata$Date <- as.Date(tempdata$Date, format="%d/%m/%Y")
tempdata$Time <- times(tempdata$Time)
tempdata = tempdata[1:3]
names(tempdata)[2:3] <- c("time","temp")

RHdata =read.csv("raw_data/GHS30_Eglob-TxCxW_ClimateRH_20110117-20110321_L1.csv")
RHdata$Date <- as.Date(RHdata$Date, format="%d/%m/%Y")
RHdata$Time <- times(RHdata$Time)
RHdata = RHdata[1:3]
names(RHdata)[2:3] <- c("time","RH")

pardata.raw = read.csv("raw_data/RefMinData Query1.txt")
pardata.raw$DateTime = as.POSIXct(pardata.raw$DateTime, format="%d/%m/%Y %H:%M:%S")
pardata.raw$Date <- as.Date(pardata.raw$DateTime, format="%d/%m/%Y %H:%M:%S")
pardata.raw = subset(pardata.raw, Date >= as.Date("2011-01-16") & Date <= as.Date("2011-03-22"))

pardata.raw$DateTime15 = cut(pardata.raw$DateTime, breaks="15 min")
pardata.raw$DateTime = NULL
pardata = aggregate(PAR ~ DateTime15, FUN=mean, data=pardata.raw)
pardata$Date <- as.Date(pardata$DateTime15)
pardata = subset(pardata, Date >= as.Date("2011-01-17") & Date <= as.Date("2011-03-21"))

# need to turn the datetime 15 into hms
pardata$time <- format(as.POSIXct(pardata$DateTime15,format="%Y-%m-%d %H:%M:%S"),"%H:%M:%S")
pardata$DateTime15 = NULL

# Merge all met data into one dataframe
metdata = merge(tempdata,RHdata, by=c("Date","time"))
metdata = merge(metdata,pardata, by=c("Date","time"))

# convert RH to VPD
metdata$VPD <- RHtoVPD(metdata$RH, metdata$temp, Pa=101)

#-------------------------------------------------------------------------------------
# merge ps data to met data
Asat <- merge(metdata, Asat, by="Date")

# Rdark Q10 equations by treatment
rdarkq10 <- read.csv("processed_data/rdarkq10.csv")
#need to calculate Rdark through time using rdarkq10 equation by volume
Asat <- merge(Asat, rdarkq10, by=c("Water","Date"))

q10_crous <- 1.95
q25_drake <- 1.86
#refit RD
Asat$Rd_pred <- with(Asat, rd18.5 * q25_drake^((temp-18.5)/10))

Asat$period <- ifelse(Asat$PAR > 2, "Day", "Night")
Asat$Rd_pred <- ifelse(Asat$period=="Day", 0.7*Asat$Rd_pred, Asat$Rd_pred) # If it's day, subtract 30% from the leaf R fraction

# Calculate Vcmax
Asat.day = subset(Asat, period == "Day") # consider only day values
Asat.final = Vcmax1pt(Asat.day,Tcorr=T) # function to calculate Vcmax

Vcmax25.df = summaryBy(Vcmax25 ~ Date+Water, data=Asat.final, FUN=c(mean,standard.error))
names(Vcmax25.df) = c("Date", "Water", "Vcmax25", "Vcmax25_se")

# Test whether there are treatment and temporal effect on Vcmax25
Vcmax25.anova <- lm(Vcmax25 ~ Water+Date, data = Vcmax25.df)
anova(Vcmax25.anova)
# Neither treatment effect nor temporal effect on Vcmax25

#---------------------------------------------------------------------------------------------
# Plot Vcmax25 over time and treatment
font.size = 12
pd <- position_dodge(2)
p2 = ggplot(Vcmax25.df, aes(x=Date, y=Vcmax25, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=Vcmax25-Vcmax25_se, ymax=Vcmax25+Vcmax25_se), colour="grey", width=3) +
  geom_line(position=pd, data = Vcmax25.df, aes(x = Date, y = Vcmax25, group = Water, colour=Water)) +
  ylab(expression(V[cmax] ~ (mu~mol ~ m^{-2} ~ s^{-1}))) +
  # scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Treatment") +
  # scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.75,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

png(file = "output/2.Vcmax.png")
print (p2)
dev.off()
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------


