
# calculate g1 and compare with SWC

#-------------------------------------------------------------------------------------
# Read in Asat spot measurements
Asat <- read.csv("raw_data/GHS30_Eglob-TxCxW_GEasat_20110117-20110321_L1.csv")
# Asat = subset(Asat, Cond >= 0) # Filter Licor data having conductance less than 0.05
Asat = subset(Asat, Cond >= 0 & Photo >=0) # Filter Licor data having conductance >= 0 and Photo >= 0

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

#-------------------------------------------------------------------------------------
Asat$g1 = (((Asat$Cond * Asat$CO2S) / (1.6 * Asat$Photo)) - 1) * (Asat$VpdL)^0.5
# g1.summary = summaryBy(g1 ~ Date+Water, data=Asat.final, FUN=c(mean,standard.error))
# boxplot(g1.mean~ Water, data=g1.summary)
g1.df = summaryBy(g1 ~ Date+Water, data=Asat, FUN=c(mean,standard.error))
names(g1.df) = c("Date", "Water", "g1", "g1_se")
# g1.df = subset(g1.df, Water %in% as.factor(c("Rewatered  drought", "Sustained drought",  "Well watered")))

# Test whether there are treatment and temporal effect on g1
g1.anova <- lm(g1 ~ Water+Date, data = g1.df)
anova(g1.anova)
# No treatment effect on g1, might have slight temporal effect though
g1.anova <- lm(g1 ~ Date, data = g1.df)
anova(g1.anova)

#---------------------------------------------------------------------------------------------
# Plot g1 over time and treatment
font.size = 12
pd <- position_dodge(2)
p3 = ggplot(g1.df, aes(x=Date, y=g1, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=g1-g1_se, ymax=g1+g1_se), colour="grey", width=3) +
  geom_line(position=pd, data = g1.df, aes(x = Date, y = g1, group = Water, colour=Water)) +
  ylab(expression(g1 ~ (mmol ~ m^{-2} ~ s^{-1}))) +
  # scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Treatment") +
  # scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.8,0.6), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

png(file = "output/3.g1.png")
print (p3)
dev.off()

#---------------------------------------------------------------------------------------------
# Compare g1 with SWC
# Read in SWC measurements
swc <- read.csv("raw_data/GHS30_Eglob-TxCxW_SWC,SWP_20110117-20110321_L1.csv")
swc$Date <- as.Date(swc$Date,format="%d/%m/%Y")
swc <- subset(swc, Temp == "Amb" & CO2 == 400) # only consider the ambient drought treatment
swc = subset(swc, Date >= as.Date("2011-01-17") & Date <= as.Date("2011-03-21"))
swc$Time <- format(as.POSIXct(swc$Time,format="%H:%M:%S"),"%H:%M:%S")
swc = summaryBy(VWC ~ Date+Water, data=swc, FUN=c(mean,standard.error))
names(swc) = c("Date", "Water", "VWC", "VWC_se")
swc$treatment[swc$Water %in% as.factor("Rewatered Drought")] = "Rewatered  drought"
swc$treatment[swc$Water %in% as.factor("Sustained Drought")] = "Sustained drought"
swc$treatment[swc$Water %in% as.factor("Well Watered")] = "Well watered"
swc$Water = NULL
colnames(swc)[which(names(swc) == "treatment")] <- "Water"

g1.swc = merge(g1.df, swc, by=c("Date","Water"))
font.size = 12
pd <- position_dodge(0)
p4 = ggplot(g1.swc, aes(x=VWC, y=g1, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  # geom_errorbar(position=pd, aes(ymin=g1-g1_se, ymax=g1+g1_se), colour="grey", width=3) +
  geom_line(position=pd, data = g1.swc, aes(x = VWC, y = g1, group = Water, colour=Water)) +
  ylab(expression(g1 ~ (mmol ~ m^{-2} ~ s^{-1}))) +
  xlab("VWC (%)") +
  # scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Treatment") +
  # scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.8,0.6), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_text(size = font.size, vjust=0.3)) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

png(file = "output/4.g1_vs_VWC.png")
print (p4)
dev.off()

#---------------------------------------------------------------------------------------------

