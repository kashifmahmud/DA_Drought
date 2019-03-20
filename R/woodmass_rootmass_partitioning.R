# Estimate the partitioning of woodmass to stem and branches   

treeMass$BMratio = treeMass$branchMass / treeMass$woodMass
treeMass$SMratio = treeMass$boleMass / treeMass$woodMass
treeMass$Treatment <- as.factor(paste(treeMass$T_treatment, treeMass$chamber_type))

# Test for any significant difference between the treatment groups
summary(aov(BMratio ~ Treatment, data = treeMass)) # NO, there is no significant difference accross the treatments
summary(aov(SMratio ~ Treatment, data = treeMass)) # NO, there is no significant difference accross the treatments

woodmass.ratio = data.frame(matrix(ncol = 1, nrow = 252))
names(woodmass.ratio) = c("Date")
woodmass.ratio$Date = as.Date(as.Date("2013-09-17"):as.Date("2014-05-26"))

woodmass.ratio.data = summaryBy(BMratio+SMratio ~ Date, data=treeMass, FUN=c(mean,standard.error))
names(woodmass.ratio.data)[2:ncol(woodmass.ratio.data)] = c("BMratio","SMratio","BMratio_SE","SMratio_SE")

woodmass.ratio = merge(woodmass.ratio, woodmass.ratio.data, by="Date", all = TRUE)
woodmass.ratio = na.interpolation(woodmass.ratio, option = "linear")

#----------------------------------------------------------------------------------------------------------------
# Plot woodmass ratios
woodmass.ratio.melt = melt(woodmass.ratio[,c("Date","BMratio","SMratio")], id.vars = "Date")
woodmass.ratio.melt$biomass = as.factor("Wood")
woodmass.ratio.melt$SE = melt(woodmass.ratio[,c("Date","BMratio_SE","SMratio_SE")], id.vars = "Date")[,3]

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Estimate the partitioning of rootmass  
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

rootmass.harvest.ratio = data.frame(matrix(ncol = 9, nrow = 252))
names(rootmass.harvest.ratio) = c("Date","FRratio","IRratio","CRratio","BRratio","FRratio_SE","IRratio_SE","CRratio_SE","BRratio_SE")
rootmass.harvest.ratio$Date = as.Date(as.Date("2013-09-17"):as.Date("2014-05-26"))

rootmass.harvest.ratio[nrow(rootmass.harvest.ratio),] <- summaryBy(FRratio+IRratio+CRratio+BRratio ~ Date, data=rootmass.harvest, FUN=c(mean,standard.error), keep.names=TRUE , na.rm=TRUE) # Sum of all same day Rd
# names(rootmass.harvest.ratio)[2:ncol(rootmass.harvest.ratio)] = c("FRratio","IRratio","CRratio","BRratio","FRratio_SE","IRratio_SE","CRratio_SE","BRratio_SE")
# rootmass.harvest.ratio.melt = melt(rootmass.harvest.ratio[,c("Date","FRratio","IRratio","CRratio","BRratio")], id.vars = "Date")
# rootmass.harvest.ratio.melt$biomass = as.factor("Belowground")
# rootmass.harvest.ratio.melt$SE = melt(rootmass.harvest.ratio[,c("Date","FRratio_SE","IRratio_SE","CRratio_SE","BRratio_SE")], id.vars = "Date")[,3]

#----------------------------------------------------------------------------------------------------------------
# Calculate the ratios of harvest rootmass from Court's Pot experiment
rootmass.pot.exp = read.csv("raw_data/seedling_mass_pot_exp.csv") # unit = mg of DM
rootmass.pot.exp = subset(rootmass.pot.exp, volume == 1000)
rootmass.pot.exp$rootmass = rootmass.pot.exp$coarseroot + rootmass.pot.exp$fineroot
rootmass.pot.exp$IRratio = rootmass.pot.exp$coarseroot / rootmass.pot.exp$rootmass
rootmass.pot.exp$FRratio = rootmass.pot.exp$fineroot / rootmass.pot.exp$rootmass

rootmass.harvest.ratio[1,c("CRratio","BRratio", "CRratio_SE","BRratio_SE")] = 0
rootmass.harvest.ratio[1,c("FRratio","IRratio", "FRratio_SE","IRratio_SE")] = summaryBy(FRratio+IRratio ~ volume, data=rootmass.pot.exp, FUN=c(mean,standard.error))[2:5]
for (i in 2:ncol(rootmass.harvest.ratio)) {
  rootmass.harvest.ratio[1:nrow(rootmass.harvest.ratio),i] = seq(rootmass.harvest.ratio[1,i], rootmass.harvest.ratio[nrow(rootmass.harvest.ratio),i], length.out = nrow(rootmass.harvest.ratio))
} 

rootmass.harvest.ratio.melt = melt(rootmass.harvest.ratio[,c("Date","FRratio","IRratio","CRratio","BRratio")], id.vars = "Date")
rootmass.harvest.ratio.melt$biomass = as.factor("Root")
rootmass.harvest.ratio.melt$SE = melt(rootmass.harvest.ratio[,c("Date","FRratio_SE","IRratio_SE","CRratio_SE","BRratio_SE")], id.vars = "Date")[,3]

treeMass.ratio = merge(woodmass.ratio, rootmass.harvest.ratio, by="Date", all=TRUE)
# treeMass.ratio = subset(treeMass.ratio, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-02-12"))
treeMass.ratio = subset(treeMass.ratio, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-05-26"))

treeMass.ratio.melt = rbind(woodmass.ratio.melt, rootmass.harvest.ratio.melt)

# Merge woodmass and rootmass partitioning data with daily GPP, Ra, LA, mass pool data
data.all = merge(data.GPP.Ra.LA.mass, treeMass.ratio, by="Date", all=TRUE)

# font.size = 12
pd <- position_dodge(0) # move the overlapped errorbars horizontally
p4 = ggplot(treeMass.ratio.melt, aes(x=Date, y=value, group = interaction(variable,biomass), colour=variable, shape=as.factor(biomass))) + 
  # geom_point(position=pd, size=1) +
  # geom_errorbar(position=pd, aes(ymin=value-SE, ymax=value+SE), colour="grey", width=1) +
  geom_ribbon(data = treeMass.ratio.melt, aes(ymin=value-SE, ymax=value+SE), linetype=2, alpha=0.1,size=0.1) +
  geom_line(position=pd,data = treeMass.ratio.melt, aes(x = Date, y = value, group = interaction(variable,biomass), colour=variable, linetype=biomass)) + 
  ylab("Biomass (wood & root) ratio") + xlab("") +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date), max(data.biomass$Date))) +
  labs(colour="Biomass types", linetype="Biomass") +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size=font.size-1)) +
  theme(legend.position = c(0.7,0.67), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("output/4.Woodmass_Rootmass_ratio.png", units="px", width=2000, height=1000, res=180)
print (p4)
dev.off()

print (p4)

