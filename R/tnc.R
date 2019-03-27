#-------------------------------------------------------------------------------------
# Read TNC data
tnc <- read.csv("raw_data/GHS30_Eglob-TxCxW_carbohydrates_20110117-20110321_L1.csv")
tnc$Date <- parse_date_time(tnc$Date, orders = "dmy")
tnc$Date <- as.Date(tnc$Date)

tnc <- subset(tnc, Temp == "Amb" & CO2 == 400) # only consider the ambient drought treatment
tnc = subset(tnc, Date <= as.Date("2011-02-27"))

keeps <- c("Date", "Water", "Organ", "StarchW", "SolSugW")
tnc = tnc[ , keeps, drop = FALSE]
tnc$tnc = tnc$StarchW + tnc$SolSugW
tnc$tnc.conc = tnc$tnc / 1000 # unit conversion: (mg of tnc per g of dry weight material) to (g of tnc per g of dry weight material)

# unit conversion: 1 g tnc has 0.4 gC (12/30) and 1 g plant has 0.48 gC
tnc$tnc.conc = tnc$tnc.conc * 0.4 / c1 # gC in tnc per gC of dry biomass
# tnc$tnc.conc_SE = tnc$tnc.conc_SE * 0.4 / c1 # gC in tnc per gC of dry biomass

#---------------------------------------------------------------------------------------------
# Plot allometry over time and treatment
plots = list()
font.size = 12
pd <- position_dodge(0.75)
plots[[1]] = ggplot(subset(tnc, Organ %in% as.factor("Leaf")), aes(x=Date, y=tnc.conc, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  # geom_errorbar(position=pd, aes(ymin=Height-Height_se, ymax=Height+Height_se), colour="grey", width=3) +
  geom_line(position=pd, aes(x = Date, y = tnc.conc, group = Water, colour=Water)) + 
  ylab(expression(Leaf ~TNC ~ (g ~C ~g ~C^{-1} ~ leafmass))) +
  # scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Treatment") +
  scale_color_manual(labels = c("Rewatered drought","Sustained drought","Well watered"), values = c("orange","red","green")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.2,0.85), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(colour=FALSE)

plots[[2]] = ggplot(subset(tnc, Organ %in% as.factor("Stem")), aes(x=Date, y=tnc.conc, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  # geom_errorbar(position=pd, aes(ymin=Height-Height_se, ymax=Height+Height_se), colour="grey", width=3) +
  geom_line(position=pd, aes(x = Date, y = tnc.conc, group = Water, colour=Water)) + 
  ylab(expression(Wood ~TNC ~ (g ~C ~g ~C^{-1} ~ woodmass))) +
  # scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Treatment") +
  scale_color_manual(labels = c("Rewatered drought","Sustained drought","Well watered"), values = c("orange","red","green")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.2,0.85), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(colour=FALSE)

plots[[3]] = ggplot(subset(tnc, Organ %in% as.factor("Root")), aes(x=Date, y=tnc.conc, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  # geom_errorbar(position=pd, aes(ymin=Height-Height_se, ymax=Height+Height_se), colour="grey", width=3) +
  geom_line(position=pd, aes(x = Date, y = tnc.conc, group = Water, colour=Water)) + 
  ylab(expression(Root ~TNC ~ (g ~C ~g ~C^{-1} ~ rootmass))) +
  # scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Treatment") +
  scale_color_manual(labels = c("Rewatered drought","Sustained drought","Well watered"), values = c("orange","red","green")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.2,0.85), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# do.call(grid.arrange,  plots)
png("output/6.tnc.png", units="px", width=700, height=1000, res=100)
grid.arrange(grobs=plots, ncol=1)
dev.off()

#-------------------------------------------------------------------------------------


