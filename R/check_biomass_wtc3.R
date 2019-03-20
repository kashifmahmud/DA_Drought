# Read biomass
biomass.sub = subset(data.biomass, Date %in% as.Date(c("2013-09-17","2014-05-26")))

biomass.sub$biomass = biomass.sub$RM + biomass.sub$LM + biomass.sub$WM

biomass.sub = subset(biomass.sub, Date %in% as.Date("2014-05-26"))
biomass.sub[c(3:4),c("LM","RM","WM")] = biomass.sub[,c("LM","RM","WM")] / biomass.sub$biomass[c(1:2)]
biomass.sub[c(3:4),c("LM","RM","WM")]

#-------------------------------------------------------------------------------------
# Calculate foliage litterfall rate
# Script to process leaf litter data
litterfall = read.csv("raw_data/WTC_TEMP_CM_LEAFLITTER_20130913-20140528_L1.csv")
litterfall$startDate = as.Date(litterfall$startDate)
litterfall$collectionDate = as.Date(litterfall$collectionDate)
litterfall$Date <- (litterfall$startDate + ((litterfall$collectionDate - litterfall$startDate) / 2))
litterfall = subset(litterfall, Date >= as.Date("2013-09-14") & Date <= as.Date("2014-05-27"))
# litterfall = subset(litterfall, Date >= as.Date("2013-09-14") & Date <= as.Date("2014-02-12"))

# convert to data.table in place
litterfall = setDT(litterfall)
# dcast and do individual sums
litterfall.cast = dcast.data.table(litterfall, chamber ~ Date, value.var = 'litter', fun.aggregate = sum)
litterfall.cast <- melt(litterfall.cast, id.vars = "chamber")
litterfall.cast = merge(litterfall.cast, unique(treeMass[,c("chamber","T_treatment")]), all=TRUE)
# litterfall.cum.melt$chamber_type = as.factor( ifelse(litterfall.cum.melt$chamber %in% drought.chamb, "drought", "watered") )
names(litterfall.cast)[2:3] = c("Date","litter")
litterfall.cast$Date = as.Date(litterfall.cast$Date)
litterfall.cast = summaryBy(litter ~ Date+T_treatment, data=litterfall.cast, FUN=c(mean,standard.error))
names(litterfall.cast)[3:4] = c("litter","litter_SE")

# Foliage mass data
lm = data.biomass[,c("Date","T_treatment","LM","LM_SE")]
lm.litter = merge(lm,litterfall.cast, by=c("Date","T_treatment"), all = TRUE)
lm.litter[,c(3:4)] = na.spline(lm.litter[,c(3:4)])
lm.litter$litterrate = lm.litter$litter / lm.litter$LM / 14
lm.litter$litterrate_SE = (((lm.litter$litter_SE/lm.litter$litter)^2 + (lm.litter$LM_SE/lm.litter$LM)^2)^0.5) * lm.litter$litterrate

lm.litter = na.omit(lm.litter)

mean(lm.litter$litterrate[which(lm.litter$T_treatment %in% as.factor("ambient"))], na.rm = T)
mean(lm.litter$litterrate[which(lm.litter$T_treatment %in% as.factor("elevated"))], na.rm = T)

plots = list()
pd <- position_dodge(2)
plots[[1]] = ggplot(lm.litter, aes(x=Date, y=litterrate, group = as.factor(T_treatment), colour=as.factor(T_treatment))) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=litterrate-litterrate_SE, ymax=litterrate+litterrate_SE), colour="grey", width=2) +
  geom_smooth(method='lm', se = FALSE) +
  # geom_ribbon(data = data.biomass, aes(ymin=WM-WM_SE, ymax=WM+WM_SE), linetype=2, alpha=0.1,size=0.1) +
  # geom_line(data = data.biomass, aes(x = Date, y = WM, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) +
  # geom_hline(yintercept = 900) + ylab("Plant Height (cm)") +
  ylab("Leaf litter turnover rate (% of Leaf mass)") +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(lm.litter$Date)-2, max(lm.litter$Date)+2)) +
  labs(colour="Treatment") +
  scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  # annotate("text", x = mean(lm.litter$Date), y = 915, size = font.size-7, label = paste("Chamber Height 9m line")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.5,0.85), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  # theme(axis.text.x = element_text(angle=0, hjust = 0)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("output/Leaf_litter_turnover_rate.png", units="px", width=1200, height=800, res=150)
# pdf(file = "output/Litterfall_rates.pdf")
do.call(grid.arrange,  plots)
dev.off()

