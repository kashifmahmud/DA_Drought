
# Read in R dark spot measurements
rdark <- read.csv("raw_data/GHS30_Eglob-TxCxW_GErdark_20110117-20110321_L1.csv")
rdark$Date <- as.Date(rdark$Date,format="%d/%m/%Y")

# Clean the Licor data
rdark = subset(rdark, PARi >= 0 & Ci >= 0 & Cond >= 0)
rdark <- (subset(rdark, select = c("Date", "Potnum", "Temp", "CO2", "Water", "Photo", "Cond", 
                                                      "Ci", "Tleaf", "Area", "Trmmol")))
rdark <- subset(rdark, Temp == "Amb" & CO2 == 400) # only consider the ambient drought treatment
# rdark <- subset(rdark, CO2 == 400) # only consider the ambient drought treatment
# boxplot(Photo~ Water, data=rdark)

# rdark_clean <- rdark[rdark$Photo >= -.5,] # ???
# boxplot(Photo~ Water, data=rdark_clean)

#---------------------------------------------------------------------------------------------
### need to calculate a q10 function to model Rdark through time
## use rdark at 25c for upper end, average the data by water treatment
q25_drake <- 1.86
q10_crous <- 1.95

# water treatment means
rdark_agg = summaryBy(Photo+Cond+Ci+Trmmol+Tleaf ~ Date+Water, data=rdark, FUN=mean)
names(rdark_agg) = c("Date", "Water", "Photo", "Cond", "Ci", "Trmmol", "Tleaf")
rdark_agg$Photo <- -1*rdark_agg$Photo

rdark_eq <- cbind(subset(rdark_agg, select = c("Date", "Water", "Photo", "Tleaf")), q25_drake)
rdark_eq <- cbind(rdark_eq, q10_crous)
mean(rdark_eq$Tleaf)
names(rdark_eq)[3] <- "rd18.5"

# q10 equation 
rdark_eq$rd25_euct <- with(rdark_eq, rd18.5*(q25_drake^((abs(Tleaf-25))/10)))
rdark_eq$rd25_eucs <- with(rdark_eq, rd18.5*(q10_crous^((abs(Tleaf-25))/10)))

# rdark_eq$rd25_euct <- with(rdark_eq, (q25_drake/rd18.4)^(10/(25-Tleaf)))
# rdark_eq$rd25_eucs <- with(rdark_eq, (q10_crous/rd18.4)^(10/(25-Tleaf)))

# Test whether there are treatment and temporal effect on Rd
rdark_Anova <- lm(rd18.5 ~ Water+Date, data = rdark_eq)
anova(rdark_Anova)
# Both treatment effect and temporal effect on Rd

#---------------------------------------------------------------------------------------------
# Plot Rd over time and treatment
font.size = 12
pd <- position_dodge(2)
p1 = ggplot(rdark_eq, aes(x=Date, y=rd18.5, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  # geom_errorbar(position=pd, aes(ymin=g1-g1_se, ymax=g1+g1_se), colour="grey", width=3) +
  geom_line(position=pd, data = rdark_eq, aes(x = Date, y = rd18.5, group = Water, colour=Water)) +
  ylab(expression(R["d,18.5"] ~ (mu ~ mol ~ m^{-2} ~ s^{-1}))) +
  # scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Treatment") +
  # scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.3,0.2), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

png(file = "output/1.Rd.png")
print (p1)
dev.off()

write.csv(rdark_eq[,c(1:3)], "processed_data/rdarkq10.csv", row.names=FALSE)

#---------------------------------------------------------------------------------------------

