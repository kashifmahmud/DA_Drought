
#-------------------------------------------------------------------------------------
# Read allometry data
allometry <- read.csv("raw_data/GHS30_Eglob-TxCxW_growth_20110107-20110321_L1.csv")
allometry$Date <- as.Date(allometry$Date, format="%d/%m/%Y")

allometry <- subset(allometry, Temp == "Amb" & CO2 == 400) # only consider the ambient drought treatment

allometry$D = (allometry$Basaldia1 + allometry$Basaldia2) / 2 # Calculate mean diameter
allometry <- (subset(allometry, select = c("Date", "Water", "Height", "D")))

allometry.df = summaryBy(Height+D ~ Date+Water, data=allometry, FUN=c(mean,standard.error))
names(allometry.df) = c("Date", "Water", "Height", "D", "Height_se", "D_se")

allometry.df = subset(allometry.df, Date <= as.Date("2011-02-27"))

#---------------------------------------------------------------------------------------------
# Plot allometry over time and treatment
plots = list()
font.size = 12
pd <- position_dodge(0.75)
plots[[1]] = ggplot(allometry.df, aes(x=Date, y=Height, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=Height-Height_se, ymax=Height+Height_se), colour="grey", width=3) +
  geom_line(position=pd, data = allometry.df, aes(x = Date, y = Height, group = Water, colour=Water)) +
  ylab(expression(H ~ (cm))) +
  # scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Treatment") +
  scale_color_manual(labels = c("Rewatered drought","Sustained drought","Well watered"), values = c("orange","red","green")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.2,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

plots[[2]] = ggplot(allometry.df, aes(x=Date, y=D, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=D-D_se, ymax=D+D_se), colour="grey", width=3) +
  geom_line(position=pd, data = allometry.df, aes(x = Date, y = D, group = Water, colour=Water)) +
  ylab(expression(D ~ (mm))) +
  # scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Treatment") +
  scale_color_manual(labels = c("Rewatered drought","Sustained drought","Well watered"), values = c("orange","red","green")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.2,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

do.call(grid.arrange,  plots)
png("output/4.Allometry.png", units="px", width=700, height=750, res=100)
grid.arrange(grobs=plots, ncol=1)
dev.off()

#-------------------------------------------------------------------------------------