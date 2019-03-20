# Merge and plot all biomass data

# data.biomass = merge(rootmass, treeMass.sum, by = c("Date", "T_treatment", "chamber_type"), all=TRUE)
data.biomass = merge(rootmass, treeMass.sum, by = c("Date", "T_treatment"), all=TRUE)
# data.biomass = subset(data.biomass, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-02-12"))

plots = list()
pd <- position_dodge(3)
plots[[1]] = ggplot(data.biomass, aes(x=Date, y=WM, group = T_treatment, colour=T_treatment)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=WM-WM_SE, ymax=WM+WM_SE), colour="grey", width=1) +
  # geom_ribbon(data = data.biomass, aes(ymin=WM-WM_SE, ymax=WM+WM_SE), linetype=2, alpha=0.1,size=0.1) +
  # geom_line(data = data.biomass, aes(x = Date, y = WM, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) + 
  ylab(expression("Wood Mass"~"(g C "*plant^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Temperature") +
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

plots[[2]] = ggplot(data.biomass, aes(x=Date, y=RM, group = T_treatment, colour=T_treatment)) + 
  geom_point(position=pd,size=2) +
  # geom_ribbon(position=pd,data = data.biomass, aes(ymin=RM-RM_SE, ymax=RM+RM_SE), linetype=2, alpha=0.1,size=0.1) +
  geom_errorbar(position=pd,aes(ymin=RM-RM_SE, ymax=RM+RM_SE), colour="grey", width=1) +
  # geom_line(data = data.biomass, aes(x = Date, y = RM, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) +
  ylab(expression("Root Mass"~"(g C "*plant^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Temperature") +
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

plots[[3]] = ggplot(data.biomass, aes(x=Date, y=LM, group = T_treatment, colour=T_treatment)) + 
  geom_point(position=pd) +
  # geom_ribbon(data = data.biomass, aes(ymin=LM-LM_SE, ymax=LM+LM_SE), linetype=2, alpha=0.1,size=0.1) +
  geom_errorbar(position=pd,aes(ymin=LM-LM_SE, ymax=LM+LM_SE), colour="grey", width=1) +
  # geom_line(data = data.biomass, aes(x = Date, y = LM, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) + 
  ylab(expression("Foliage Mass"~"(g C "*plant^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Temperature") +
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

pdf(file = "output/2.Biomass.pdf")
do.call(grid.arrange,  plots)
dev.off()

do.call(grid.arrange,  plots)
