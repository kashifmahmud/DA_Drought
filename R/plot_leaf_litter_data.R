# Merge leaf litter with all biomass data and plot the leaf litterfall data

data.biomass = merge(data.biomass, litterfall.cum.melt, all=TRUE)

pd <- position_dodge(3)
p2 = ggplot(data.biomass[complete.cases(data.biomass$litter),], aes(x=Date, y=litter, group = T_treatment, colour=T_treatment)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=litter-litter_SE, ymax=litter+litter_SE), colour="grey", width=1) +
  geom_line(position=pd, data = data.biomass[complete.cases(data.biomass$litter),], aes(x = Date, y = litter, group = T_treatment, colour=T_treatment)) +
  ylab(expression("Leaf Litter"~"(g C "*plant^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data.biomass$Date)-2, max(data.biomass$Date)+2)) +
  labs(colour="Temperature") +
  scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.25,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

pdf(file = "output/3.leaf_litter.pdf")
print (p2)
dev.off()

print (p2)
