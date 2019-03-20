# Merge leaf tnc with all biomass data and plot the leaf tnc storage data

data.biomass = merge(data.biomass, tnc.final[,c("Date", "T_treatment", "TNC_tot", "TNC_tot_SE", "TNC_leaf", "TNC_leaf_SE", "TNC_wood", "TNC_wood_SE", "TNC_root", "TNC_root_SE")], all=TRUE)

data.GPP.Ra.LA.mass = merge(data, data.biomass, by=c("Date","T_treatment"), all=TRUE)

pd <- position_dodge(2)
p3 = ggplot(data.biomass[complete.cases(data.biomass$TNC_leaf),], aes(x=Date, y=TNC_leaf, group = T_treatment, colour=T_treatment)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=TNC_leaf-TNC_leaf_SE, ymax=TNC_leaf+TNC_leaf_SE), colour="grey", width=1) +
  # geom_line(position=pd, data = data.biomass[complete.cases(data.biomass$TNC),], aes(x = Date, y = TNC, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) + 
  ylab(expression(C["n,f"]~"(g C "*plant^"-1"*")")) +
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

pdf(file = "output/4.leaf_tnc.pdf")
print (p3)
dev.off()

print (p3)
