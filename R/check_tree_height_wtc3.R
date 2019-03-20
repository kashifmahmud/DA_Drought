# Read data form WTC3 experiment
height.dia.raw <- read.csv("raw_data/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV") # units: Height = cm, dia = mm
height.dia.raw$DateTime = as.Date(height.dia.raw$DateTime)
flux <- read.csv("raw_data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
chambers = unique(flux$chamber)

height.dia = data.frame(chamber = rep(chambers, each = 1),
                        Date = rep(as.Date("2014-05-27"),length(chambers)),
                        T_treatment = character(length(chambers)),
                        # W_treatment = character(length(chambers) * length(unique(height.dia.raw$DateTime))),
                        diameter = numeric(length(chambers)),
                        height = numeric(length(chambers)), stringsAsFactors=FALSE)
keeps <- c("chamber", "DateTime", "T_treatment", "Water_treatment", "Plant_height")
height.dia.sub = height.dia.raw[ , keeps, drop = FALSE]
names(height.dia.sub)[2] = "Date"
height.dia.sub = subset(height.dia.sub, Water_treatment %in% as.factor("control"))
height.dia.sub = subset(height.dia.sub, T_treatment %in% as.factor(c("ambient","elevated")))
height.dia.sub = subset(height.dia.sub, chamber %in% as.factor(c("C02","C05","C07","C09","C10","C12")))

plots = list()
pd <- position_dodge(0)
plots[[1]] = ggplot(height.dia.sub, aes(x=Date, y=Plant_height, group = interaction(chamber,T_treatment), colour=as.factor(chamber), shape=as.factor(T_treatment))) + 
  geom_point(position=pd) +
  # geom_errorbar(position=pd, aes(ymin=WM-WM_SE, ymax=WM+WM_SE), colour="grey", width=1) +
  # geom_ribbon(data = data.biomass, aes(ymin=WM-WM_SE, ymax=WM+WM_SE), linetype=2, alpha=0.1,size=0.1) +
  # geom_line(data = data.biomass, aes(x = Date, y = WM, group = interaction(T_treatment,chamber_type), colour=T_treatment, linetype=chamber_type)) +
  geom_hline(yintercept = 900) + ylab("Plant Height (cm)") +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(height.dia.sub$Date)-2, max(height.dia.sub$Date)+2)) +
  labs(colour="Chambers",shape="Treatment") +
  # scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  annotate("text", x = mean(height.dia.sub$Date), y = 915, size = font.size-7, label = paste("Chamber Height 9m line")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.3,0.75), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(file = "output/Plant_height.pdf")
do.call(grid.arrange,  plots)
dev.off()
