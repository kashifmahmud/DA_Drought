# Script to read and plot stem height and diameter for various treatment cases 

# Read data from WTC3 experiment
height.dia.raw <- read.csv("raw_data/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV")
height.dia.raw$DateTime = as.Date(height.dia.raw$DateTime)
flux <- read.csv("raw_data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
chambers = unique(flux$chamber)

height.dia = data.frame(chamber = rep(chambers, each = length(unique(height.dia.raw$DateTime))),
                        Date = rep(unique(height.dia.raw$DateTime),length(chambers)),
                        T_treatment = character(length(chambers) * length(unique(height.dia.raw$DateTime))),
                        W_treatment = character(length(chambers) * length(unique(height.dia.raw$DateTime))),
                        diameter = numeric(length(chambers) * length(unique(height.dia.raw$DateTime))),
                        height = numeric(length(chambers) * length(unique(height.dia.raw$DateTime))), stringsAsFactors=FALSE)
height.dia = subset(height.dia, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-05-27"))

for(i in 1:length(chambers)) {
  height.dia.sub = subset(height.dia.raw, chamber %in% as.factor(chambers[i]))
  # height.dia.sub = subset(height.dia, chamber %in% as.factor(chambers[i]) & Water_treatment %in% as.factor("control"))
  if (i==11) {  # Remove the reference measurements made on Chamber 11
    height.dia.sub.1 = subset(height.dia.sub, Stem_number %in% 1)
    height.dia.sub.1.2 = subset(height.dia.sub, Stem_number %in% 1.2)
    height.dia.sub.1[height.dia.sub.1$DateTime >= as.Date("2013-12-24"),"Plant_height"] = height.dia.sub.1.2[,"Plant_height"]
    height.dia.sub = height.dia.sub.1
  }
  keeps <- c("chamber", "DateTime", "T_treatment", "Water_treatment", "Plant_height", "X15", "X65")
  height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
  
  D.15 <- lm(X15 ~ X65, data=height.dia.sub)
  # visreg(D.15, "X65", overlay=TRUE)
  # summary(D.15)
  eq.D = function(x){coefficients(D.15)[1] + coefficients(D.15)[2] * x }
  
  height.dia.sub$X15 = eq.D(height.dia.sub$X65)
  height.dia.sub = subset(height.dia.sub, DateTime >= as.Date("2013-09-14") & DateTime <= as.Date("2014-05-27"))
  height.dia.sub = height.dia.sub[!is.na(height.dia.sub$X65),]
  # keeps <- c("chamber", "DateTime", "T_treatment", "Water_treatment", "X15", "Plant_height")
  keeps <- c("chamber", "DateTime", "T_treatment", "Water_treatment", "X65", "Plant_height")
  height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
  names(height.dia.sub) <- c("chamber","Date","T_treatment","W_treatment","diameter","height")
  height.dia.sub$T_treatment = as.character(height.dia.sub$T_treatment)
  height.dia.sub$W_treatment = as.character(height.dia.sub$W_treatment)
  
  height.dia[(1+(i-1)*length(unique(height.dia$Date))) : (i*length(unique(height.dia$Date))), 
             c("T_treatment","W_treatment","diameter","height")] = height.dia.sub[,c("T_treatment","W_treatment","diameter","height")]
}
height.dia$T_treatment = as.factor(height.dia$T_treatment)
height.dia$W_treatment = as.factor(height.dia$W_treatment)

# Average the ambient and elevated temperature treatments regardless of drought/watered treatment
# n=6 for whole period
height.dia.idn <- summaryBy(height+diameter ~ Date+T_treatment, data=height.dia, FUN=c(mean,standard.error,length))
height.dia.idn$W_treatment = as.factor("all")
height.dia.idn$case = as.factor("1")

# Average the ambient and elevated temperature treatments considering the drought/watered treatment from February 2014
# 1st haf: n=6, 2nd half: n=3
height.dia.final <- summaryBy(height+diameter ~ Date+T_treatment+W_treatment, data=height.dia, FUN=c(mean,standard.error,length))
height.dia.final$case = as.factor("2")
height.dia.final = rbind(height.dia.idn, height.dia.final)

# Average the ambient and elevated temperature treatments considering the drought/watered treatment seperated from the start of the experiment
# n=3 for whole period
drought.chamb = unique(height.dia$chamber[ height.dia$W_treatment %in% as.factor("drydown")])
height.dia.sub = height.dia
height.dia.sub$chamber_type = as.factor( ifelse(height.dia.sub$chamber %in% drought.chamb, "drought_treat", "watered_treat") )
height.dia.idn <- summaryBy(height+diameter ~ Date+T_treatment+chamber_type, data=height.dia.sub, FUN=c(mean,standard.error,length))
names(height.dia.idn)[3] = "W_treatment"
height.dia.idn$case = as.factor("3")
height.dia.final = rbind(height.dia.final, height.dia.idn)
names(height.dia.final)[3:6] = c("height", "diameter", "height_SE", "diameter_SE")


# Plot both height and diameter for various scenarios
melted.height.dia.final = melt(height.dia.final[,c("diameter","height","Date","T_treatment","W_treatment","case")], id.vars=c("Date","T_treatment","W_treatment","case"))
melted.height.dia.error = melt(height.dia.final[,c("diameter_SE","height_SE","Date","T_treatment","W_treatment","case")], id.vars=c("Date","T_treatment","W_treatment","case"))
# melted.height.dia.final = melt(height.dia.final[,c("diameter","height","Date","T_treatment","case")], id.vars=c("Date","T_treatment","case"))
# melted.height.dia.error = melt(height.dia.final[,c("diameter_SE","height_SE","Date","T_treatment","case")], id.vars=c("Date","T_treatment","case"))
melted.height.dia.final$SE = melted.height.dia.error$value

# add.points = subset(melted.height.dia.final,Date %in% as.Date("2014-02-04") & W_treatment %in% as.factor("control"))
add.points = subset(melted.height.dia.final, Date <= as.Date("2014-02-04") & W_treatment %in% as.factor("control"))
add.points$W_treatment = as.factor("drydown")
melted.height.dia.final = rbind(melted.height.dia.final, add.points)

font.size = 12
pd <- position_dodge(0) # move the overlapped errorbars horizontally
cbPalette = c("black", "green3", "red", "magenta", "blue")
melted.height.dia.final.trait = subset(melted.height.dia.final,variable %in% as.factor("diameter") & T_treatment %in% as.factor("elevated"))

p1 = ggplot(melted.height.dia.final.trait, aes(x=Date, y=value, group = W_treatment, colour=W_treatment)) + 
  geom_point(position=pd) +
  geom_line(position=pd,data = melted.height.dia.final.trait, aes(x = Date, y = value, group = W_treatment, colour=W_treatment)) +
  ylab("Diameter (mm)") +
  xlab("Month") +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
  labs(colour="Treatments") +
  scale_colour_manual(breaks=c("all","control","watered_treat","drydown","drought_treat"), labels=c("All Warm chambers (n=6)","Wet (n=6:predrought to n=3:postdrought)","Wet (n=3)",
                                                                                                    "Dry (n=6:predrought to n=3:postdrought)","Dry (n=3)"),values=cbPalette) +
  ggtitle(" Elevated temperature treatments") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size-2)) +
  theme(legend.position = c(0.2,0.75)) + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size, family="Arial")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
print(p1)

png("output/1.Diameter_elevated_predrought_difference.png", units="px", width=2200, height=1600, res=220)
print(p1)
dev.off()
# dev.copy2pdf(file="output/6.Diameter_elevated_predrought_difference.pdf")

i = 0
font.size = 12
plots = list() 
meas = as.factor(c("diameter","height"))
error = as.factor(c("diameter_SE","height_SE"))
temp = as.factor(c("ambient","elevated"))
pd <- position_dodge(0) # move the overlapped errorbars horizontally
for (p in 1:length(meas)) {
  for (q in 1:length(temp)) {
    melted.height.dia.final.trait = subset(melted.height.dia.final,variable %in% meas[p] & T_treatment %in% temp[q])
    # melted.height.dia.final.trait = subset(melted.height.dia.final,variable %in% meas[p] & T_treatment %in% as.factor("elevated") & case %in% as.factor("1"))
    # melted.height.dia.final.trait = subset(melted.height.dia.final,variable %in% meas[p] & T_treatment %in% as.factor("ambient") & case %in% as.factor("3"))
    
    i = 2*(p-1) + q
    plots[[i]] = ggplot(melted.height.dia.final.trait, aes(x=Date, y=value, group = W_treatment, colour=W_treatment)) + 
      geom_point(position=pd) +
      geom_errorbar(position=pd,aes(ymin=value-SE, ymax=value+SE), colour="grey", width=2) +
      geom_line(position=pd,data = melted.height.dia.final.trait, aes(x = Date, y = value, group = W_treatment, colour=W_treatment)) +
      ylab(paste(as.character(meas[p]),"(mm)")) +
      xlab("Month") +
      scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
      labs(colour="Treatments") +
      # scale_colour_manual(breaks=c("all","control","watered_treat","drydown","drought_treat"), labels=c("All Warm chambers (n=6)","Wet (n=6:predrought to n=3:postdrought)","Wet (n=3)",
      #                                                "Dry (n=6:predrought to n=3:postdrought)","Dry (n=3)"),values=cbPalette) +
      ggtitle(paste(as.character(temp[q]), "temperature")) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.title = element_text(colour="black", size=font.size)) +
      theme(legend.text = element_text(colour="black", size = font.size-2)) +
      theme(legend.position = c(0.35,0.88)) + theme(legend.key.height=unit(0.8,"line")) +
      theme(legend.key = element_blank()) +
      theme(text = element_text(size=font.size, family="Arial")) +
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    if (p==2) {
      plots[[i]] = plots[[i]] + ylab(paste(as.character(meas[p]),"(cm)"))
    }
    if (q==1) {
      plots[[i]] = plots[[i]] + guides(colour=FALSE)
    }
  }
}

for (r in 1:length(meas)) {
  melted.height.dia.final.trait = subset(melted.height.dia.final,variable %in% meas[r])
  
  i = i + 1
  plots[[i]] = ggplot(melted.height.dia.final.trait, aes(x=Date, y=value, group = interaction(T_treatment,W_treatment), colour=T_treatment, shape=W_treatment)) + 
    geom_point(position=pd) +
    # geom_errorbar(position=pd,aes(ymin=value-SE, ymax=value+SE), colour="grey", width=2) +
    geom_line(position=pd,data = melted.height.dia.final.trait, aes(x = Date, y = value, group = interaction(T_treatment,W_treatment), colour=T_treatment, linetype=W_treatment)) +
    ylab(paste(as.character(meas[r]),"(mm)")) +
    xlab("Month") +
    scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
    labs(colour="Temperature", shape="Water", linetype="Water") +
    scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
    theme_bw() +
    theme(legend.title = element_text(colour="black", size=font.size)) +
    theme(legend.text = element_text(colour="black", size = font.size)) +
    theme(legend.position = c(0.9,0.3)) + theme(legend.key.height=unit(0.7,"line")) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size, family="Arial")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  if (r==2) {
    plots[[i]] = plots[[i]] + ylab(paste(as.character(meas[r]),"(cm)"))
  }
}

png("output/1.dia_over_time.png", units="px", width=2200, height=1600, res=220)
lay <- rbind(c(1,1),c(2,3))
grid.arrange(grobs = plots[c(5,1,2)], layout_matrix = lay)
# lay <- rbind(c(1,1),c(1,1),c(2,2),c(3,3))
# grid.arrange(grobs = plots[c(5,1,2)], layout_matrix = lay)
dev.off()
png("output/1.height_over_time.png", units="px", width=2200, height=1600, res=220)
grid.arrange(grobs = plots[c(6,3,4)], layout_matrix = lay)
dev.off()

