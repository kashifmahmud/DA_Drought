# Check Rabove from data
data.all = read.csv("processed_data/data_all.csv") 
keeps <- c("Date", "Treatment", "Ra", "LM", "WM", "Rd.foliage.mean", "Rd.stem.mean", "SMratio", "Rd.branch.mean", "BMratio")
data.all = data.all[ , keeps, drop = FALSE]

for (j in 1:length(treat.group)) {
  data = subset(data.all,(Treatment %in% treat.group[j]))
  data$Ra = cumsum(data$Ra)
  
  data.set = data
  data.set[,c("LM","WM")] = na.spline(data.set[,c("LM","WM")])
  for (i in 1:nrow(data.set)) {
    data.set$Resp.leaf[i] = data.set$Rd.foliage.mean[i]*data.set$LM[i]
    data.set$Resp.stem[i] = data.set$Rd.stem.mean[i]*data.set$WM[i]*data.set$SMratio[i] 
    data.set$Resp.branch[i] = data.set$Rd.branch.mean[i]*data.set$WM[i]*data.set$BMratio[i]
  }
  data.set$Resp.growth[1] = 0
  for (i in 2:nrow(data.set)) {
    data.set$Resp.growth[i] = 0.3*(data.set$LM[i]-data.set$LM[i-1] + data.set$WM[i]-data.set$WM[i-1])
  }
  data.set$Resp.leaf = cumsum(data.set$Resp.leaf); data.set$Resp.stem = cumsum(data.set$Resp.stem); 
  data.set$Resp.branch = cumsum(data.set$Resp.branch); data.set$Resp.growth = cumsum(data.set$Resp.growth); 
  data.set$Rabove = data.set$Resp.leaf + data.set$Resp.stem + data.set$Resp.branch + data.set$Resp.growth
  
  if (j == 1) {
    data.set.final = data.set
  }
  if (j > 1) {
    data.set.final = rbind(data.set,data.set.final)
  }
}
keeps = c("Date","Treatment","Ra","Resp.leaf","Resp.stem","Resp.branch","Resp.growth")
data.set.final = data.set.final[ , keeps, drop = FALSE]
data.set.final$Resp.leaf.stem = data.set.final$Resp.leaf + data.set.final$Resp.stem
data.set.final$Resp.leaf.stem.branch = data.set.final$Resp.leaf.stem + data.set.final$Resp.branch
data.set.final$Resp.leaf.stem.branch.growth = data.set.final$Resp.leaf.stem.branch + data.set.final$Resp.growth
# data.set.final[data.set.final <= 0] <- 0

# data.set.final.melt = melt(data.set.final[,c("Date","Treatment","Ra","Rabove")], id.vars=c("Treatment","Date"))
# data.set.final.melt$Date = as.Date(data.set.final.melt$Date)

data.set.final.melt = melt(data.set.final[,c("Date","Treatment","Ra","Resp.leaf","Resp.leaf.stem","Resp.leaf.stem.branch","Resp.leaf.stem.branch.growth")], id.vars=c("Treatment","Date"))
data.set.final.melt$Date = as.Date(data.set.final.melt$Date)

plot = list() 
font.size = 12
pd <- position_dodge(0) # move the overlapped errorbars horizontally
# cbPalette = c("firebrick", "green3")
cbPalette = c("firebrick", "orange", "green3", "yellow3", "#0072B2", "#D55E00", "gray")
title = as.character(c("Room#1","Room#2","Room#3","Room#4","Room#5","Room#6"))
for (v in 1:length(treat.group)) {
  data.set.final.melt.set = subset(data.set.final.melt,(Treatment %in% treat.group[v]))
  plot[[v]] = ggplot() + 
    geom_line(position=pd,data = data.set.final.melt.set, aes(x = Date, y = value, group = variable, colour=variable)) +
    ylab(expression(R[above]~"Cumulative (g C "*plant^"-1"*")")) + xlab("Date") +
    labs(colour="Measurements") + 
    # scale_y_continuous(trans = 'log10') +
    scale_colour_manual(breaks=c("Ra","Resp.leaf","Resp.leaf.stem","Resp.leaf.stem.branch","Resp.leaf.stem.branch.growth"), labels=c(expression(~R[a(from~flux)]),"Leaf resp",
                "Leaf + Stem resp","Leaf + Stem + branch resp","Leaf + Stem + Branch + Growth resp"),
                        values=cbPalette) +
    # scale_colour_manual(breaks=c("Ra","Rabove"), labels=c(expression(~R[a(from~flux)]),expression(~R[a(from~components)])),
    #                     values=cbPalette) +
    # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g C plant-1)")) +
    theme_bw() +
    annotate("text", x = mean(data.set.final.melt.set$Date), y = max(data.set.final.melt.set$value)-2, size = font.size-9, label = paste("Considering 30% Growth Respiration")) +
    annotate("text", x = max(data.set.final.melt.set$Date)-5, y = min(data.set.final.melt.set$value)+5, size = font.size-8, label = paste(as.character(treat.group[v]))) +
    # theme(plot.title = element_text(size = 20, face = "bold")) +
    theme(legend.title = element_text(colour="black", size=font.size-2)) +
    theme(legend.text = element_text(colour="black", size = font.size-3)) +
    theme(legend.key.height=unit(0.6,"line")) +
    theme(legend.position = c(0.25,0.75),legend.direction = "vertical", legend.text.align = 0) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    # theme(plot.title = element_text(hjust = 0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
}

png("output/Figure_Rabove_balance_v2.png", units="px", width=1600, height=1300, res=220)
print (do.call(grid.arrange,  plot))
dev.off()


