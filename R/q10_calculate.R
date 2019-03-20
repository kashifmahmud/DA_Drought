#- process leaf-scale R vs. T data to find q10 value for WTC3
rvt <- read.csv("raw_data/WTC_TEMP_CM_GX-RdarkVsT_20140207-20140423_L1.csv")
rvt$Date <- as.Date(rvt$Date)
rvt.sub <- subset(rvt,Date==as.Date("2014-02-07")) #- just pull out the data measured prior to the drought

rvt.45 <- subset(rvt.sub, Tleaf > 18 & Tleaf<=40)
rvt.45$lnRmass <- log(rvt.45$Rmass)
rvt.45$varT <- (rvt.45$Tleaf-25)/10

rvt.45$Tleaf_bin <- cut(rvt.45$Tleaf,breaks=seq(from=18,to=40,length=25))
rvt.45$Tleaf_bin_mid <- sapply(strsplit(gsub("^\\W|\\W$", "", rvt.45$Tleaf_bin), ","), function(x)sum(as.numeric(x))/2) 
rvt.treat.bin <- summaryBy(Rmass+lnRmass+varT~date+T_treatment+Water_treatment+Tleaf_bin_mid,data=rvt.45,keep.names=T,FUN=mean)

#- fit R vs T equation  
lmRleaf <- lm(lnRmass~varT+T_treatment+Water_treatment,data=rvt.treat.bin)
summary(lmRleaf)
coef(lmRleaf) 
# q10_amb_wet <- unname(exp(coef(lmRleaf)[2]))
# q10_ele_wet <- unname(exp(coef(lmRleaf)[2]+coef(lmRleaf)[3]))
# q10_amb_dry <- unname(exp(coef(lmRleaf)[2]+coef(lmRleaf)[4]))
# q10_ele_dry <- unname(exp(coef(lmRleaf)[2]+coef(lmRleaf)[3]+coef(lmRleaf)[4]))

# No significant difference for drought/control treatments, so only consider temperature treatment effect
lmRleaf_amb <- lm(lnRmass~varT,data=subset(rvt.treat.bin,T_treatment %in% as.factor("ambient")))
q10_amb <- unname(exp(coef(lmRleaf_amb)[2]))
rd25_amb = unname(exp(coef(lmRleaf_amb)[1]))
lmRleaf_ele <- lm(lnRmass~varT,data=subset(rvt.treat.bin,T_treatment %in% as.factor("elevated")))
q10_ele <- unname(exp(coef(lmRleaf_ele)[2]))
rd25_ele = unname(exp(coef(lmRleaf_ele)[1]))

# get model predictions of R across all temperatures
xvals_Rleaf <- seq(18,40, length=101)
# predRleaf_amb <- exp(log(rd25_amb) + (log(q10_amb) * (xvals_Rleaf-25)/10))
# predRleaf_ele <- exp(log(rd25_ele) + (log(q10_ele) * (xvals_Rleaf-25)/10))

# No significant differences between treatments (according to John's paper)
lmRleaf <- lm(lnRmass~varT,data=rvt.treat.bin)
q10 <- unname(exp(coef(lmRleaf)[2])) ########## this is final q10 value for all treatments
predRleaf_amb <- exp(log(rd25_amb) + (log(q10) * (xvals_Rleaf-25)/10))
predRleaf_ele <- exp(log(rd25_ele) + (log(q10) * (xvals_Rleaf-25)/10))

#----------------------------------------------------------------------------------------------------------------
#- plot Rleaf vs. Tair
rvt.treat <- summaryBy(Rmass~Date+T_treatment+Tleaf_bin_mid,data=rvt.treat.bin,keep.names=F,FUN=c(mean,standard.error))
png("output/5.Rdark_vs_T.png", units="px", width=1200, height=1000, res=300)
par(mfrow=c(1,1), mar = c(4,4.5,0.2,0.2))
xlims=c(18,40)

# #- plot MASS BASED leaf R over T for the first date of high resolution T-response curves
plotBy(Rmass.mean~Tleaf_bin_mid|T_treatment,data=rvt.treat,
       xlab=expression(Temperature~(degree*C)),ylab=expression(R[leaf]~(mu*mol~CO[2]~m^-2~s^-1)),col=c("black","red"),pch=1,
       xlim=xlims,type="p",lwd=3,cex=0.3,cex.lab=1,legend=F,
       panel.first=adderrorbars(x=rvt.treat$Tleaf_bin_mid,y=rvt.treat$Rmass.mean,SE=rvt.treat$Rmass.standard.error,direction="updown"))
# magaxis(side=c(1,2,3,4),labels=c(1,1,0,1))
# title(xlab=expression(Temperature~(degree*C)),xpd=NA)
lines(x=xvals_Rleaf,y=predRleaf_amb,col="black",lwd=2)
lines(x=xvals_Rleaf,y=predRleaf_ele,col="red",lwd=2)
legend("topleft",legend=unique(rvt.treat$T_treatment),col=c('black','red'),lty=1,bty="n",cex=0.8,pt.cex=1.2)

# dev.copy2pdf(file="output/5.Rdark_vs_T.pdf")
dev.off()

img <- readPNG("output/5.Rdark_vs_T.png")
# h<-dim(img)[1]
# w<-dim(img)[2]
# 
# # #open new file for saving the image in "output" folder
# # png("output/Figure_1_CBM_wtc3.png", width=w, height=h)
# par(mar=c(0,0,0,0), xpd=NA, mgp=c(0,0,0), oma=c(0,0,0,0), ann=F)
# plot.new()
# plot.window(0:1, 0:1)
# 
# #fill plot with image
# usr<-par("usr")
# rasterImage(img, usr[1], usr[3], usr[2], usr[4])
grid.raster(img)

# plot.q10()

