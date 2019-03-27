
#-------------------------------------------------------------------------------------
# Read the biomass models (linear regression from Duan et al 2013 Tree Phys)
biomass = allometry.df
biomass$D = biomass$D / 10 # unit conversion: mm to cm
biomass$D_se = biomass$D_se / 10 # unit conversion: mm to cm

# Leaf dry mass model (No treatment effect)
biomass$LM = exp(0.615 * log((biomass$D)^2 * biomass$Height) + 0.36)
# biomass$LM_se = exp(0.615 * log(biomass$D_se * biomass$D_se * biomass$Height_se) + 0.36)
biomass$LM_se = ((biomass$D_se/biomass$D)^2 + (biomass$Height_se/biomass$Height)^2)^0.5 * biomass$LM

# Wood dry mass model (No treatment effect)
biomass$WM = exp(0.758 * log((biomass$D)^2 * biomass$Height) - 0.735)
# biomass$WM_se = exp(0.758 * log(biomass$D_se * biomass$D_se * biomass$Height_se) - 0.735)
biomass$WM_se = ((biomass$D_se/biomass$D)^2 + (biomass$Height_se/biomass$Height)^2)^0.5 * biomass$WM

# Root dry mass model (No treatment effect)
biomass$RM = exp(0.456 * log((biomass$D)^2 * biomass$Height) + 0.774)
# biomass$RM_se = exp(0.456 * log(biomass$D_se * biomass$D_se * biomass$Height_se) + 0.774)
biomass$RM_se = ((biomass$D_se/biomass$D)^2 + (biomass$Height_se/biomass$Height)^2)^0.5 * biomass$RM

# Calculate total plant mass
biomass$TM = biomass$LM + biomass$WM + biomass$RM
biomass$TM_se = ((biomass$LM_se)^2 + (biomass$WM_se)^2 + (biomass$RM_se)^2 )^0.5
biomass[,c(7:14)] = biomass[,c(7:14)] * c1 # unit conversion: g DM to g C

#---------------------------------------------------------------------------------------------
# Plot biomass over time and treatment
plots = list()
font.size = 12
pd <- position_dodge(0.75)
plots[[1]] = ggplot(biomass, aes(x=Date, y=LM, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=LM-LM_se, ymax=LM+LM_se), colour="grey", width=3) +
  geom_line(position=pd, data = biomass, aes(x = Date, y = LM, group = Water, colour=Water)) +
  ylab(expression(Leaf ~ mass ~ (g ~ C))) +
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

plots[[2]] = ggplot(biomass, aes(x=Date, y=WM, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=WM-WM_se, ymax=WM+WM_se), colour="grey", width=3) +
  geom_line(position=pd, data = biomass, aes(x = Date, y = WM, group = Water, colour=Water)) +
  ylab(expression(Wood ~ mass ~ (g ~ C))) +
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

plots[[3]] = ggplot(biomass, aes(x=Date, y=RM, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=RM-RM_se, ymax=RM+RM_se), colour="grey", width=3) +
  geom_line(position=pd, data = biomass, aes(x = Date, y = RM, group = Water, colour=Water)) +
  ylab(expression(Root ~ mass ~ (g ~ C))) +
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

plots[[4]] = ggplot(biomass, aes(x=Date, y=TM, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=TM-TM_se, ymax=TM+TM_se), colour="grey", width=3) +
  geom_line(position=pd, data = biomass, aes(x = Date, y = TM, group = Water, colour=Water)) +
  ylab(expression(Total ~ mass ~ (g ~ C))) +
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
png("output/4.Biomass_Duan_regression.png", units="px", width=1200, height=1000, res=100)
grid.arrange(grobs=plots, ncol=2)
dev.off()

#-------------------------------------------------------------------------------------
# Read in harvest data
harvest <- read.csv("raw_data/GHS30_Eglob-TxCxW_harvest_20110117-20110321_L1.csv")
harvest$Date <- as.Date(harvest$Date, format="%Y/%m/%d")
harvest$X = NULL
harvest = harvest[complete.cases(harvest), ]
harvest$D = (harvest$Basaldia1 + harvest$Basaldia2) / 2 # Calculate mean diameter
harvest$D = harvest$D / 10 # unit conversion: mm to cm
harvest$D2 = harvest$D * harvest$D
harvest$TotDW = harvest$LeafDW + harvest$RootDW + harvest$StemDW
harvest$D2H = harvest$D2 * harvest$Height
harvest = subset(harvest, Date <= as.Date("2011-02-27"))
#---------------------------------------------------------------------------------------------
# Plot harvest biomass over time and treatment
plots = list()
font.size = 12
pd <- position_dodge(1)
plots[[1]] = ggplot(data = harvest, aes(x = Date, y = LeafDW, group = Water, colour=Water)) + geom_point(position=pd) +
  ylab(expression(Leaf ~ mass ~ (g))) + labs(colour="Treatment") + theme_bw() +
  theme(legend.position = c(0.2,0.85), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(panel.grid = element_blank())

plots[[2]] = ggplot(data = harvest, aes(x = Date, y = StemDW, group = Water, colour=Water)) + geom_point(position=pd) +
  ylab(expression(Wood ~ mass ~ (g))) + theme_bw() +
  theme(legend.position = c(0.2,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(panel.grid = element_blank()) + guides(colour=FALSE)

plots[[3]] = ggplot(data = harvest, aes(x = Date, y = RootDW, group = Water, colour=Water)) + geom_point(position=pd) +
  ylab(expression(Root ~ mass ~ (g))) + theme_bw() +
  theme(legend.position = c(0.2,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(panel.grid = element_blank()) + guides(colour=FALSE)

plots[[4]] = ggplot(data = harvest, aes(x = Date, y = TotDW, group = Water, colour=Water)) + geom_point(position=pd) +
  ylab(expression(Total ~ mass ~ (g))) + theme_bw() +
  theme(legend.position = c(0.2,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(panel.grid = element_blank()) + guides(colour=FALSE)

# do.call(grid.arrange,  plots)
png("output/4.Harvest_over_time.png", units="px", width=1200, height=1000, res=100)
grid.arrange(grobs=plots, ncol=2)
dev.off()

#-------------------------------------------------------------------------------------
# Fit a linear regression for Leafmass by stem dia and height (ignoring temperature variation)
lm2 <- lm(log(LeafDW) ~ log(D) + log(Height), data=harvest)
lm1 <- lm(log(LeafDW) ~ log(D2H), data=harvest)

# Fit a linear regression by stem dia and height (ignoring temperature variation)
wm2 <- lm(log(StemDW) ~ log(D) + log(Height), data=harvest)
wm1 <- lm(log(StemDW) ~ log(D2H), data=harvest)
# sm2 <- lm(log(StemDW) ~ log(D) + log(Height), data=subset(harvest, Water %in% as.factor("Well watered")))

# Fit a linear regression for Rootmass by stem dia and height (ignoring temperature variation)
rm2 <- lm(log(RootDW) ~ log(D) + log(Height), data=harvest)
rm1 <- lm(log(RootDW) ~ log(D2H), data=harvest)
# rm3 <- lm(log(RootDW) ~ log(D), data=harvest)

# # Fit a linear regression for Leafarea by stem dia and height (ignoring temperature variation)
# la1 <- lm(log(Leafarea) ~ log(D) + log(Height), data=harvest)

#-----------------------------------------------------------------------------------------
# Plot predictions
layout(matrix(c(1,2,0,0,3,4),3,2,byrow=TRUE), widths=c(1,1), heights=c(10,1,10))
visreg(lm2, "D", overlay=TRUE)
visreg(lm2, "Height", overlay=TRUE)
visreg(lm1, "D2H", overlay=TRUE,legend=FALSE)
# visreg(lm2, "Height", overlay=TRUE)
plots[[1]] = recordPlot()

layout(matrix(c(1,2,0,0,3,4),3,2,byrow=TRUE), widths=c(1,1), heights=c(10,1,10))
visreg(wm2, "D", overlay=TRUE)
visreg(wm2, "Height", overlay=TRUE)
visreg(wm1, "D2H", overlay=TRUE,legend=FALSE)
# visreg(sm2, "Height", overlay=TRUE)
plots[[2]] = recordPlot()

layout(matrix(c(1,2,0,0,3,4),3,2,byrow=TRUE), widths=c(1,1), heights=c(10,1,10))
visreg(rm2, "D", overlay=TRUE)
visreg(rm2, "Height", overlay=TRUE)
visreg(rm1, "D2H", overlay=TRUE,legend=FALSE)
# visreg(rm3, "D", overlay=TRUE)
plots[[3]] = recordPlot()

pdf(file = "output/5.biomass_model_comparison.pdf")
plots[[1]]; plots[[2]]; plots[[3]]
dev.off()

#-----------------------------------------------------------------------------------------
# Save model sumary and stat comparison
sink("output/5.biomass_model_comparison.txt")
cat("\n\nLeafmass models:\n----------------\n### Linear regression with H and D:"); summary(lm1)
cat("\n### Linear regression with H and D^2:"); summary(lm2)
cat("### Comparison between both models:\n")
AIC(lm1, lm2); BIC(lm1, lm2)

cat("Stemmass models:\n----------------\n### Linear regression with H and D:"); summary(sm1)
cat("\n### Linear regression with H and D^2:"); summary(sm2)
cat("### Comparison between both models:\n")
AIC(wm1, wm2); BIC(wm1, wm2)

cat("\n\nRootmass models:\n----------------\n### Linear regression with H and D:");  summary(rm1)
cat("\n### Linear regression with H * D^2:"); summary(rm2)
cat("### Comparison between both models:\n")
AIC(rm1, rm2); BIC(rm1, rm2)
sink()
#-----------------------------------------------------------------------------------------

# Estimate the leaf mass from the fitted linear regression equation
eq = function(x,y){exp(coefficients(lm1)[1] + coefficients(lm1)[2] * log(x))  }
# Calculate all seedling leaf mass from height and diameter using the linear model
harvest$LM.modelled = eq(harvest$D2H)

# Estimate the leaf mass from the fitted linear regression equation
eq = function(x,y){exp(coefficients(wm1)[1] + coefficients(wm1)[2] * log(x))  }
# Calculate all seedling leaf mass from height and diameter using the linear model
harvest$WM.modelled = eq(harvest$D2H)

# Estimate the leaf mass from the fitted linear regression equation
eq = function(x,y){exp(coefficients(rm1)[1] + coefficients(rm1)[2] * log(x))  }
# Calculate all seedling leaf mass from height and diameter using the linear model
harvest$RM.modelled = eq(harvest$D2H)

harvest$TM.modelled = harvest$LM.modelled + harvest$WM.modelled + harvest$RM.modelled

#-----------------------------------------------------------------------------------------
# Plotting observation and modelled data
plots = list() 
par(mfrow = c(2, 2))
# Plotting observation and modelled Leaf mass data
plot(harvest$Height,harvest$LM.modelled,col="red",main="Height vs Leafmass", pch=15, xlab="Height (cm)", ylab="Leafmass (g DM)")
lines(harvest$Height,harvest$LeafDW,type="p", col="green", pch=20)
legend('topleft', c("Measurements", "Predicted"), lty=1, col=c('green','red'), bty='n', cex=0.75)

plot(harvest$D,harvest$LM.modelled, col="red", main="Diameter vs Leafmass", pch=15, xlab="Diameter (cm)", ylab="Leafmass (g DM)")
lines(harvest$D,harvest$LeafDW,type="p",col="green", pch=20)
legend('topleft', c("Measurements", "Predicted"), lty=1, col=c('green','red'), bty='n', cex=0.75)

# Plotting observation and modelled wood mass data
plot(harvest$Height,harvest$WM.modelled,col="red",main="Height vs Wood mass", pch=15, xlab="Height (cm)", ylab="Wood mass (g DM)")
lines(harvest$Height,harvest$StemDW,type="p", col="green", pch=20)
legend('topleft', c("Measurements", "Predicted"), lty=1, col=c('green','red'), bty='n', cex=0.75)

plot(harvest$D,harvest$WM.modelled, col="red", main="Diameter vs Wood mass", pch=15, xlab="Diameter (cm)", ylab="Wood mass (g DM)")
lines(harvest$D,harvest$StemDW,type="p",col="green", pch=20)
legend('topleft', c("Measurements", "Predicted"), lty=1, col=c('green','red'), bty='n', cex=0.75)
plots[[1]] = recordPlot()

# Plotting observation and modelled Root mass data
par(mfrow = c(2, 2))
plot(harvest$Height,harvest$RM.modelled,col="red",main="Height vs Rootmass", pch=15, xlab="Height (cm)", ylab="Rootmass (g DM)")
lines(harvest$Height,harvest$RootDW,type="p", col="green", pch=20)
legend('topleft', c("Measurements", "Predicted"), lty=1, col=c('green','red'), bty='n', cex=0.75)

plot(harvest$D,harvest$RM.modelled, col="red", main="Diameter vs Rootmass", pch=15, xlab="Diameter (cm)", ylab="Rootmass (g DM)")
lines(harvest$D,harvest$RootDW,type="p",col="green", pch=20)
legend('topleft', c("Measurements", "Predicted"), lty=1, col=c('green','red'), bty='n', cex=0.75)

# Plotting observation and modelled total mass data
plot(harvest$Height,harvest$TM.modelled,col="red",main="Height vs Total mass", pch=15, xlab="Height (cm)", ylab="Total mass (g DM)")
lines(harvest$Height,harvest$TotDW,type="p", col="green", pch=20)
legend('topleft', c("Measurements", "Predicted"), lty=1, col=c('green','red'), bty='n', cex=0.75)

plot(harvest$D,harvest$TM.modelled, col="red", main="Diameter vs Total mass", pch=15, xlab="Diameter (cm)", ylab="Total mass (g DM)")
lines(harvest$D,harvest$TotDW,type="p",col="green", pch=20)
legend('topleft', c("Measurements", "Predicted"), lty=1, col=c('green','red'), bty='n', cex=0.75)
plots[[2]] = recordPlot()


# Plotting observation and modelled data
# pdf(file = "output/5.biomass_measured_vs_predicted.pdf")
par(mfrow = c(2, 2))
# Plotting observation vs modelled Leaf mass data
plot(harvest$LeafDW,harvest$LM.modelled,col="green",main="LM", pch=15, xlab="LM data (g DM)", ylab="LM predicted (g DM)")
abline(0,1,col="red")

plot(harvest$StemDW,harvest$WM.modelled,col="green",main="WM", pch=15, xlab="WM data (g DM)", ylab="WM predicted (g DM)")
abline(0,1,col="red")

plot(harvest$RootDW,harvest$RM.modelled,col="green",main="RM", pch=15, xlab="RM data (g DM)", ylab="RM predicted (g DM)")
abline(0,1,col="red")

plot(harvest$TotDW,harvest$TM.modelled,col="green",main="TM", pch=15, xlab="TM data (g DM)", ylab="TM predicted (g DM)")
abline(0,1,col="red")
plots[[3]] = recordPlot()

# pdf(file = "output/5.biomass_measured_vs_predicted.pdf")
# plots[[1]]
# plots[[2]]
# plots[[3]]
# dev.off()


# Plot the original biomass data and compare with predicted ones
data.df = summaryBy(LeafDW+StemDW+RootDW+TotDW ~ Date+Water, data=harvest, FUN=c(mean,standard.error))
names(data.df) = c("Date", "Water", "LM", "WM", "RM", "TM", "LM_se", "WM_se", "RM_se", "TM_se")
melted.data.df = melt(data.df, id.vars=c("Date","Water"))
melted.data.df$Group = as.factor("Measured")

pred.df = summaryBy(LM.modelled+WM.modelled+RM.modelled+TM.modelled ~ Date+Water, data=harvest, FUN=c(mean,standard.error))
names(pred.df) = c("Date", "Water", "LM", "WM", "RM", "TM", "LM_se", "WM_se", "RM_se", "TM_se")
melted.pred.df = melt(pred.df, id.vars=c("Date","Water"))
melted.pred.df$Group = as.factor("Predicted")
melted.data = rbind(melted.data.df,melted.pred.df)

# plot all biomass data
plot = list()
meas = as.factor(c("LM", "WM", "RM", "TM"))
error = as.factor(c("LM_se", "WM_se", "RM_se", "TM_se"))
pd <- position_dodge(1) # move the overlapped errorbars horizontally
for (p in 1:length(meas)) {
  summary.data.Cpool = subset(melted.data,variable %in% meas[p])
  summary.error.Cpool = subset(melted.data,variable %in% error[p])
  summary.error.Cpool$parameter = summary.data.Cpool$value
  
  plot[[p]] = ggplot(summary.error.Cpool, aes(x=Date, y=parameter, group = interaction(Water,Group), colour=as.factor(Water), shape=as.factor(Group))) + 
    geom_point(position=pd,size=2.5) +
    geom_errorbar(position=pd,aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=2) +
    geom_line(position=pd,data = summary.error.Cpool, aes(x = Date, y = parameter, group = interaction(Water,Group), colour=as.factor(Water), linetype=as.factor(Group))) +
    ylab(paste(as.character(meas[p]),"(g DM)")) + 
    # xlab("Month") +
    # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g DM)")) +
    # ggtitle("C pools - Measured (points) vs Modelled (lines)") +
    labs(colour="Treatment",shape="Data Type",linetype="Data Type") +
    theme_bw() +
    # annotate("text", x = min(summary.error.Cpool$Date), y = max(summary.error.Cpool$value), size = 14, label = paste(title[p])) +
    # theme(plot.title = element_text(size = 20, face = "bold")) +
    theme(legend.title = element_text(colour="black", size=12)) +
    theme(legend.text = element_text(colour="black", size = 12)) +
    # theme(legend.key.height=unit(0.9,"line")) +
    theme(legend.position = c(0.2,0.80)) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=12)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = 14, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  if (p!=1) {
    plot[[p]] = plot[[p]] + guides(colour=FALSE, shape=FALSE, linetype=FALSE)
  }
}

pdf(file = "output/5.biomass_measured_vs_predicted.pdf",width=12, height=15)
# pdf(file = "output/5.biomass_measured_vs_predicted.pdf")
print (do.call(grid.arrange,  plot))
plots[[3]]
plots[[1]]
plots[[2]]
# grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]])
dev.off() 

#-----------------------------------------------------------------------------------------
# predict biomass using Kashif's models
biomass.final = allometry.df
biomass.final$D = biomass.final$D / 10 # unit conversion: mm to cm
biomass.final$D_se = biomass.final$D_se / 10 # unit conversion: mm to cm
biomass.final$D2H = biomass.final$D * biomass.final$D * biomass.final$Height
biomass.final$D2H_se = ((biomass.final$D_se / biomass.final$D)^2 + (biomass.final$D_se / biomass.final$D)^2 + 
                          (biomass.final$Height_se / biomass.final$Height)^2 )^0.5 * biomass.final$D2H

# Leaf dry mass model (No treatment effect)
eq = function(x,y){exp(coefficients(lm1)[1] + coefficients(lm1)[2] * log(x))  }
# Calculate all seedling leaf mass from height and diameter using the linear model
biomass.final$LM = eq(biomass.final$D2H)
biomass.final$LM_se = biomass.final$D2H_se/biomass.final$D2H * biomass.final$LM

# Wood dry mass model (No treatment effect)
eq = function(x,y){exp(coefficients(wm1)[1] + coefficients(wm1)[2] * log(x))  }
# Calculate all seedling leaf mass from height and diameter using the linear model
biomass.final$WM = eq(biomass.final$D2H)
biomass.final$WM_se = biomass.final$D2H_se/biomass.final$D2H * biomass.final$WM

# Root dry mass model (No treatment effect)
eq = function(x,y){exp(coefficients(rm1)[1] + coefficients(rm1)[2] * log(x))  }
# Calculate all seedling leaf mass from height and diameter using the linear model
biomass.final$RM = eq(biomass.final$D2H)
biomass.final$RM_se = biomass.final$D2H_se/biomass.final$D2H * biomass.final$RM

# Calculate total plant mass
biomass.final$TM = biomass.final$LM + biomass.final$WM + biomass.final$RM
biomass.final$TM_se = ((biomass.final$LM_se)^2 + (biomass.final$WM_se)^2 + (biomass.final$RM_se)^2 )^0.5
biomass.final[,c(9:16)] = biomass.final[,c(9:16)] * c1 # unit conversion: g DM to g C

#-----------------------------------------------------------------------------------------

# Plot biomass over time and treatment
plots = list()
font.size = 12
pd <- position_dodge(0.75)
plots[[1]] = ggplot(biomass.final, aes(x=Date, y=LM, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=LM-LM_se, ymax=LM+LM_se), colour="grey", width=3) +
  geom_line(position=pd, data = biomass.final, aes(x = Date, y = LM, group = Water, colour=Water)) +
  ylab(expression(Leaf ~ mass ~ (g ~ C))) +
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

plots[[2]] = ggplot(biomass.final, aes(x=Date, y=WM, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=WM-WM_se, ymax=WM+WM_se), colour="grey", width=3) +
  geom_line(position=pd, data = biomass.final, aes(x = Date, y = WM, group = Water, colour=Water)) +
  ylab(expression(Wood ~ mass ~ (g ~ C))) +
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

plots[[3]] = ggplot(biomass.final, aes(x=Date, y=RM, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=RM-RM_se, ymax=RM+RM_se), colour="grey", width=3) +
  geom_line(position=pd, data = biomass.final, aes(x = Date, y = RM, group = Water, colour=Water)) +
  ylab(expression(Root ~ mass ~ (g ~ C))) +
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

plots[[4]] = ggplot(biomass.final, aes(x=Date, y=TM, group = Water, colour=Water)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=TM-TM_se, ymax=TM+TM_se), colour="grey", width=3) +
  geom_line(position=pd, data = biomass.final, aes(x = Date, y = TM, group = Water, colour=Water)) +
  ylab(expression(Total ~ mass ~ (g ~ C))) +
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
png("output/5.Biomass_final.png", units="px", width=1200, height=1000, res=100)
grid.arrange(grobs=plots, ncol=2)
dev.off()

#-------------------------------------------------------------------------------------
# Save the biomass data for MCMC CBM
# write.csv(pred.data, file = "processed_data/modelled_data.csv", row.names = FALSE)
write.csv(biomass.final, file = "processed_data/biomass.final.csv", row.names = FALSE)

#-------------------------------------------------------------------------------------



