#- find a linear relationship between leaf-scale R and Tair
resp.fol <- read.csv("raw_data/WTC_TEMP_CM_GX-Rdark25_20130617-20140402_L2.csv")
names(resp.fol)[2] = "Date"
resp.fol$Date <- as.Date(resp.fol$Date)
resp.fol$T_treatment = as.factor( ifelse(resp.fol$T_treatment %in% as.factor("amb"), "ambient", "elevated") )
keeps <- c("Date", "chamber", "T_treatment", "Water_treatment", "Rmass", "Photo")
resp.fol = resp.fol[ , keeps, drop = FALSE]
resp.fol$Photo = -(resp.fol$Photo)
resp.fol = summaryBy(Rmass+Photo ~ Date+T_treatment, data=resp.fol, FUN=c(mean,standard.error))

# #- plot MASS BASED leaf R over time for the first date of high resolution T-response curves
font.size = 10
plots = list() 
pd <- position_dodge(0) # move the overlapped errorbars horizontally
plots[[1]] = ggplot(resp.fol, aes(x=Date, y=Rmass.mean, group = T_treatment, colour=T_treatment)) + 
  geom_point(position=pd, size=2) +
  geom_errorbar(position=pd, aes(ymin=Rmass.mean-Rmass.standard.error, ymax=Rmass.mean+Rmass.standard.error), colour="grey", width=1) +
  geom_line(position=pd,data = resp.fol, aes(x = Date, y = Rmass.mean, group = T_treatment)) +
  ylab(expression(R[leaf25]~(nmol~CO[2]~g^-1~s^-1))) + xlab("") +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(resp.fol$Date)-2, max(resp.fol$Date)+2)) +
  scale_color_manual(labels = c("ambient", "elevated"), values = c("blue", "red")) +
  labs(colour="Treatment") +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size=font.size)) +
  theme(legend.position = c(0.8,0.8), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size, family="Arial")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

#----------------------------------------------------------------------------------------------------------------
# import site weather data, take only Tair, format date stuff
files <- list.files(path = "raw_data/WTC_TEMP_CM_WTCMET", pattern = ".csv", full.names = TRUE)
temp <- lapply(files, fread, sep=",")
Tair <- rbindlist( temp )

Tair <- Tair[ , c("chamber","DateTime","Tair_al")]
Tair$Date <- as.Date(Tair$DateTime)

# need to turn the datetime into hms
Tair$DateTime <- ymd_hms(Tair$DateTime)
Tair$time <- format(Tair$DateTime, format='%H:%M:%S')

# met.data.na = met.data[is.na(met.data$SoilTemp),] # Check any NA values for soil temperature

# subset by Date range of experiment
# Tair <- subset(Tair[, c("chamber","Date","time","Tair_al")], Date  >= "2013-09-17" & Date  <= "2014-05-26")
Tair$chamber = as.factor(Tair$chamber)
Tair = merge(Tair, unique(height.dia[,c("chamber","T_treatment")]), by="chamber")

###----------------------------------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------------------------------
# Temperature dependant parameter variability
# Find the previous 3-day mean daily Tair
Tair.sub = subset(Tair[, c("chamber","T_treatment","Date","time","Tair_al")], Date  >= "2013-09-14" & Date  <= "2014-05-26")

prec.day.no = 3
day = data.frame(Date = as.Date(1:(prec.day.no*(length(unique(Tair.sub$Date))-3)), origin=Sys.Date()), iter=numeric(prec.day.no*(length(unique(Tair.sub$Date))-3)),
                 Date.resp = as.Date(1:(prec.day.no*(length(unique(Tair.sub$Date))-3)), origin=Sys.Date()))
for (i in 4:length(unique(Tair.sub$Date))) {
  day$Date[(1+(i-4)*prec.day.no):((i-3)*prec.day.no)] = as.Date((unique(Tair.sub$Date)-prec.day.no)[i] : (unique(Tair.sub$Date)-1)[i])
  day$iter[(1+(i-4)*prec.day.no):((i-3)*prec.day.no)] = i-3
  day$Date.resp[(1+(i-4)*prec.day.no):((i-3)*prec.day.no)] = as.Date(unique(Tair.sub$Date)[i])
}

# Tair.sub = subset(Tair.sub[, c("chamber","T_treatment","Date","time","Tair_al")], Date %in% day$Date)
Tair.sub = merge(Tair.sub, day, by="Date", allow.cartesian=TRUE)
Tair.sum = summaryBy(Tair_al ~ Date.resp+T_treatment, data=Tair.sub, FUN=c(mean), na.rm=TRUE)
names(Tair.sum)[c(1,3)] = c("Date","Tair")

# write csv file with daily inputs of GPP, Ra, LA, Tair
write.csv(Tair.sum, file = "processed_data/Tair.csv", row.names = FALSE)

# Merge Tair with other data set
data.all = merge(data.all, Tair.sum, by=c("Date","T_treatment"), all=TRUE)

###----------------------------------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------------------------------

prec.day.no = 3
day = data.frame(Date = as.Date(1:(prec.day.no*length(unique(resp.fol$Date))), origin=Sys.Date()), iter=numeric(prec.day.no*length(unique(resp.fol$Date))),
                 Date.resp = as.Date(1:(prec.day.no*length(unique(resp.fol$Date))), origin=Sys.Date()))
for (i in 1:length(unique(resp.fol$Date))) {
  day$Date[(1+(i-1)*prec.day.no):(i*prec.day.no)] = as.Date((unique(resp.fol$Date)-prec.day.no)[i] : (unique(resp.fol$Date)-1)[i])
  day$iter[(1+(i-1)*prec.day.no):(i*prec.day.no)] = i
  day$Date.resp[(1+(i-1)*prec.day.no):(i*prec.day.no)] = as.Date(unique(resp.fol$Date)[i])
}

Tair.sub = subset(Tair[, c("chamber","T_treatment","Date","time","Tair_al")], Date %in% day$Date)
Tair.sub = merge(Tair.sub, day, by="Date")
Tair.sum = summaryBy(Tair_al ~ Date.resp+T_treatment, data=Tair.sub, FUN=c(mean), na.rm=TRUE)
names(Tair.sum)[1] = "Date"
resp.fol.final = merge(resp.fol, Tair.sum, by=c("Date","T_treatment"))

Rleaf25 <- lm(Rmass.mean~Tair_al.mean,data=resp.fol.final)
summary(Rleaf25)

plots[[2]] = ggplot(data = resp.fol.final, aes(x=Tair_al.mean, y=Rmass.mean)) + 
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_smooth(method="lm",se=FALSE) + geom_point(size=0.75) + 
  ylab(expression(R[leaf25]~(nmol~CO[2]~g^-1~s^-1))) + xlab(expression("Previous 3-day mean daily"~T[air]~(degree*C))) + 
  theme_bw() + theme(legend.key.width=unit(0.9,"line")) +
  theme(text = element_text(size=font.size, family="Arial")) +
  theme(axis.title.x = element_text(size = font.size, vjust=0.3)) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

prec.day.no = 7
day = data.frame(Date = as.Date(1:(prec.day.no*length(unique(resp.fol$Date))), origin=Sys.Date()), iter=numeric(prec.day.no*length(unique(resp.fol$Date))),
                 Date.resp = as.Date(1:(prec.day.no*length(unique(resp.fol$Date))), origin=Sys.Date()))
for (i in 1:length(unique(resp.fol$Date))) {
  day$Date[(1+(i-1)*prec.day.no):(i*prec.day.no)] = as.Date((unique(resp.fol$Date)-prec.day.no)[i] : (unique(resp.fol$Date)-1)[i])
  day$iter[(1+(i-1)*prec.day.no):(i*prec.day.no)] = i
  day$Date.resp[(1+(i-1)*prec.day.no):(i*prec.day.no)] = as.Date(unique(resp.fol$Date)[i])
}

Tair.sub = subset(Tair[, c("chamber","T_treatment","Date","time","Tair_al")], Date %in% day$Date)
Tair.sub = merge(Tair.sub, day, by="Date")
Tair.sum = summaryBy(Tair_al ~ Date.resp+T_treatment, data=Tair.sub, FUN=c(mean), na.rm=TRUE)
names(Tair.sum)[1] = "Date"
resp.fol.final = merge(resp.fol, Tair.sum, by=c("Date","T_treatment"))

Rleaf25 <- lm(Rmass.mean~Tair_al.mean,data=resp.fol.final)
summary(Rleaf25)

plots[[3]] = ggplot(data = resp.fol.final, aes(x=Tair_al.mean, y=Rmass.mean)) + 
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_smooth(method="lm",se=FALSE) + geom_point(size=0.75) + 
  ylab(expression(R[leaf25]~(nmol~CO[2]~g^-1~s^-1))) + xlab(expression("Previous 7-day mean daily"~T[air]~(degree*C))) + 
  theme_bw() + theme(legend.key.width=unit(0.9,"line")) +
  theme(text = element_text(size=font.size, family="Arial")) +
  theme(axis.title.x = element_text(size = font.size, vjust=0.3)) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


prec.day.no = 1
day = data.frame(Date = as.Date(1:(prec.day.no*length(unique(resp.fol$Date))), origin=Sys.Date()), iter=numeric(prec.day.no*length(unique(resp.fol$Date))),
                 Date.resp = as.Date(1:(prec.day.no*length(unique(resp.fol$Date))), origin=Sys.Date()))
for (i in 1:length(unique(resp.fol$Date))) {
  day$Date[(1+(i-1)*prec.day.no):(i*prec.day.no)] = as.Date((unique(resp.fol$Date)-prec.day.no)[i] : (unique(resp.fol$Date)-1)[i])
  day$iter[(1+(i-1)*prec.day.no):(i*prec.day.no)] = i
  day$Date.resp[(1+(i-1)*prec.day.no):(i*prec.day.no)] = as.Date(unique(resp.fol$Date)[i])
}

Tair.sub = subset(Tair[, c("chamber","T_treatment","Date","time","Tair_al")], Date %in% day$Date)
Tair.sub = merge(Tair.sub, day, by="Date")
Tair.sum = summaryBy(Tair_al ~ Date.resp+T_treatment, data=Tair.sub, FUN=c(mean), na.rm=TRUE)
names(Tair.sum)[1] = "Date"
resp.fol.final = merge(resp.fol, Tair.sum, by=c("Date","T_treatment"))

Rleaf25 <- lm(Rmass.mean~Tair_al.mean,data=resp.fol.final)
summary(Rleaf25)

plots[[4]] = ggplot(data = resp.fol.final, aes(x=Tair_al.mean, y=Rmass.mean)) + 
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_smooth(method="lm",se=FALSE) + geom_point(size=0.75) + 
  ylab(expression(R[leaf25]~(nmol~CO[2]~g^-1~s^-1))) + xlab(expression("Previous 1-day mean daily"~T[air]~(degree*C))) + 
  theme_bw() + theme(legend.key.width=unit(0.9,"line")) +
  theme(text = element_text(size=font.size, family="Arial")) +
  theme(axis.title.x = element_text(size = font.size, vjust=0.3)) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

Tair.sub = subset(Tair[, c("chamber","T_treatment","Date","time","Tair_al")], Date %in% resp.fol$Date)
Tair.sum = summaryBy(Tair_al ~ Date+T_treatment, data=Tair.sub, FUN=c(mean), na.rm=TRUE)
resp.fol.final = merge(resp.fol, Tair.sum, by=c("Date","T_treatment"))

Rleaf25 <- lm(Rmass.mean~Tair_al.mean,data=resp.fol.final)
summary(Rleaf25)

plots[[5]] = ggplot(data = resp.fol.final, aes(x=Tair_al.mean, y=Rmass.mean)) + 
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_smooth(method="lm",se=FALSE) + geom_point(size=0.75) + 
  ylab(expression(R[leaf25]~(nmol~CO[2]~g^-1~s^-1))) + xlab(expression("Mean daily"~T[air]~(degree*C))) + 
  theme_bw() + theme(legend.key.width=unit(0.9,"line")) +
  theme(text = element_text(size=font.size, family="Arial")) +
  theme(axis.title.x = element_text(size = font.size, vjust=0.3)) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

prec.day.no = 3
day = data.frame(Date = as.Date(1:(prec.day.no*length(unique(resp.fol$Date))), origin=Sys.Date()), iter=numeric(prec.day.no*length(unique(resp.fol$Date))),
                 Date.resp = as.Date(1:(prec.day.no*length(unique(resp.fol$Date))), origin=Sys.Date()))
for (i in 1:length(unique(resp.fol$Date))) {
  day$Date[(1+(i-1)*prec.day.no):(i*prec.day.no)] = as.Date((unique(resp.fol$Date)-prec.day.no+1)[i] : (unique(resp.fol$Date))[i])
  day$iter[(1+(i-1)*prec.day.no):(i*prec.day.no)] = i
  day$Date.resp[(1+(i-1)*prec.day.no):(i*prec.day.no)] = as.Date(unique(resp.fol$Date)[i])
}

Tair.sub = subset(Tair[, c("chamber","T_treatment","Date","time","Tair_al")], Date %in% day$Date)
Tair.sub = merge(Tair.sub, day, by="Date")
Tair.sum = summaryBy(Tair_al ~ Date.resp+T_treatment, data=Tair.sub, FUN=c(mean), na.rm=TRUE)
names(Tair.sum)[1] = "Date"
resp.fol.final = merge(resp.fol, Tair.sum, by=c("Date","T_treatment"))

Rleaf25 <- lm(Rmass.mean~Tair_al.mean,data=resp.fol.final)
summary(Rleaf25)

plots[[6]] = ggplot(data = resp.fol.final, aes(x=Tair_al.mean, y=Rmass.mean)) + 
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_smooth(method="lm",se=FALSE) + geom_point(size=0.75) + 
  ylab(expression(R[leaf25]~(nmol~CO[2]~g^-1~s^-1))) + xlab(expression("3-day mean daily"~T[air]~(degree*C))) + 
  theme_bw() + theme(legend.key.width=unit(0.9,"line")) +
  theme(text = element_text(size=font.size, family="Arial")) +
  theme(axis.title.x = element_text(size = font.size, vjust=0.3)) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("output/7.Rf25_vs_Tair.png", units="px", width=2000, height=1000, res=130)
do.call(grid.arrange,  plots)
dev.off()

do.call(grid.arrange,  plots)

#----------------------------------------------------------------------------------------------------------------
prec.day.no = 3
day = data.frame(Date = as.Date(rep(as.Date("2013-09-17"):as.Date("2014-05-26"),prec.day.no)), iter=numeric(252*prec.day.no),
                 Date.resp = as.Date(rep(as.Date("2013-09-17"):as.Date("2014-05-26"),prec.day.no)))
day = day[with(day, order(Date.resp)), ]
for (i in 1:252) {
  day$Date[(1+(i-1)*prec.day.no):(i*prec.day.no)] = as.Date((day$Date.resp[i*prec.day.no]-prec.day.no+1) : day$Date.resp[i*prec.day.no])
  day$iter[(1+(i-1)*prec.day.no):(i*prec.day.no)] = i
  # day$Date.resp[(1+(i-1)*prec.day.no):(i*prec.day.no)] = as.Date(unique(resp.fol$Date)[i])
}
Tair.sub = subset(Tair[, c("chamber","T_treatment","Date","time","Tair_al")], Date %in% day$Date)
Tair.sub$repeats = ifelse(Tair.sub$Date %in% as.Date(c("2013-09-15","2014-05-26")),1,
                          ifelse(Tair.sub$Date %in% as.Date(c("2013-09-16","2014-05-25")),2,3))
Tair.sub = Tair.sub[rep(seq(nrow(Tair.sub)), Tair.sub$repeats),]

Tair.sub = unique(merge(Tair.sub, day, by="Date", allow.cartesian=TRUE))
Tair.final = summaryBy(Tair_al ~ Date.resp+T_treatment, data=Tair.sub, FUN=c(mean), na.rm=TRUE)
names(Tair.final)[1] = "Date"
# Tair.final = merge(resp.fol, Tair.sum, by=c("Date","T_treatment"))

# Tair.final = subset(Tair[, c("chamber","T_treatment","Date","time","Tair_al")], Date  >= "2013-09-17" & Date  <= "2014-05-26")
# Tair.final = summaryBy(Tair_al ~ Date+T_treatment, data=Tair.final, FUN=c(mean), na.rm=TRUE)
Tair.final$rd25.foliage = coef(Rleaf25)[1] + (coef(Rleaf25)[2] * Tair.final$Tair_al.mean)
Tair.final$rd25.foliage = Tair.final$rd25.foliage * (10^-9 * 12) * (3600 * 24) # unit conversion from nmolCO2 g-1 s-1 to gC gC-1 d-1
# Tair.final$rd25.foliage = Tair.final$rd25.foliage * (10^-9 * 12) * (3600 * 24) * (1/c1) # unit conversion from nmolCO2 g-1 s-1 to gC gC-1 d-1
