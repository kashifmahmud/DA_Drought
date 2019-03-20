# Script to read and setup a model for GPP, Ra (aboveground respiration) and leaf area
#  read the hourly flux data
data.hr <- read.csv("raw_data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
data.hr$DateTime <- as.POSIXct(data.hr$DateTime,format="%Y-%m-%d %H:%M:%S",tz="GMT")
data.hr$Date <- as.Date(data.hr$DateTime)

#- partition the net fluxes into GPP and Ra components. Specify the actitation energy and the number
#-  of prior nights to include for the estimate of the basal respiration rate (lagdates)
data.hr.p <- partitionHourlyFluxCUE_arr(data.hr.gf=data.hr,Ea=57.69,lagdates=3)

#- get daily sums from the partitioned hourly data
cue.list <- returnCUE.day(dat=data.hr.p) # get daily sums from hourly data
cue.day <- cue.list[[1]]                # extract chamber values on each day
cue.day.trt <- cue.list[[2]]            # extract treatment averages

data=cue.day.trt
names(data)[3:8] = c("GPP","Ra","LA","GPP_SE","Ra_SE","LA_SE")
# ignore the first few days of data to start on 2013-09-17 from where we have H & D measurements
# data = subset(data, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-02-12"))
data = subset(data, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-05-27"))

# # write csv file with daily inputs of GPP, Ra, LA
# write.csv(data, file = "processed_data/data_GPP_Ra_LA.csv", row.names = FALSE)

#- plot GPP, Ra, and LA data over time for various treatments
font.size = 12
plots = list()
pd <- position_dodge(0)
plots[[1]] = ggplot(data=data, aes(x=Date, y=GPP, group = interaction(T_treatment), colour=T_treatment)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=GPP-GPP_SE, ymax=GPP+GPP_SE), colour="grey", width=1) +
  geom_line(position=pd, data = data, aes(x = Date, y = GPP, group = interaction(T_treatment), colour=T_treatment)) +
  ylab(expression(GPP~"(g C "*d^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data$Date), max(data$Date))) +
  labs(colour="Temperature") +
  scale_color_manual(labels = c("ambient", "warmed"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.2,0.85), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

plots[[2]] = ggplot(data=data, aes(x=Date, y=Ra, group = interaction(T_treatment), colour=T_treatment)) + 
  geom_point(position=pd) +
  geom_errorbar(position=pd, aes(ymin=Ra-Ra_SE, ymax=Ra+Ra_SE), colour="grey", width=1) +
  geom_line(position=pd, data = data, aes(x = Date, y = Ra, group = interaction(T_treatment), colour=T_treatment)) +
  ylab(expression(R[above]~"(g C "*d^"-1"*")")) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month",limits = c(min(data$Date), max(data$Date))) +
  labs(colour="Temperature") +
  scale_color_manual(labels = c("ambient", "warmed"), values = c("blue", "red")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size)) +
  theme(legend.position = c(0.2,0.85), legend.box = "horizontal") + theme(legend.key.height=unit(0.8,"line")) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

png("output/1.GPP_Ra_over_time.png", units="px", width=2000, height=1000, res=130)
do.call(grid.arrange,  plots)
dev.off()

do.call(grid.arrange,  plots)
