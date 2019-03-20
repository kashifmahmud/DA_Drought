
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Statistical analysis of Ra, GPP, and Ra/GPP (Table 1).
#- Note that many of these analyses take considerable time to run. The entire script may take 10-20 minutes on a typical machine.

#- This script assumes that the dataframe "cue.day" is available. Run "main_script.R" up to line 102 to create this dataframe.


#- Here is a manual calculation of the expected degrees of freedom.

# T = number of timepoints. w= number of warming treatments. c= number of replicate chambers within each treatment

# Source                df calculation               df  
# Warming                  w-1                        1
# chamber(Warming)         w*(c-1)                    10   # this is the random whole-plot error term to test warming main effect
# Time                     T-1                     254     # in this case I have 255 timepoints
# Time*Warming            (T-1)*(w-1)              254
# Time*chamber[warming]   (T-1)*(c-1)*w             2540   # this is the random sub-plot error term. it's the same as the residual
#                                                               , as there is no sub-replication here.. 254*5*2
#- note the random sub-plot error term (i.e, the residual) will actually have fewer df than this, as I am excluding the drought data



dat2 <- subset(cue.day,Water_treatment=="control")
dat2$T_treatment <- as.factor(dat2$T_treatment)
dat2$DateFac <- as.factor(dat2$Date)










#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Ra

#- find the "best" Box-Cox transformation for Ra_la
sp.Ra.simple <- lm(Ra_la~T_treatment*DateFac,data=dat2)
value.r <- boxcox(object=sp.Ra.simple,lambda=seq(-1,0,1/20))
exponent.r <- value.r$x[which.max(value.r$y)]
#dat2$Ra_latrans <- BoxCox(dat2$Ra_la,lambda=exponent.r)
dat2$Ra_latrans <- log(exp(dat2$Ra_la^0.5))  # extensive transformation to make data amenable to ANOVA assumptions
hist(dat2$Ra_latrans)

dat2$leafArea_level <- cut(dat2$leafArea,breaks=6)#- variance is far from constant- it varies with leaf area.  Let's model it.
boxplot(Ra_latrans~leafArea_level,data=dat2)
dat2$leafArea_inv <- 1/dat2$leafArea              #- the inverse of leaf area is used to model the variance for Ra below
dat2$Tlevel <- cut(dat2$Tair_24hrs,breaks=6)      #- create 6 bins of temperature. Used to model variance for Ra/GPP below



#- smaller subset of data, just to make the models run faster for testing
dat3 <- subset(dat2,Date < as.Date("2013-12-1"))
dat3$DateFac <- as.factor(dat3$Date)
sp.test <- lme(Ra_latrans~T_treatment*DateFac,random=list(~1|chamber),
             weights=varPower(value=0.2,form=~leafArea_inv),
             data=dat3,method="REML")
plot(sp.test,Ra_latrans~fitted(.),abline=c(0,1))                   #predicted vs. fitted
qqnorm(sp.test, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot. Not perfect, but a huge improvement.


#- fit model to all data, and
#- compare models with and without autocorrelation
sp.ra <- lme(Ra_latrans~T_treatment*DateFac,random=list(~1|chamber),
             weights=varPower(value=0.2,form=~leafArea_inv),data=dat2,method="ML")
sp.ra.ar1 <- update(sp.ra,correlation=corAR1(form=~1|chamber),method="ML")
AIC(sp.ra,sp.ra.ar1)
anova(sp.ra,sp.ra.ar1) # model with autocorrelation is immensely better!

#- refit best model with REML
sp.ra.ar1.reml <- update(sp.ra.ar1,method="REML")

#look at model diagnostics
plot(sp.ra.ar1.reml,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.ra.ar1.reml,Ra_latrans~fitted(.)|chamber,abline=c(0,1))           #predicted vs. fitted for each chamber
plot(sp.ra.ar1.reml,Ra_latrans~fitted(.),abline=c(0,1))                   #predicted vs. fitted
qqnorm(sp.ra.ar1.reml, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot. Departure at high values.
anova(sp.ra.ar1.reml,type="marginal")
Anova(sp.ra.ar1.reml)

#- compare models with and without interaction terms
sp.ra.ar1.noint <- update(sp.ra.ar1,.~.-T_treatment:DateFac)
sp.ra.ar1.reml.noint <- update(sp.ra.ar1.reml,.~.-T_treatment:DateFac)

AIC(sp.ra.ar1,sp.ra.ar1.noint)   # dropping the interaction results in a more parsimonious model (huge reduction in df)
anova(sp.ra.ar1,sp.ra.ar1.noint) # The interaction term should be removed. 

# ANOVA for Ra
anova(sp.ra.ar1.reml,type="marginal") 
anova(sp.ra.ar1.reml.noint,type="marginal") 

lsRa <- summary(lsmeans::lsmeans(sp.ra.ar1.reml.noint,"T_treatment"))
(lsRa$lsmean[1]-lsRa$lsmean[2])/lsRa$lsmean[1]*100 # percentage change in response to warming

#- estimate explainatory power (r2 values) of models with and without the warming by time interaction

sp.ra.ar1.reml.noTrt <- update(sp.ra.ar1.noint,.~.-T_treatment)
sp.ra.ar1.reml.noDate <- update(sp.ra.ar1.noint,.~.-DateFac)

ra.full <- r.squaredGLMM(sp.ra.ar1.reml)[1]      # full model
ra.noint <- r.squaredGLMM(sp.ra.ar1.reml.noint)[1] # model lacking interaction
ra.notrt <- r.squaredGLMM(sp.ra.ar1.reml.noTrt)[1] # model lacking treatment and interaction
ra.nodate <- r.squaredGLMM(sp.ra.ar1.reml.noDate)[1]

ra.full
ra.full-ra.noint
ra.full-ra.notrt
ra.full-ra.nodate
#- Warming by date interaction for Ra_la should be removed, and excluding it drops the explanatory power of the model
#    very modestly (marginal r2 was 0.504 with the interaction, 0.507 without the interaction).. 
#- Therefore the warming by date interaction is not important in a quantatiative sense, and it should be removed.
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------












#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- GPP

#- again, fit small subset of the data for testing purposes
sp.test2 <- lme(log(GPP_la)~T_treatment*DateFac,random=list(~1|chamber),
               weights=varPower(value=0.2,form=~leafArea_inv),
               data=dat3,method="REML")
plot(sp.test2,GPP_la~fitted(.),abline=c(0,1))                       #predicted vs. fitted
qqnorm(sp.test2, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot. Not perfect, but the log gives an improvement



#- fit model to all data, and
#- compare models with and without autocorrelation
sp.gpp <- lme(log(GPP_la)~T_treatment*DateFac,random=list(~1|chamber),
              weights=varPower(value=0.2,form=~leafArea_inv),data=dat2,method="ML")
sp.gpp.ar1 <- update(sp.gpp,correlation=corAR1(0.5,form=~1|chamber),method="ML")
AIC(sp.gpp,sp.gpp.ar1)
anova(sp.gpp,sp.gpp.ar1) # model with autocorrelation is immensely better!

#- refit best model with REML
sp.gpp.ar1.reml <- update(sp.gpp.ar1,method="REML")

#look at model diagnostics
plot(sp.gpp.ar1.reml,resid(.,type="p")~fitted(.) | T_treatment,abline=0)        #resid vs. fitted for each treatment
plot(sp.gpp.ar1.reml,log(GPP_la)~fitted(.)|chamber,abline=c(0,1))               #predicted vs. fitted for each chamber
plot(sp.gpp.ar1.reml,log(GPP_la)~fitted(.),abline=c(0,1))                       #predicted vs. fitted
qqnorm(sp.gpp.ar1.reml, ~ resid(., type = "p"), abline = c(0, 1))               #qqplot
hist(sp.gpp.ar1.reml$residuals)



#- compare models with and without interaction terms
sp.gpp.ar1.noint <- update(sp.gpp.ar1,.~.-T_treatment:DateFac)
sp.gpp.ar1.reml.noint <- update(sp.gpp.ar1.reml,.~.-T_treatment:DateFac)

anova(sp.gpp.ar1,sp.gpp.ar1.noint) # dropping the interaction results in a poorer model. The interaction should stay.
anova(sp.gpp.ar1.reml,type="marginal")
anova(sp.gpp.ar1.reml.noint,type="marginal")  # note that there is a huge main effect of warming on GPP if the interaction is removed
                                              # so on average, warming reduced GPP. But not equally on ALL dates.

#- estimate explainatory power (r2 values) of models with and without the warming by time interaction
sp.gpp.ar1.reml.noTrt <- update(sp.gpp.ar1.reml.noint,.~.-T_treatment)
sp.gpp.ar1.reml.noDate <- update(sp.gpp.ar1.reml.noint,.~.-DateFac)

gpp.full <- r.squaredGLMM(sp.gpp.ar1.reml)[1]      # full model
gpp.noint <- r.squaredGLMM(sp.gpp.ar1.reml.noint)[1] # model lacking interaction
gpp.notrt <- r.squaredGLMM(sp.gpp.ar1.reml.noTrt)[1] # model lacking treatment and interaction
gpp.nodate <- r.squaredGLMM(sp.gpp.ar1.reml.noDate)[1] # model lacking date eraction

gpp.full
gpp.full-gpp.noint
gpp.full-gpp.notrt
gpp.full-gpp.nodate
#- so the interaction between T_treatment and date is "significant" for GPP, and it is somewhat important.
#- Marginal r2 drops a small amount when excluding the interaction (0.902 with interaction, 0.894 without)

lsGpp <- summary(lsmeans::lsmeans(sp.gpp.ar1.reml,"T_treatment"))
(lsGpp$lsmean[1]-lsGpp$lsmean[2])/lsGpp$lsmean[1]*100 # percentage change in response to warming


# extract and attempt to understand the GPP warming x date interaction
lsGpp2 <- summary(lsmeans::lsmeans(sp.gpp.ar1.reml,specs=c("T_treatment","DateFac")))
lsGppa <- subset(lsGpp2,T_treatment=="ambient")
lsGppe <- subset(lsGpp2,T_treatment=="elevated")
diff <- data.frame(DateFac = lsGppa$DateFac,diff = (lsGppa$lsmean-lsGppe$lsmean)/lsGppa$lsmean)
diff2 <- merge(diff,subset(dat2,T_treatment=="ambient"))
plot(diff~PAR,data=diff2);abline(h=0)  # warming effect is quite consistent on high light days, variable on low light days. 
               

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------







#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- Ra/GPP. 

#- find the "right" transformation
sp.CUE.simple <- lm(RtoA~T_treatment*DateFac,data=dat2)
value <- boxcox(object=sp.CUE.simple,lambda=seq(0,2,1/20))
exponent.cue <- value$x[which.max(value$y)]

dat3$RtoAtrans <- BoxCox(dat3$RtoA,lambda=exponent.cue)               # box-cox
dat2$RtoAtrans <- BoxCox(dat2$RtoA,lambda=exponent.cue)
#dat3$RtoAtrans <- with(dat3,log(RtoA/(1-RtoA))) # log ratio
#dat3$RtoAtrans <- with(dat3,asin(RtoA^0.5))      # arcsin-root transformation


hist(dat3$RtoAtrans)
boxplot(RtoA~Tlevel,data=dat3)

sp.test3<- lme(RtoAtrans~T_treatment*DateFac,random=list(~1|chamber),
               weights=varIdent(form=~Tlevel),
                data=dat3,method="REML")
plot(sp.test3,RtoAtrans~fitted(.),abline=c(0,1))                   #predicted vs. fitted
qqnorm(sp.test3, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot. Not perfect, but a huge improvement.
hist(sp.test3$residuals)






#- compare models with and without autocorrelation. These models take quite a long time to fit.
sp.cue <- lme(RtoAtrans~T_treatment*DateFac,random=list(~1|chamber),
              weights=varIdent(form=~1|Tlevel),
              data=dat2,method="ML")
sp.cue.ar1 <- update(sp.cue,correlation=corAR1(form=~1|chamber),method="ML")
AIC(sp.cue,sp.cue.ar1)
anova(sp.cue,sp.cue.ar1) # model with autocorrelation is immensely better!


#- refit best model with REML
sp.cue.ar1.reml <- update(sp.cue.ar1,method="REML")

#look at model diagnostics
plot(sp.cue.ar1.reml,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.cue.ar1.reml,RtoAtrans~fitted(.)|chamber,abline=c(0,1))            #predicted vs. fitted for each chamber
plot(sp.cue.ar1.reml,RtoAtrans~fitted(.),abline=c(0,1))                    #predicted vs. fitted
qqnorm(sp.cue.ar1.reml, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot. Ugh. Departure at both ends
hist(sp.cue.ar1.reml$residuals)


#- exclude outliers and refit
outliers <- unname(which(abs(residuals(sp.cue.ar1.reml,type="normalized"))>qnorm(0.985)))
noout <- dat2[-outliers,]
sp.cue.ar1.reml.noout <- lme(RtoAtrans~T_treatment*DateFac,random=list(~1|chamber),
                             correlation=corAR1(form=~1|chamber),method="REML",
                             weights=varIdent(form=~1|Tlevel),
                             data=noout)
plot(sp.cue.ar1.reml.noout,resid(.,type="p")~fitted(.) | T_treatment,abline=0)  #resid vs. fitted for each treatment
plot(sp.cue.ar1.reml.noout,RtoAtrans~fitted(.)|chamber,abline=c(0,1))           #predicted vs. fitted for each chamber
plot(sp.cue.ar1.reml.noout,RtoAtrans~fitted(.),abline=c(0,1))                   #predicted vs. fitted
qqnorm(sp.cue.ar1.reml.noout, ~ resid(., type = "p"), abline = c(0, 1))         #qqplot. A bit better.

#- we arrive at the same inference whether outliers are included or removed.
#- So, the non-normal residuals may not have strongly affected the statistical inference.
anova(sp.cue.ar1.reml,type="marginal")
anova(sp.cue.ar1.reml.noout,type="marginal")

lsRtoA <- summary(lsmeans::lsmeans(sp.cue.ar1.reml,"T_treatment"))
(lsRtoA$lsmean[1]-lsRtoA$lsmean[2])/lsRtoA$lsmean[1]*100 # percentage change in response to warming


# extract and attempt to understand the GPP warming x date interaction
lsRtoA2 <- summary(lsmeans::lsmeans(sp.cue.ar1.reml,specs=c("T_treatment","DateFac")))
lsRtoAa <- subset(lsRtoA2,T_treatment=="ambient")
lsRtoAe <- subset(lsRtoA2,T_treatment=="elevated")
diffRtoA <- data.frame(DateFac = lsRtoAa$DateFac,diff = (lsRtoAa$lsmean-lsRtoAe$lsmean)/lsRtoAa$lsmean)
diffRtoA2 <- merge(diffRtoA,subset(dat2,T_treatment=="ambient"))

plot(diff~VPDair,data=subset(diffRtoA2,diff>-2));abline(h=0)  # warming effect tends to increase with temperature
                  

#- compare models with and without interaction terms
sp.cue.ar1.noint <- update(sp.cue.ar1,.~.-T_treatment:DateFac)
sp.cue.ar1.reml.noint <- update(sp.cue.ar1.reml,.~.-T_treatment:DateFac)
anova(sp.cue.ar1,sp.cue.ar1.noint) # dropping the interaction results in a poorer model. The interaction is important.

#- estimate explainatory power (r2 values) of models with and without the warming by time interaction
sp.cue.ar1.reml.noTrt <- update(sp.cue.ar1.reml.noint,.~.-T_treatment)
sp.cue.ar1.reml.noDate <- update(sp.cue.ar1.reml.noint,.~.-DateFac)

cue.full <- r.squaredGLMM(sp.cue.ar1.reml)[1]      # full model
cue.noint <- r.squaredGLMM(sp.cue.ar1.reml.noint)[1] # model lacking interaction
cue.notrt <- r.squaredGLMM(sp.cue.ar1.reml.noTrt)[1] # model lacking treatment and interaction
cue.nodate <- r.squaredGLMM(sp.cue.ar1.reml.noDate)[1] # model lacking date or interaction

cue.full
cue.full-cue.noint
cue.full-cue.notrt
cue.full-cue.nodate

#- so the warming by date interaction is "significant" and leads to a better model, but it is quantitatively of modest 
#-  importance, as the marginal r2 drops just a little when excluding the interaction.

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------






#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#- merge models together to make Table 1
table.r1 <- as.matrix(anova(sp.ra.ar1.reml,type="marginal"))
table.r2 <- cbind(table.r1,as.matrix(anova(sp.gpp.ar1.reml,type="marginal"))[1:4,1:4])
table1 <- cbind(table.r2,as.matrix(anova(sp.cue.ar1.reml,type="marginal"))[1:4,1:4])[2:4,]
row.names(table1) <- c("Warming","Date","Warming x Date")
table1 

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------







#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#-- analysis of just the exceedingly hot days
hotDates <- unique(dat.hr.p[which(dat.hr.p$Tair_al>40),"Date"]) # find dates with temperatures exceeding 40 deg C
dat.hr.p.hot <- subset(dat.hr.p,Date %in% hotDates)

#- daily climate metrics
hotDates_met <- summaryBy(Tair_al~T_treatment+Date,data=dat.hr.p.hot,FUN=c(mean,min,max),na.rm=T)
summaryBy(Tair_al.max~T_treatment,data=hotDates_met) # average maximum temperature on these hot dates

#- re-analyze on hot dates only. Capture heteroskedasticity across dates
sp.CUE.hot <- lme(RtoA~T_treatment*DateFac,random=list(~1|chamber),
                  weights=varIdent(form=~1|DateFac),data=subset(dat2,Date %in% hotDates))

#look at model diagnostics
plot(sp.CUE.hot,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.CUE.hot,RtoA~fitted(.)|chamber,abline=c(0,1))         #predicted vs. fitted for each species
plot(sp.CUE.hot,RtoA~fitted(.),abline=c(0,1))              #predicted vs. fitted
qqnorm(sp.CUE.hot, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot
anova(sp.CUE.hot)
anova(sp.CUE.hot,type="marginal")
lsmeans(sp.CUE.hot,"T_treatment") 
CIs <- intervals(sp.CUE.hot,which="fixed")
CIs$fixed[2,] # 95% confidence intervals for T_treatment effect

#-----
#- alternatively, just analyze data above the inflection point between Ra/GPP and T (24-hr average T > 22),
hotDates2 <- unique(cue.day[which(cue.day$Tair_24hrs>22 & cue.day$T_treatment=="ambient"),"Date"]) # find dates with temperatures exceeding 40 deg C

# summarize these hot dates
dat.hr.p.hot2 <- subset(dat.hr.p,Date %in% hotDates2)
hotDates_met2 <- summaryBy(Tair_al~T_treatment+Date,data=dat.hr.p.hot2,FUN=c(mean,min,max),na.rm=T)
summaryBy(Tair_al.mean+Tair_al.min+Tair_al.max~T_treatment,data=hotDates_met2) # average maximum temperature on these hot dates


#- re-analyze on hot dates only. Capture heteroskedasticity across dates
hot2 <- subset(dat2,Date %in% hotDates2)
hot2$DateFac <- factor(hot2$DateFac)
sp.CUE.hot2 <- lme(RtoA~T_treatment*DateFac,random=list(~1|chamber),
                  weights=varIdent(form=~1|DateFac),
                  data=hot2)


#look at model diagnostics
plot(sp.CUE.hot2,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.CUE.hot2,RtoA~fitted(.)|chamber,abline=c(0,1))         #predicted vs. fitted for each chamber
plot(sp.CUE.hot2,RtoA~fitted(.),abline=c(0,1))              #predicted vs. fitted
qqnorm(sp.CUE.hot2, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot. not bad
anova(sp.CUE.hot2)
lsmeans(sp.CUE.hot2,"T_treatment") 
CIs2 <- intervals(sp.CUE.hot2,which="fixed")
CIs2$fixed[2,] # 95% confidence intervals for T_treatment effect






#-----
#- alternatively, just analyze all data where maximum air T exceeded 33 in warmed.
#- note that most of tehse are the same dates where mean 24-hour temperature exceeded 22 deg C
hotDates3 <- unique(dat.hr.p[which(dat.hr.p$Tair_al>33),"Date"]) # find dates with temperatures exceeding 33 deg C # find dates with temperatures exceeding 40 deg C

# summarize these hot dates
dat.hr.p.hot3 <- subset(dat.hr.p,Date %in% hotDates3)
hotDates_met3 <- summaryBy(Tair_al~T_treatment+Date,data=dat.hr.p.hot3,FUN=c(mean,min,max),na.rm=T)
summaryBy(Tair_al.mean+Tair_al.min+Tair_al.max~T_treatment,data=hotDates_met3) # average maximum temperature on these hot dates


#- re-analyze on hot dates only. Capture heteroskedasticity across dates
hot3 <- subset(dat2,Date %in% hotDates3)
hot3$DateFac <- factor(hot3$DateFac)
sp.CUE.hot3 <- lme(RtoA~T_treatment*DateFac,random=list(~1|chamber),
                   weights=varIdent(form=~1|DateFac),
                   data=hot3)

#look at model diagnostics
plot(sp.CUE.hot3,resid(.,type="p")~fitted(.) | T_treatment,abline=0)   #resid vs. fitted for each treatment
plot(sp.CUE.hot3,RtoA~fitted(.)|chamber,abline=c(0,1))         #predicted vs. fitted for each chamber
plot(sp.CUE.hot3,RtoA~fitted(.),abline=c(0,1))              #predicted vs. fitted
qqnorm(sp.CUE.hot3, ~ resid(., type = "p"), abline = c(0, 1))     #qqplot. not bad
anova(sp.CUE.hot3)
lsmeans(sp.CUE.hot3,"T_treatment") 
CIs3 <- intervals(sp.CUE.hot3,which="fixed")
CIs3$fixed[2,] # 95% confidence intervals for T_treatment effect
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#- back of the envelope calculation for one of the reviewers regarding their comment "(5)	In a back of the envelope 
#      calculation – how much would the “hot day” increase in Ra/GPP affect overall carbon balance ? "
#-----
#- alternatively, just analyze data above the inflection point between Ra/GPP and T (24-hr average T > 22),
hotDates2 <- unique(cue.day[which(cue.day$Tair_24hrs>22 & cue.day$T_treatment=="ambient"),"Date"]) # find dates with temperatures exceeding 40 deg C
dat2$hotorcold <- ifelse(dat2$Date %in% hotDates2,"hot","cold")

dat3 <- summaryBy(Ra_la+GPP_la+RtoA~Date+T_treatment+hotorcold,data=dat2,FUN=mean,keep.names=T)
summaryBy(Ra_la+GPP_la+RtoA~hotorcold,FUN=sum,data=dat3)
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

