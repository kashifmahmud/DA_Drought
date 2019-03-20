# run coupled photo-gs model, will retrun Aleaf at 15min interval

# # convert RH to VPD
# Asat$VPD <- RHtoVPD(Asat$RH, Asat$temp, Pa=101)


# run coupled photo-gs model, assume Rd = 0
A_pred <- Photosyn(VPD=Asat$VPD, Ca=400, PPFD=Asat$PPFD, Tleaf=Asat$temp, 
                   Jmax=Asat$Jmax.mean, Vcmax=Asat$Vcmax.mean, Rd=0, g1=Asat$g1)


# need a new dfr with Aleaf and Anet across the day
Aleaf <- A_pred[,c(1:4, 7:11)]
Aleaf_15min <- cbind(Aleaf, Asat[,c(1, 5:6)])
write.csv(Aleaf_15min, "calculated data/Aleaf_15min_gross.csv", row.names=FALSE)


Aleaf_15min$Date <- as.Date(Aleaf_15min$Date)
Aleaf_15min$volume <- as.factor(Aleaf_15min$volume)
Aleaf_15min$photo15gc <- with(Aleaf_15min, ALEAF*15*60*10^-6*12)

Aleaf <- summaryBy(photo15gc ~ Date+Water, data=Aleaf_15min, FUN=sum, keep.names=TRUE )
names(Aleaf)[3] <- "carbon_day"

write.csv(Aleaf, "processed_data/Aleaf_model/cday_120_clean_gross.csv", row.names=FALSE)
Aleaf_agg <- summaryBy(carbon_day ~ Water, data=Aleaf, FUN=mean, keep.names=TRUE )
write.csv(Aleaf_agg, "processed_data/Aleaf_model/gCday_means_clean_gross.csv", row.names=FALSE)


