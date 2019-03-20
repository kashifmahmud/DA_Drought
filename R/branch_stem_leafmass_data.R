####### From John's script #########
#- Get an estimate of branch, stem, and leaf mass as well as leaf area for each day of the experiment
treeMass.daily <- returnTreeMass()
# treeMassFlux <- merge(dat.hr.gf,treeMass,by=c("chamber","Date","T_treatment","Water_treatment"))
treeMass.daily = subset(treeMass.daily, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-05-26"))

treeMass.daily$Measurement[treeMass.daily$Date %in% as.Date("2014-05-26")] = 40
treeMass = na.omit(treeMass.daily) # consider only fortnightly direct measurements of H and D
# treeMass = treeMass.daily # consider daily data interpolated from fortnightly direct measurements of H and D

treeMass$woodMass = treeMass$branchMass + treeMass$boleMass
# treeMass$chamber_type = as.factor( ifelse(treeMass$chamber %in% drought.chamb, "drought", "watered") )
# treeMass.sum = summaryBy(leafMass+woodMass ~ Date+T_treatment+chamber_type, data=treeMass, FUN=c(mean,standard.error))
# names(treeMass.sum)[4:7] = c("LM","WM","LM_SE","WM_SE")
# treeMass.sum[,c(4:7)] = treeMass.sum[,c(4:7)] * c1 # unit conversion from gDM to gC
treeMass.sum = summaryBy(leafMass+woodMass ~ Date+T_treatment, data=treeMass, FUN=c(mean,standard.error))
names(treeMass.sum)[3:6] = c("LM","WM","LM_SE","WM_SE")
treeMass.sum[,c(3:6)] = treeMass.sum[,c(3:6)] * c1 # unit conversion from gDM to gC

