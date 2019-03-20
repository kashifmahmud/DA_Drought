# Script to process leaf storage data 
# read the leaf tnc data over time
tnc = read.csv("raw_data/WTC_TEMP_CM_LEAFCARB_20130515-20140402_L2.csv") # unit = mg of tnc per g of dry leafmass
tnc = na.omit(tnc)
tnc$Date = as.Date(tnc$Date)
tnc = subset(tnc, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-05-27"))
# tnc$chamber_type = as.factor( ifelse(tnc$chamber %in% drought.chamb, "drought", "watered") )
keeps <- c("Date", "chamber", "T_treatment", "TNC_mgg")
tnc = tnc[ , keeps, drop = FALSE]
tnc.mean <- summaryBy(TNC_mgg ~ Date+T_treatment, data=tnc, FUN=c(mean,standard.error))
names(tnc.mean)[3:4] = c("tnc.conc","tnc.conc_SE")

# read the diurnal leaf tnc data - 2-hourly intervals during a sunny (20 February) and an overcast/rainy day (26 March)
tnc.diurnal = read.csv("raw_data/WTC_TEMP_CM_LEAFCARB-DIURNAL_20140220-20140326_R.csv") # unit = mg of tnc per g of dry leafmass
tnc.diurnal$Sample.Date <- parse_date_time(tnc.diurnal$Sample.Date,"d m y")
tnc.diurnal$Sample.Date = as.Date(tnc.diurnal$Sample.Date, format = "%d/%m/%Y")
keeps <- c("Sample.Date", "Chamber", "TNC")
tnc.diurnal = tnc.diurnal[ , keeps, drop = FALSE]
names(tnc.diurnal) = c("Date","chamber","TNC")
tnc.diurnal = merge(tnc.diurnal, unique(height.dia[,c("chamber","T_treatment")]), by="chamber")
tnc.diurnal.mean <- summaryBy(TNC ~ Date+T_treatment, data=tnc.diurnal, FUN=c(mean,standard.error))
names(tnc.diurnal.mean)[3:4] = c("tnc.conc","tnc.conc_SE")

# read the diurnal leaf tnc data - trees were girdled
tnc.girdled = read.csv("raw_data/WTC_TEMP_CM_PETIOLEGIRDLE-LEAFMASS-AREA-CARB_20140507-20140512_L2.csv") # unit = mg of tnc per g of dry leafmass
tnc.girdled$Date <- parse_date_time(tnc.girdled$Date,"y m d")
tnc.girdled$Date = as.Date(tnc.girdled$Date, format = "%Y/%m/%d")
tnc.girdled = subset(tnc.girdled, treatment %in% as.factor("control"))
tnc.girdled$TNC = tnc.girdled$starch + tnc.girdled$soluble_sugars
tnc.girdled = merge(tnc.girdled, unique(height.dia[,c("chamber","T_treatment")]), by="chamber")
tnc.girdled.mean <- summaryBy(TNC ~ Date+T_treatment, data=tnc.girdled, FUN=c(mean,standard.error))
names(tnc.girdled.mean)[3:4] = c("tnc.conc","tnc.conc_SE")

# bind all available tnc data into one dataframe
tnc.mean = rbind(tnc.mean, tnc.diurnal.mean)
tnc.mean = rbind(tnc.mean, tnc.girdled.mean)

# unit conversion: 1 g tnc has 0.4 gC (12/30) and 1 g plant has 0.48 gC
tnc.mean$tnc.conc = tnc.mean$tnc.conc / 1000 # g of tnc per g of dry leafmass
tnc.mean$tnc.conc_SE = tnc.mean$tnc.conc_SE / 1000 # g of tnc per g of dry leafmass
tnc.mean$tnc.conc = tnc.mean$tnc.conc * 0.4 / c1 # gC in tnc per gC of dry leafmass
tnc.mean$tnc.conc_SE = tnc.mean$tnc.conc_SE * 0.4 / c1 # gC in tnc per gC of dry leafmass

# get the daily interpotalerd LM
# treeMass.daily$chamber_type = as.factor( ifelse(treeMass.daily$chamber %in% drought.chamb, "drought", "watered") )
treeMass.daily.sum = summaryBy(leafMass ~ Date+T_treatment, data=treeMass.daily, FUN=c(mean,standard.error))
names(treeMass.daily.sum)[3:4] = c("LM","LM_SE")
treeMass.daily.sum[,c(3:4)] = treeMass.daily.sum[,c(3:4)] * c1 # unit conversion from gDM to gC

tnc.final = merge(tnc.mean, treeMass.daily.sum[,c("Date","T_treatment","LM","LM_SE")], by=c("Date","T_treatment"), all=FALSE)
tnc.final$TNC_leaf = tnc.final$tnc.conc * tnc.final$LM # Unit = gC
# tnc.final$TNC_leaf_SE = tnc.final$tnc.conc_SE * tnc.final$LM # Unit = gC

tnc.final$TNC_leaf_SE = (((tnc.final$tnc.conc_SE/tnc.final$tnc.conc)^2 + (tnc.final$LM_SE/tnc.final$LM)^2)^0.5) * tnc.final$TNC_leaf


tnc.final = merge(tnc.final,tnc.partitioning, by="Date")

# Estimate total plant TNC from the partitioning data of WTC-4
tnc.final$TNC_tot = tnc.final$TNC_leaf / tnc.final$foliage # Unit = gC
tnc.final$TNC_tot_SE = tnc.final$TNC_leaf_SE / tnc.final$foliage # Unit = gC

tnc.final$TNC_wood = tnc.final$TNC_tot * tnc.final$wood # Unit = gC
tnc.final$TNC_wood_SE = tnc.final$TNC_tot_SE * tnc.final$wood # Unit = gC

tnc.final$TNC_root = tnc.final$TNC_tot * tnc.final$root # Unit = gC
tnc.final$TNC_root_SE = tnc.final$TNC_tot_SE * tnc.final$root # Unit = gC



# tnc.final$TNC_tot = tnc.final$TNC_leaf / tnc.partitioning$average_ratio[tnc.partitioning$organs == "foliage"] # Unit = gC
# tnc.final$TNC_tot_SE = tnc.final$TNC_leaf_SE / tnc.partitioning$average_ratio[tnc.partitioning$organs == "foliage"] # Unit = gC
# 
# tnc.final$TNC_wood = tnc.final$TNC_tot * tnc.partitioning$average_ratio[tnc.partitioning$organs == "wood"] # Unit = gC
# tnc.final$TNC_wood_SE = tnc.final$TNC_tot_SE * tnc.partitioning$average_ratio[tnc.partitioning$organs == "wood"] # Unit = gC
# 
# tnc.final$TNC_root = tnc.final$TNC_tot * tnc.partitioning$average_ratio[tnc.partitioning$organs == "root"] # Unit = gC
# tnc.final$TNC_root_SE = tnc.final$TNC_tot_SE * tnc.partitioning$average_ratio[tnc.partitioning$organs == "root"] # Unit = gC




