# This script take the subset of processed data for particular treatment group
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
# Subset the data for different treatments (e.g. 5L, 10L, ....., 1000L)
# Consider one treatment group at a time to run MCMC with the CBM
GPP.data = subset(GPP.data.processed,(volume %in% vol[v])) 
names(GPP.data)[3] = "GPP"
Rd.data = subset(Rd.data.processed,(volume %in% vol[v]))
Mleaf.data = subset(Mleaf.data.processed,(volume %in% vol[v]))
Mstem.data = subset(Mstem.data.processed,(volume %in% vol[v]))
Mroot.data = subset(Mroot.data.processed,(volume %in% vol[v]))
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# Merge all GPP, Rd, Sleaf, Cleaf, Cstem, Croot data
data = merge(GPP.data,Rd.data, all = TRUE)
data = merge(data,Mleaf.data, all = TRUE)
data = merge(data,Mstem.data, all = TRUE)
data = merge(data,Mroot.data, all = TRUE)
names(data)[4:ncol(data)] = c("Rd","Mleaf","Mleaf_SD","Mstem","Mstem_SD","Mroot","Mroot_SD")
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
if (with.storage==T) {
  Sleaf.data = tnc.data = subset(tnc.data.processed,(volume %in% vol[v]))
  keeps <- c("Date", "volume", "tnc", "tnc_SE")
  Sleaf.data = Sleaf.data[ , keeps, drop = FALSE]
  
  data = merge(data,Sleaf.data, all = TRUE)
  names(data)[4:ncol(data)] = c("Rd","Mleaf","Mleaf_SD","Mstem","Mstem_SD","Mroot","Mroot_SD","Sleaf","Sleaf_SD")
# data = data[with(data, order(volume)), ]
# row.names(data) = c(1:nrow(data))
}

# write.csv(lm.daily.melt[,c("Date","volume","leafmass","leafmass_SE")], file = "processed_data/Cleaf_daily_data.csv", row.names = FALSE)
#-------------------------------------------------------------------------------------
