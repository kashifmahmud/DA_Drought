#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- Central analysis script for DA with WTC3 experiment manuscript, entitled
#  "".

#  The idea is to keep this script nice and tidy, but reproducibly do all the
#  analysis and make all of the figures for the manuscript. Raw and processed data will be 
#  placed in the "raw_data" and "processed_data" folders respectively, while figures 
#  and tables will be placed in the "output" folder.
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# # Clear the workspace (if needed)
# rm(list=ls())
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# Load required libraries. There are quite a few, including some non-standard functions that are not on CRAN. 
# This script will check for required libraries and install any that are missing	
source('R/load_packages_drought.R')	

# Load the custom analysis and plotting functions that do all of the actual work	
source("R/functions_drought.R")	
source("R/functions_drought_CBM.R")	
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# #- Download data files for WTC3 experiment. This downloads the zipfile from figshare
# download.file("https://ndownloader.figshare.com/files/4857112", "data.zip", mode="wb")
# # Extract data to different folders.
# unzip("data.zip")

#-------------------------------------------------------------------------------------

# #-------------------------------------------------------------------------------------
#- This script imports and processes the raw WTC3 experiment data to model the carbon pools and fluxes using DA
# source("R/initial_data_processing_wtc3.R")
rmd2rscript("report_initial_data_processing_drought.Rmd")
source("report_initial_data_processing_drought.R")
# #-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# #- Make figure 1. Model representation of storage, allocation and autotrophic respiration processes and 
# # pathways in the CBM with storage pool, separate growth and maintenance respiration components.
# plot.model.wtc3()
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Read the processed data and clean the workspace from data pre-processing
#  Processed data are placed in "processed_data" folder
# source("R/read_data_wtc3.R")
#-------------------------------------------------------------------------------------
# Model run for WTC3 dataset (without LA feedback)
# treat.group = as.factor(c("ambient")) # Assign 1 treatment to check the model
# treat.group = as.factor(c("ambient drought","elevated drought")) # Assign 1 treatment to check the model
# treat.group = as.factor(c("ambient drought","ambient watered")) # Assign 2 treatments to compare the results

treat.group = as.factor(c("ambient","elevated")) # Assign all treatments
data.all = read.csv("processed_data/data_all.csv") 
data.all$treatment.no [data.all$Treatment %in% as.factor("ambient") ] = 1
data.all$treatment.no [data.all$Treatment %in% as.factor("elevated") ] = 2
tnc.partitioning = read.csv("processed_data/tnc_partitioning_data.csv")
# data.all[,c("GPP_SE","Ra_SE","LA_SE","LM_SE","WM_SE","RM_SE","litter_SE","TNC_tot_SE","TNC_leaf_SE","TNC_wood_SE","TNC_root_SE")] =
#   0.5*data.all[,c("GPP_SE","Ra_SE","LA_SE","LM_SE","WM_SE","RM_SE","litter_SE","TNC_tot_SE","TNC_leaf_SE","TNC_wood_SE","TNC_root_SE")]
# data.all[,c("TNC_tot_SE","TNC_leaf_SE","TNC_wood_SE","TNC_root_SE")] = 3 * data.all[,c("TNC_tot_SE","TNC_leaf_SE","TNC_wood_SE","TNC_root_SE")]

# #-------------------------------------------------------------------------------------
# #- Matching C balance of the entire experiment considering C inputs and outputs
# source("R/C_balance_wtc3.R")
# 
# #-------------------------------------------------------------------------------------
# # # 3000 chain length is sufficient for the convergance
chainLength = 100
no.param.par.var = 2
with.storage = T
model.comparison=F
model.optimization=F
treat.group=c(list(1,2))
# treat.group=c(list(list(c(1,2),c(1,2))))
# # start <- proc.time() # Start clock
# # # result = CBM.wtc3(chainLength = 3000, no.param.par.var=(nrow(data.all)/4)/30, treat.group=treat.group, with.storage, model.comparison=F, model.optimization=F) # Monthly parameters
# result = CBM.wtc3(chainLength, no.param.par.var, treat.group, with.storage, model.comparison, model.optimization) # Quadratic/Cubic parameters
# 
# # Run the model with constant k and Y parameters
# result = CBM.wtc3_const_k_Y(chainLength, no.param.par.var, treat.group, with.storage, model.comparison, model.optimization) # Quadratic/Cubic parameters
# 
# # time_elapsed_series <- proc.time() - start # End clock
# # result[[6]]
# # write.csv(result[[6]], "output/bic.csv", row.names=FALSE) # unit of respiration rates: gC per gC plant per day
# # 
# # Plot parameters and biomass data fit
# plot.Modelled.parameters(result,with.storage)
# plot.Modelled.biomass(result,with.storage)
#-------------------------------------------------------------------------------------
source("R/functions_wtc3.R")	
source("R/functions_wtc3_CBM.R")	

# Model run for WTC3 dataset with clustering
cluster <- makeCluster(detectCores()-1)
# clusterEvalQ(cluster, library(xts))
clusterExport(cl=cluster, list("data.all","tnc.partitioning","treat.group"))
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterExport(cluster, ex)
result.cluster = list()
bic.cluster = list()

start <- proc.time() # Start clock
result <- clusterMap(cluster, CBM.wtc3, treat.group=c(list(list(1,2))),
                     MoreArgs=list(chainLength=1000, no.param.par.var=2, with.storage=T, model.comparison=F, model.optimization=F))

# # Test whether parameters need to be seperate for both ambient and warmed treatments
# result <- clusterMap(cluster, CBM.wtc3, treat.group=c(list(list(1,2,c(1,2)))),
#                      MoreArgs=list(chainLength=300, no.param.par.var=2, with.storage=T, model.comparison=F, model.optimization=F))

time_elapsed_series <- proc.time() - start # End clock
stopCluster(cluster)

# listOfDataFrames <- vector(mode = "list", length = length(treat.group[[1]]))
# for (i in 1:length(treat.group[[1]])) {
#   listOfDataFrames[[i]] <- data.frame(result[[i]][[6]])
# }
# bic = do.call("rbind", listOfDataFrames)
# write.csv(bic, "output/bic.csv", row.names=FALSE)
write.csv(result[[1]][[6]], "output/bic.csv", row.names=FALSE)

# Plot parameters and biomass data fit
plot.Modelled.parameters.wtc3(result,with.storage=T)
plot.Modelled.biomass.wtc3(result,with.storage=T)
result[[1]][[6]]
#-------------------------------------------------------------------------------------
# Calculate total C partitioning for individual treatments 
# and make figure 7 and Table S1
# source("R/functions_CBM.R")
source("R/C_partitioning_wtc3.R")
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
# Check Rabove from data
source("R/Rabove_balance_wtc3.R")

#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
# Check Tree height whether they hit the chamber top
source("R/check_tree_height_wtc3.R")

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# Check initial and final biomass variation
source("R/check_biomass_wtc3.R")

#-------------------------------------------------------------------------------------

# Check the C input and output balance
data.set = subset(data.all,(Treatment %in% treat.group[v]))
input = sum(data.set$GPP)
Rm = subset(result[[4]], variable %in% as.factor("Rm"))
Rm.sum = sum(Rm$value)

Rm.daily = mean(data.set$Rd.foliage.mean)*(data.set$LM[1]+data.set$LM[nrow(data.set)])/2 + mean(data.set$Rd.stem.mean)*(data.set$WM[1]+data.set$WM[nrow(data.set)])/2*mean(data.set$SMratio) + 
  mean(data.set$Rd.branch.mean)*(data.set$WM[1]+data.set$WM[nrow(data.set)])/2*mean(data.set$BMratio) + 
  mean(data.set$Rd.fineroot.mean)*(data.set$RM[1]+data.set$RM[nrow(data.set)])/2*mean(data.set$FRratio_SE) + mean(data.set$Rd.intermediateroot.mean)*(data.set$RM[1]+data.set$RM[nrow(data.set)])/2*mean(data.set$IRratio) + 
  mean(data.set$Rd.coarseroot.mean)*(data.set$RM[1]+data.set$RM[nrow(data.set)])/2*mean(data.set$CRratio_SE) + mean(data.set$Rd.boleroot.mean)*(data.set$RM[1]+data.set$RM[nrow(data.set)])/2*mean(data.set$BRratio) 
Rm.sum = Rm.daily * 252

output = Rm.sum + (data.set$LM[nrow(data.set)] - data.set$LM[1]) + (data.set$WM[nrow(data.set)] - data.set$WM[1]) + (data.set$RM[nrow(data.set)] - data.set$RM[1]) + 
  (data.set$TNC_tot[max(which(complete.cases(data.set$TNC_tot)))] - data.set$TNC_tot[min(which(complete.cases(data.set$TNC_tot)))]) + data.set$litter[max(which(complete.cases(data.set$litter)))]

Y.modelled = subset(result[[2]], variable %in% as.factor("Y"))

#-------------------------------------------------------------------------------------
# Check Rabove from data
Rabove.ini = data.set$Rd.foliage.mean[1]*data.set$LM[1] + data.set$Rd.stem.mean[1]*data.set$WM[1]*data.set$SMratio[1] + data.set$Rd.branch.mean[1]*data.set$WM[1]*data.set$BMratio[1] +
  Y.modelled$Parameter[1]*((data.set$LM[nrow(data.set)]-data.set$LM[1])/251 + (data.set$WM[nrow(data.set)]-data.set$WM[1])/251)

Rabove.end = data.set$Rd.foliage.mean[nrow(data.set)]*data.set$LM[nrow(data.set)] + data.set$Rd.stem.mean[nrow(data.set)]*data.set$WM[nrow(data.set)]*data.set$SMratio[nrow(data.set)] + data.set$Rd.branch.mean[nrow(data.set)]*data.set$WM[nrow(data.set)]*data.set$BMratio[nrow(data.set)] +
  Y.modelled$Parameter[nrow(data.set)]*((data.set$LM[nrow(data.set)]-data.set$LM[1])/251 + (data.set$WM[nrow(data.set)]-data.set$WM[1])/251)

# mean(data.set$Ra)/mean(data.set$GPP)

# Rabove = mean(data.set$Rd.foliage.mean)*(data.set$LM[1]+data.set$LM[nrow(data.set)])/2 + mean(data.set$Rd.stem.mean)*(data.set$WM[1]+data.set$WM[nrow(data.set)])/2*mean(data.set$SMratio)
# + mean(data.set$Rd.branch.mean)*(data.set$WM[1]+data.set$WM[nrow(data.set)])/2*mean(data.set$BMratio) +
#   mean(Y.modelled$Parameter)*(data.set$LM[nrow(data.set)]-data.set$LM[1] + data.set$WM[nrow(data.set)]-data.set$WM[1])
# Ra = sum(data.set$Ra)

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------

# cluster <- makeCluster(detectCores()-1)
# # clusterEvalQ(cluster, library(xts))
# clusterExport(cl=cluster, list("treat.group","data.all","tnc.partitioning"))
# ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
# clusterExport(cluster, ex)
# result.cluster = list()
# bic.cluster = list()
# 
# start <- proc.time() # Start clock
# result.cluster <- clusterMap(cluster, CBM.wtc3, no.param.par.var=c((nrow(data.all)/4)/30,(nrow(data.all)/4)/7), 
#             MoreArgs=list(chainLength=100, treat.group=treat.group, with.storage=c(T,T,T,T), model.comparison=c(F,F,F,F), model.optimization=c(F,F,F,F)))
# 
# time_elapsed_series <- proc.time() - start # End clock
# bic.without.storage = result.cluster[[1]][[6]]
# bic.with.storage = result.cluster[[2]][[6]]
# bic.group1 = result.cluster[[3]]
# bic.group2 = result.cluster[[4]]
# bic.group3 = result.cluster[[5]]
# bic.group4 = result.cluster[[6]]
# result = result.cluster[[7]]
# stopCluster(cluster)

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------



