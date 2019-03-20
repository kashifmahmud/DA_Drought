#- This sript reads the processed data from Sink limited container volume experiment 
#-----------------------------------------------------------------------------------------
# Import daily GPP, daily Rd data
data.processed = read.csv("processed_data/data_GPP_Ra_LA.csv") 
GPP.data.processed = data.processed[,c("Date","T_treatment","chamber_type","GPP")] # Units gC d-1

Ra.data.processed = data.processed[,c("Date","T_treatment","chamber_type","Ra","Ra_SE")] # Unit gC per gC plant per day

data.biomass.processed = read.csv("processed_data/data_biomass_litter_tnc.csv")
tnc.data.processed = data.biomass.processed[,c("Date","T_treatment","chamber_type","TNC","TNC_SE")] # Unit gC per gC plant

# Import weekly Cleaf, weekly Cstem, initial/harvest Croot data with Mean and SD
Mleaf.data.processed = data.biomass.processed[,c("Date","T_treatment","chamber_type","LM","LM_SE")] # Unit gC
Mstem.data.processed = data.biomass.processed[,c("Date","T_treatment","chamber_type","WM","WM_SE")] # Unit gC
Mroot.data.processed = data.biomass.processed[,c("Date","T_treatment","chamber_type","RM","RM_SE")] # Unit gC
LA.data.processed = data.processed[,c("Date","T_treatment","chamber_type","LA","LA_SE")] # Unit m^2

# # Import harvest data
# sla.harvest.processed = read.csv("processed_data/sla.harvest.csv") # Unit of SLA = m2 leaf area per gC of leaf biomass
# 
# # Import the self shading factors
# # sigma.data.processed <- read.csv("processed_data/M_leafarea_model.csv") # Unitless
# sigma.data.processed <- read.csv("processed_data/M_leafarea_model_prev.csv") # Unitless

