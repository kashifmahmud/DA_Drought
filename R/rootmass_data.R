#- Get the root mass (initial - estimated and final - harvested) for each of the treatment

# Script to setup a model for the root mass from sink limited pot experiment and WTC3 harvest data to predict 
# initial root mass from Height and Diameter
plot.summary.pot = read.csv("data_pot_experiment/plot_summary.csv")
initial.data.pot <- read.csv("data_pot_experiment/seedling_initial.csv")
keeps <- c("root_mass", "diameter_15", "height")
data.pot = initial.data.pot[ , keeps, drop = FALSE]
names(data.pot) <- c("rootmass","diameter","height")
data.pot$volume = as.integer(1000)

harvest.data.pot <- read.csv("data_pot_experiment/seedling_mass.csv")
harvest.data.pot$rootmass = harvest.data.pot$coarseroot + harvest.data.pot$fineroot
keeps <- c("rootmass", "plot", "pot", "volume")
harvest.data.pot = harvest.data.pot[ , keeps, drop = FALSE]

height.dia.pot <- read.csv("data_pot_experiment/height_diameter.csv")
height.dia.pot$Date <- parse_date_time(height.dia.pot$Date,"d m y")
height.dia.pot$Date = as.Date(height.dia.pot$Date, format = "%d/%m/%Y")
height.dia.pot = subset(height.dia.pot,Date=="2013-05-21")

height.dia.pot <- merge(plot.summary.pot, height.dia.pot, by=c("pot","plot"))
harvest.data.pot <- merge(height.dia.pot, harvest.data.pot, by=c("pot","plot","volume"))
keeps <- c("rootmass", "diameter", "height","volume")
harvest.data.pot = harvest.data.pot[ , keeps, drop = FALSE]
data.pot = rbind(data.pot, harvest.data.pot)
# data.pot[,"rootmass"] = data.pot[,"rootmass"] * c1 # unit conversion: gDM to gC 
data.free.pot = subset(data.pot, volume %in% 1000)
data.free.pot = data.free.pot[,-4]

# Read harvest rootmass data form WTC3 experiment
rootmass.harvest = read.csv("raw_data/WTC_TEMP_CM_HARVEST-ROOTS_20140529-20140606_L1_v1.csv")
# rootmass.harvest$chamber_type = as.factor( ifelse(rootmass.harvest$chamber %in% drought.chamb, "drought", "watered") )
rootmass.harvest$Date = as.Date("2014-05-26")
# rootmass.harvest = merge(unique(height.dia[,c("chamber","T_treatment","chamber_type")]), rootmass.harvest, by=c("chamber"))
rootmass.harvest = merge(unique(height.dia[,c("chamber","T_treatment")]), rootmass.harvest, by=c("chamber"))

# height.dia.sub = subset(height.dia, Date %in% as.Date("2014-05-27") & T_treatment %in% as.factor("ambient") & chamber_type %in% as.factor("watered"))
height.dia.sub = subset(height.dia, Date %in% as.Date("2014-05-27") & T_treatment %in% as.factor("ambient"))
keeps <- c("chamber", "height", "diameter")
height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
# height.dia.sub <- merge(height.dia.sub, rootmass.harvest[,c("chamber","RootDMtotal","chamber_type")], by=c("chamber"))
height.dia.sub <- merge(height.dia.sub, rootmass.harvest[,c("chamber","RootDMtotal")], by=c("chamber"))
keeps <- c("RootDMtotal", "diameter", "height")
height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
names(height.dia.sub) <- c("rootmass","diameter","height")

data.free.pot = rbind(data.free.pot, height.dia.sub) # Merge pot exp data (initial and harvest of free seedlings) and WTC3 harvest data (ambient, watered)

# Fit a linear regression to estimate root mass by stem dia and height (ignoring temperature variation)
# the regression is fitted using the Sink limited pot experiemnt of Court Campany for similar seedlings (Euc.Teri.)
rm1 <- lm(log(rootmass) ~ log(diameter) + log(height), data=data.free.pot)
rm2 <- lm(log(rootmass) ~ log(diameter) * log(height), data=data.free.pot)

# Estimate the WTC3 root mass from the fitted linear regression equation
# Go for the regression with interaction (rm1) that suits the pot experiment (Euc.Teri.) data best
# otherwise we could only use diameter to predict the rootmass (as intercept and height aren't significant)
# eq.rm = function(x,y){exp(coefficients(rm1)[1] + coefficients(rm1)[2] * log(x)  + coefficients(rm1)[3] * log(y))}
eq.rm = function(x,y){exp(coefficients(rm2)[1] + coefficients(rm2)[2] * log(x)  + coefficients(rm2)[3] * log(y)
                          + coefficients(rm2)[4] * log(x) * log(y) )}

# estimate the initial rootmass from H and D data
# height.dia.initial = subset(height.dia.final, Date %in% as.Date(c("2013-09-17","2014-02-04")))
height.dia.initial = subset(height.dia.final, Date %in% as.Date("2013-09-17"))
height.dia.initial$RM = eq.rm(height.dia.initial$diameter, height.dia.initial$height)
# height.dia.initial$RM_SE = (((height.dia.initial$diameter_SE/height.dia.initial$diameter)^2 + (height.dia.initial$height_SE/height.dia.initial$height)^2 )^0.5) * height.dia.initial$RM
height.dia.initial$RM_SE = ( (((coefficients(rm2)[2]*height.dia.initial$diameter_SE/height.dia.initial$diameter)^2 + (coefficients(rm2)[3]*height.dia.initial$height_SE/height.dia.initial$height)^2 +
                                 coefficients(rm2)[4]*((height.dia.initial$diameter_SE/height.dia.initial$diameter)^2 + (height.dia.initial$height_SE/height.dia.initial$height)^2) )^0.5) ) * height.dia.initial$RM


# processing the harvest rootmass
rootmass.harvest.mean <- summaryBy(RootDMtotal ~ Date+T_treatment+chamber_type, data=rootmass.harvest, FUN=c(mean,standard.error))
names(rootmass.harvest.mean)[3:4] = c("RM","RM_SE")

rootmass = merge(height.dia.initial[,c("Date","T_treatment","RM","RM_SE")], rootmass.harvest.mean, all = TRUE)
# rootmass = height.dia.initial[,c("Date","T_treatment","chamber_type","RM","RM_SE")]
# rootmass = height.dia.initial[,c("Date","T_treatment","RM","RM_SE")]
rootmass$RM = rootmass$RM * c1 # unit conversion from gDM to gC
rootmass$RM_SE = rootmass$RM_SE * c1 # unit conversion from gDM to gC

