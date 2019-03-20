# Script to read and plot stem height and diameter for various treatment cases 
# Read data form WTC3 experiment
height.dia.raw <- read.csv("raw_data/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV") # units: Height = cm, dia = mm
height.dia.raw$DateTime = as.Date(height.dia.raw$DateTime)
flux <- read.csv("raw_data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
chambers = unique(flux$chamber)

height.dia = data.frame(chamber = rep(chambers, each = length(unique(height.dia.raw$DateTime))),
                        Date = rep(unique(height.dia.raw$DateTime),length(chambers)),
                        T_treatment = character(length(chambers) * length(unique(height.dia.raw$DateTime))),
                        # W_treatment = character(length(chambers) * length(unique(height.dia.raw$DateTime))),
                        diameter = numeric(length(chambers) * length(unique(height.dia.raw$DateTime))),
                        height = numeric(length(chambers) * length(unique(height.dia.raw$DateTime))), stringsAsFactors=FALSE)
height.dia = subset(height.dia, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-05-27"))
# height.dia = subset(height.dia, Date >= as.Date("2013-09-17") & Date <= as.Date("2014-02-12"))

for(i in 1:length(chambers)) {
  height.dia.sub = subset(height.dia.raw, chamber %in% as.factor(chambers[i]))
  # height.dia.sub = subset(height.dia, chamber %in% as.factor(chambers[i]) & Water_treatment %in% as.factor("control"))
  if (i==11) {  # Remove the reference measurements made on Chamber 11
    height.dia.sub.1 = subset(height.dia.sub, Stem_number %in% 1)
    height.dia.sub.1.2 = subset(height.dia.sub, Stem_number %in% 1.2)
    height.dia.sub.1[height.dia.sub.1$DateTime >= as.Date("2013-12-24"),"Plant_height"] = height.dia.sub.1.2[,"Plant_height"]
    height.dia.sub = height.dia.sub.1
  }
  keeps <- c("chamber", "DateTime", "T_treatment", "Plant_height", "X15", "X65")
  # keeps <- c("chamber", "DateTime", "T_treatment", "Water_treatment", "Plant_height", "X15", "X65")
  height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
  
  D.15 <- lm(X15 ~ X65, data=height.dia.sub)
  # visreg(D.15, "X65", overlay=TRUE)
  # summary(D.15)
  eq.D = function(x){coefficients(D.15)[1] + coefficients(D.15)[2] * x }
  
  height.dia.sub$X15 = eq.D(height.dia.sub$X65)
  height.dia.sub = subset(height.dia.sub, DateTime >= as.Date("2013-09-14") & DateTime <= as.Date("2014-05-27"))
  height.dia.sub = height.dia.sub[!is.na(height.dia.sub$X65),]
  # keeps <- c("chamber", "DateTime", "T_treatment", "Water_treatment", "X15", "Plant_height")
  keeps <- c("chamber", "DateTime", "T_treatment", "X15", "Plant_height")
  height.dia.sub = height.dia.sub[ , keeps, drop = FALSE]
  # names(height.dia.sub) <- c("chamber","Date","T_treatment","W_treatment","diameter","height")
  names(height.dia.sub) <- c("chamber","Date","T_treatment","diameter","height")
  height.dia.sub$T_treatment = as.character(height.dia.sub$T_treatment)
  # height.dia.sub$W_treatment = as.character(height.dia.sub$W_treatment)
  
  # height.dia[(1+(i-1)*length(unique(height.dia$Date))) : (i*length(unique(height.dia$Date))), 
  #            c("T_treatment","W_treatment","diameter","height")] = height.dia.sub[,c("T_treatment","W_treatment","diameter","height")]
  height.dia[(1+(i-1)*length(unique(height.dia$Date))) : (i*length(unique(height.dia$Date))), 
             c("T_treatment","diameter","height")] = height.dia.sub[,c("T_treatment","diameter","height")]
}
height.dia$T_treatment = as.factor(height.dia$T_treatment)
# height.dia$W_treatment = as.factor(height.dia$W_treatment)

# Average the ambient and elevated temperature treatments considering the drought/watered treatment seperated from the start of the experiment
# n=3 for whole period
# drought.chamb = unique(height.dia$chamber[ height.dia$W_treatment %in% as.factor("drydown")])
# height.dia$chamber_type = as.factor( ifelse(height.dia$chamber %in% drought.chamb, "drought", "watered") )
# height.dia.final <- summaryBy(height+diameter ~ Date+T_treatment+chamber_type, data=height.dia, FUN=c(mean,standard.error))
# names(height.dia.final)[4:7] = c("height", "diameter", "height_SE", "diameter_SE")
height.dia.final <- summaryBy(height+diameter ~ Date+T_treatment, data=height.dia, FUN=c(mean,standard.error))
names(height.dia.final)[3:6] = c("height", "diameter", "height_SE", "diameter_SE")

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
