#- function to load a package, and install it if necessary
Library <- function(pkg, ...){
  
  PACK <- .packages(all.available=TRUE)
  pkgc <- deparse(substitute(pkg))
  
  if(pkgc %in% PACK){
    library(pkgc, character.only=TRUE)
  } else {
    install.packages(pkgc, ...)
    library(pkgc, character.only=TRUE)
  }
  
}
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#- load all the libraries (and install them if needed)
Library(mvtnorm) # Creates candidate parameter vector as a multivariate normal jump away from the current candidate
Library(reshape2)
Library(ggplot2)
Library(lubridate)
Library(rio)
Library(dplyr)
Library(zoo)
Library(doBy)
Library(corrplot)
Library(png)
Library(grid)
Library(gridExtra)
Library(RColorBrewer)
Library(sciplot)
Library(reshape)
Library(nlme)
Library(visreg)
Library(multcomp)
Library(lme4)
Library(lmerTest)
Library(plotrix)
Library(plantecophys)
Library(scales)
Library(effects)
Library(plyr)
Library(car)
Library(lattice)
Library(data.table)
Library(grDevices)
Library(shape)
Library(forecast)
Library(mgcv)
Library(scales)
Library(gplots)
Library(magicaxis)
Library(Hmisc)
Library(hexbin)
Library(lsmeans)
Library(calibrate)
Library(MuMIn)
Library(rmarkdown)
Library(knitr)
Library(imputeTS)
Library(gdata)
Library(parallel)
Library(snow)
Library(broom)
Library(chron)

#- the following libraries aren't on CRAN, but can be installed from github or bitbucket with devtools
if (require("devtools")==F) {install.packages("devtools")
  library(devtools,quietly=T)}
if (require("plantecophys")==F) {
  install_bitbucket("remkoduursma/plantecophys")
  library(plantecophys,quietly=T)}
if (require("HIEv")==F) {
  install_bitbucket("remkoduursma/HIEv")
  library(HIEv,quietly=T)}
if (require("plotBy")==F){ 
  install_bitbucket("remkoduursma/plotBy")
  library(plotBy,quietly=T)}

# font_import(pattern="[A/a]rial", prompt=FALSE)

#- check if the data and output directories exist. If they don't, create them.
dir.create(file.path("data"),showWarnings=F)
dir.create(file.path("processed_data"),showWarnings=F)
dir.create(file.path("output"),showWarnings=F)
