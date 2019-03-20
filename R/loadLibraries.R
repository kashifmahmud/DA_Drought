#- required libraries
if (require("shape",quietly=T)==F) {
  install.packages("shape")
  library(shape,quietly=T)}

if (require("forecast",quietly=T)==F) {
  install.packages("forecast")
  library(forecast,quietly=T)}

if (require("mgcv",quietly=T)==F) {
  install.packages("mgcv")
  library(mgcv,quietly=T)}

if (require("scales",quietly=T)==F) {
  install.packages("scales")
  library(scales,quietly=T)}

if (require("gplots",quietly=T)==F) {
  install.packages("gplots")
  library(gplots,quietly=T)}

if (require("magicaxis",quietly=T)==F) {
  install.packages("magicaxis")
  library(magicaxis,quietly=T)}

if (require("lubridate",quietly=T)==F){
  install.packages("lubridate")
  library(lubridate,quietly=T)} 

if (require("doBy",quietly=T)==F) {
  install.packages("doBy")
  library(doBy,quietly=T)}

if (require("Hmisc",quietly=T)==F) {
  install.packages("Hmisc")
  library(Hmisc,quietly=T)}

if (require("zoo",quietly=T)==F) {
  install.packages("zoo")
  library(zoo,quietly=T)}

if (require("hexbin",quietly=T)==F) {
  install.packages("hexbin")
  library(hexbin,quietly=T)}

if (require("nlme",quietly=T)==F) {
  install.packages("nlme")
  library(nlme,quietly=T)}

if (require("lsmeans",quietly=T)==F) {
  install.packages("lsmeans")
  library(lsmeans,quietly=T)}

if (require("car",quietly=T)==F) {
  install.packages("car")
  library(car,quietly=T)}

if (require("data.table",quietly=T)==F) {
  install.packages("data.table")
  library(data.table,quietly=T)}

if (require("calibrate",quietly=T)==F) {
  install.packages("calibrate")
  library(calibrate,quietly=T)}

if (require("MuMIn",quietly=T)==F) {
  install.packages("MuMIn")
  library(MuMIn,quietly=T)}

if (require("dplyr",quietly=T)==F) {
  install.packages("dplyr")
  library(dplyr,quietly=T)}

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


#- check if the data and output directories exist. If they don't, create them.
dir.create(file.path("data"),showWarnings=F)
dir.create(file.path("output"),showWarnings=F)
