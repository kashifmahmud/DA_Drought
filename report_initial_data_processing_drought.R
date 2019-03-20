#' ---	
#' title: "Processing WTC-3 data to infer temperature effect on Carbon Balance"	
#' author: "Kashif Mahmud"	
#' date: "27 March 2018"	
#' output:	
#'   html_document: default	
#'   word_document: default	
#' ---	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' ### *R script to process the raw WTC-3 data to estimate the carbon pools and fluxes for Data Assimilation (DA)*  	
#' 	
#' <!-- Load required libraries. There are quite a few, including some non-standard functions that are not on CRAN. This script will check for required libraries and install any that are missing.   -->	
#' 	
source('R/load_packages_wtc3.R')	
#' 	
#' 	
#' <!-- Load the custom analysis and plotting functions that do all of the actual work   -->	
#' 	
source("R/functions_wtc3.R")	
source("R/functions_wtc3_CBM.R")	
#' 	
#' 	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' ##### Briefly describe the Carbon Balance Model (CBM)    	
#' We planned to use a DA-modelling framework, similar to that used by Richardson et al. (2013) and Mahmud et al. (in preparation). This approach uses a simple carbon balance model shown in **Figure 1**. The model is driven by daily inputs of gross primary production (GPP). Total maintenance respiration (R~m,tot~) is immediately subtracted, and the remainder enters a non-structural C pool (C~n~). This pool is utilized for growth at a rate *k*. Of the utilization flux, a fraction *Y* is used in growth respiration (R~g~), and the remaining fraction (1-*Y*) is allocated to structural C pools (C~s~): among foliage, wood and root (C~s,f~, C~s,w~, C~s,r~). The foliage and root pools are assumed to turn over with rate *s~f~* and *s~r~* respectively. We assume there is no wood turnover as evidenced from the experiment.   	
#' 	
#' **Figure 1**: Structure of the Carbon Balance Model (CBM). Pools, shown as solid boxes: C~n~, non-structural storage C; C~s,f~, structural C in foliage; C~s,r~, structural C in roots; C~s,w~, structural C in wood. Dummy pools, shown as dashed boxes: C~f,lit~, total C in foliage litterfall; C~r,lit~, total C in root turnover. Fluxes, denoted by arrows, include: GPP, gross primary production; R~m,tot~, maintenance respiration; Rg, growth respiration. Fluxes are governed by seven key parameters: *k*, storage utilization coefficient; *Y*, growth respiration fraction; *a~f~*, allocation to foliage; *a~w~*, allocation to wood; $a_r = (1-a_f-a_w)$, allocation to roots; *s~f~*, foliage turnover rate; *s~r~*, root turnover rate.    	
#' 	
#' 	
plot.model.wtc3()	
#' 	
#' 	
#' The dynamics of the four carbon pools (C~n~, C~s,f~, C~s,w~, and C~s,r~) are described by four difference equations:   	
#' $$\Delta C_n = GPP - R_{m,tot} - kC_n$$	
#' $$\Delta C_{s,f} = kC_n(1-Y)a_f - s_fC_{s,f}$$	
#' $$\Delta C_{s,w} = kC_n(1-Y)a_w$$	
#' $$\Delta C_{s,r} = kC_n(1-Y)(1-a_f-a_w) - s_rC_{s,r}$$	
#' 	
#' Where *k* is the storage utilization coefficient; *Y* is the growth respiration fraction; *a~f~*, *a~w~*, *a~r~* are the allocation to foliage, wood and root respectively; *s~f~* and *s~r~* are the leaf and root turnover rates respectively; $\sum s_fC_{s,f} = C_{f,lit}$ and $\sum s_rC_{s,r} = C_{r,lit}$ are the dummy foliage and root litter pools respectively.   	
#' 	
#' Total maintenance respiration, R~m,tot~ was calculated as a temperature-dependent respiration rates for foliage, wood and root (R~m,f~, R~m,w~ and R~m,r~ respectively), multiplied by plant organ C masses (C~t,f~, C~t,w~, and C~t,r~ are the total C in foliage, wood and root respectively). Wood respiration was further partitioned into stem and branches. Similarly root respiration was partitioned into different root size classes (fine, intermediate, coarse and bole roots). Growth respiration (R~g~) at each time step was calculated as a modeled respiration rate (*Y*) multiplied by plant biomass changes at that time step ($\Delta C_{t,f} + \Delta C_{t,w} + \Delta C_{t,r}$).   	
#' 	
#' $$R_{m,tot} = R_{m,f} C_{t,f} + R_{m,w} C_{t,w} + R_{m,r}C_{t,r}$$	
#' $$R_{m,w} C_{t,w} = R_{m,s} C_{t,s} + R_{m,b} C_{t,b}$$	
#' $$R_{m,r} C_{t,r} = R_{m,fr} C_{t,fr} + R_{m,ir} C_{t,ir} + R_{m,cr} C_{t,cr} + R_{m,br} C_{t,br}$$	
#' $$R_g = Y(\Delta C_{t,f} + \Delta C_{t,w} + \Delta C_{t,r})$$	
#' 	
#' Where R~m,s~, R~m,b~, R~m,fr~, R~m,ir~, R~m,cr~ and R~m,br~ are the maintenance respiration rates with C~t,s~, C~t,b~, C~t,fr~, C~t,ir~, C~t,cr~ and C~t,br~ are the total C in stem, branches, fine roots, intermediate roots, coarse roots and bole roots respectively.    	
#' 	
#' The non-structural (storage) C pool (C~n~) is assumed to be divided among foliage, wood and root tissues (C~n,f~, C~n,w~, C~n,r~) according to empirically-determined fractions.   	
#'   - *However, WTC-3 experiment only measured leaf non-structural C (C~n,f~), and therefore to estimate the partitioning of the non-structural C among different organs, we used data from WTC-4 experiment on similar-sized seedlings of a related species (Eucalyptus parramattensis). We will consider different treatments from the experimental dataset, to find the C~n~ partitioning to foliage, wood and roots.*   	
#'   - *NOT NEEDED: Another possibility (if WTC-4 data are not available!!) would be to assume another two parameters to estimate the C~n~ partitioning to foliage, wood and roots (p~f~, C~n~ partitioning to foliage; p~w~, C~n~ partitioning to wood; $p_r = (1-p_f-p_w)$, C~n~ partitioning to roots).*   	
#'   	
#' Total carbon in each tissue (C~t~) is then calculated as the sum of non-structural carbon (C~n~) and structural carbon (C~s~) for that tissue.    	
#' 	
#' $$C_{t,f} = C_{n,f} + C_{s,f}$$	
#' $$C_{t,w} = C_{n,w} + C_{s,w}$$	
#' $$C_{t,r} = C_{n,r} + C_{s,r}$$	
#' We plan to estimate seven parameters (*k*, *Y*, *a~f~*, *a~w~*, *a~r~*, *s~f~*, *s~r~*) of the CBM for this experiment using DA. GPP and maintenance respiration rates will be used as model inputs in the DA framework, whereas the measurements of total aboveground respiration (R~above~), total C masses (C~t,f~, C~t,w~, C~t,r~), foliage litterfall (C~f,lit~) and foliage NSC (C~n,f~) will serve to constrain the parameters. Aboveground respiration (R~above~) was estimated by summing both maintenance and growth respiration components of foliage and wood.      	
#' $$R_{above} = R_{m,f} C_{t,f} + R_{m,w} C_{t,w} + Y\Delta C_{t,f} + Y\Delta C_{t,w}$$	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' #### Step 1:  	
#' ##### Select treatment strategy  	
#' - Check whether there are any pre-existing differences between the chambers of droughted treatments  	
#' - Should we consider:  	
#'     1. the droughted treatments separately from the start of the flux measurements (Sept 2013): n = 3 or   	
#'     2. from actual drought implementation on 12 Feb 2014: n = 6 predrought and n = 3 postdrought ??   	
#' 	
#' ##### Data used: 	
#' 1. [WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV](https://hiev.uws.edu.au/data_files/56022)   	
#' 	
#' **Figure 2**: Growth (Diameter) of Eucalyptus tereticornis trees exposed to warming and drought    	
#' 	
#' 	
source('R/treatment_selection_strategy.R')	
#' 	
#' 	
#' #### Finding:  	
#' It is obvious from Figure 1 that there were pre-existing differences between the warmed watered (black vs blue, Figure 1) and warmed drought (black vs pink, Figure 1) trees. It seems like the warmed-drought trees actually started the drought smaller than the warmed-watered trees (blue vs pink, Figure 1). If we go for the above option 2 (i.e. separate all 4 treatments from actual drought implementation on 12 Feb 2014), the diameter time series (red line, Figure 1) makes it look like there was a strong effect of the drought on the tree growth but it was more a function of pre-existing differences in tree size. The height data showed the similar pattern.    	
#' 	
#' #### So the idea is:    	
#' Based on both diameter and height data, we decided to split the whole experiment into two time frames:   	
#' 1. *Case 1*: First one tries to infer solely the temperature effect taking account of the measurements from the start of flux measurements (17 Sept 2013) till the start of drought treatment (11 Feb 2014). So total number of treatments = 2 (ambient and warmed) with sample size of n = 6 for each treatment.      	
#' 2. *Case 2*: Second part examines the effect of both temperature and drought simultaneously, considering the data from the start of drought treatment (12 Feb 2014) to the end of the experiment (26 May 2014). So total number of treatments = 4 (ambient+watered, ambient+droughted, warmed+watered, warmed+droughted) with sample size of n = 3 for each treatment.    	
#' 	
#' #### Here we focus on *Case 1*: processing the WTC-3 data from the start of flux measurements (17 Sept 2013) till the start of drought treatment (11 Feb 2014)          	
#' 	
#' 	
#' ##### Alternative thought:    	
#' Split the whole experiment into two time frame:   	
#' 1. Case 1: First one tries to infer solely the temperature effect taking account of the measurements from the start of flux measurements (17 Sept 2013) to the start of drought treatment (12 Feb 2014)   	
#' 2. Case 2: Second part examines the effect of both temperature and drought simultaneously, considering the data from the start of drought treatment (12 Feb 2014) to the end of the experiment (26 May 2014)    	
#' 	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#'    	
#' #### Step 2:    	
#' ##### Inputs:   	
#' 	
c1 = 0.48 # (unit conversion from gDM to gC: 1 gDM = 0.48 gC)   	
#' 	
#' 	
#' ##### Estimate the partitioning of the non-structural C among different organs using data from WTC-4 experiment on similar-sized seedlings of a related species (Eucalyptus parramattensis).    	
#' 	
#' 1. We considered different treatments (Ambient/Warmed) from the experimental dataset, to find the C~n~ partitioning to foliage, wood and roots.       	
#' 2. There was no statistically significant difference in TNC partitioning across the treatments.     	
#' 	
#' 	
source('R/tnc_partitioning_wtc4.R')	
#' 	
#' 	
#' <!-- write the file of TNC partitioning from WTC-4 dataset  -->	
#' 	
write.csv(tnc.partitioning, file = "processed_data/tnc_partitioning_data.csv", row.names = FALSE)	
#' 	
#' 	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#'    	
#' #### Step 3: Read and process GPP and R~above~ (aboveground respiration) data   	
#' 	
#' ##### Partitioning CO~2~ fluxes into GPP and R~above~	
#' Partitioning of the hourly net CO~2~ fluxes into the components of GPP and R~above~ were done using a technique common to eddy-covariance research (Reichstein et al. 2005); described thoroughly in Drake et al. (2016, 2017-in preparation). We assumed GPP to be zero at night when PPFD = 0, indicating the measured net C flux in such conditions was used as the measure of R~above~. We utilized the direct measurements of whole-canopy R~above~ and its temperature dependence at night to predict R~above~ for each hourly measurement as a function of air temperature. We then calculated GPP as the sum of the measured net CO~2~ flux and the predicted R~above~, given the measured air temperature.    	
#' 	
#' ##### Data used: 	
#' 1. [WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv](https://hiev.uws.edu.au/data_files/115778)   	
#' 	
#' **Figure 3**: Daily GPP and aboveground respiration (R~above~) of Eucalyptus tereticornis trees for all 4 treatments. The grey bars represent standard errors for sample size, n=3.    	
#' 	
#' 	
source('R/GPP_Ra_LA_data.R')	
#' 	
#' 	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#'     	
#' #### Step 4: Estimate biomass pools for each treatment   	
#' 1. Process stem height and diameter for various treatment cases  	
#'   - Tree height (cm) and stem diameter (mm) were measured fortnightly. Stem diameter was measured at 30-cm-intervals along each tree stem from a basal height of 15-cm (prior to floor installation) or 65-cm (after floor installation) to the tree apex. However, the diameters at 15-cm height are essential to estimate the initial root mass (at the start of the flux measurement), based on the allometry and harvest root mass of Court's experiment. Therefore, we gap-filled the 15-cm diameters by fitting a linear regression between diameters at 15 and 65 cm heights.    	
#'   - Finally we averaged both the height and diameter (at 15-cm) across the treatments (n=3) to get the mean fortnightly allometry data with standard errors.       	
#'   	
#' 2. Estimate the root mass for each treatment    	
#'   - We estimated initial root mass (g C) using allometry (height and diameter at 15-cm height) and harvest root mass of control (free) seedlings from Court's experiment on similar seedlings (Euc. Teri.). We read and merged all necessary data from Court's experiment to get the allometry and harvest root mass. We considered both data from Court's experiment and WTC-3 harvest for the control (ambient and watered) treatments, and fitted a linear regression (with interaction) to estimate root mass (with standard errors) from the initial data of stem diameter and tree height at the start of the flux measurements:	
#'   lm(log10(rootmass) ~ log10(diameter) * log10(height), P < 0.001, r^2^ = 0.99.     	
#'   - We averaged WTC-3 final harvest root mass across the treatments (n=3) to get the mean and standard errors.   	
#'   	
#' 3. Get an estimate of stem, branch and leaf mass for each treatment following the script of John for WTC-3 flux partitioning paper (in progress)  	
#'   - Fortnightly stem mass (g C) was calculated from geometry using stem diameter at different heights. We considered stem as i) frustum of a cone above 65 cm, and ii) cylinder below 65 cm to ground.    	
#'   - Branch mass was estimated using harvest branch mass, branch census and stem volume. Using the harvest data, we setup a log-log regression for branch mass and diameter:	
#'   (log10(branch mass) = -1.299 + 2.722 × log10(branch diameter), P < 0.001, r^2^ = 0.91, branch mass in g, branch diameter in mm, n = 48 branches).   	
#'   - Then we applied this allometry-mass regression to the branch census dataset to estimate the branch mass on three intermediate dates (24 Oct 2013, 15 Jan 2014, and 22 May 2014). Total branch mass and stem volume were strongly correlated in a chamber specific manner (log-log ANCOVA, P < 0.001, r^2^ = 0.95), which was used to estimate branch mass (g C) as a function of stem volume for each fortnightly growth interval.   	
#'   - We estimated fortnightly leaf mass using canopy-weighted SLA (Specific leaf area) from harvest and leaf area (LA) estimates according to Barton et al. (2012) and Drake et al. (2016). LA was based on stem height and leaf litterfall at any time. Standing LA was measured twice on 9 Sept 2013 and 10 Feb 2014 for each tree by counting all the leaves and multiplying by a tree-specific mean leaf size measured across the canopy of each tree with a handheld leaf area meter. A third direct measurement of standing canopy leaf area was calculated from the final harvest data (26 May 2014) by multiplying total canopy leaf dry mass by SLA weighted by the leaf dry mass in each canopy layer.       	
#' 	
#' ##### Data used:   	
#' 1. [WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV](https://hiev.uws.edu.au/data_files/56022)   	
#' 2. [WTC_TEMP_CM_HARVEST-ROOTS_20140529-20140606_L1_v1.csv](https://hiev.uws.edu.au/data_files/55651)   	
#' 3. [WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv](https://hiev.uws.edu.au/data_files/115778)	
#' 4. [WTC_TEMP_CM_HARVEST-CANOPY_20140526-20140528_L1_v1.csv](https://hiev.uws.edu.au/data_files/55645)	
#' 5. [WTC_TEMP_CM_BRANCHCENSUS_20130910-20140516_L0_v1.csv](https://hiev.uws.edu.au/data_files/58633)  	
#' 6. [Below-ground Sink limited container experiment (Court's exp) - Harvest diameter, height and rootmass data of free seedlings](https://hiev.uws.edu.au/data_files/174409)   	
#' 	
#' ##### Assumptions:  	
#' 1. No LMA variation over time: We plotted the LMA variation over time to see whether there was any time effect on LMA. We used all available data from Mike, Court and harvest to investigate the trend, however did not find any concrete relationship. Therefore, we ended up with the constant LMA from harvest data as used by John and suggested by Remko, Court.  	
#' 	
#' <!-- Script to process stem height and diameter for various treatment cases  -->	
#' 	
source('R/stem_height_diameter_data.R')	
#' 	
#' 	
#' <!-- Get the root mass (initial - estimated and final - harvested) for each treatment  -->	
#' 	
source('R/rootmass_data.R')	
#' 	
#' 	
#' <!-- Get an estimate of branch, stem, and leaf mass as well as leaf area for each treatment  -->	
#' 	
source('R/branch_stem_leafmass_data.R')	
#' 	
#' 	
#' <!-- Merge and plot all biomass data  -->	
#' **Figure 4**: Estimated biomass pool of WTC-3 Eucalyptus tereticornis trees for all 4 treatments. We predicted total C mass present in both foliage and wood with fortnightly interval, however only the initial and final C mass for roots. Note that, the points are jittered in all figures to see the grey standard error bars (n=3).        	
#' 	
#' 	
source('R/plot_biomass_data.R')	
#' 	
#' 	
#' 	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' #### Step 5: Process the leaf litter data   	
#' Litterfall was collected, dried, and weighed fortnightly for each tree. The dummy foliage litter pool, C~f,lit~ (g C) was estimated by cumulative summing of fortnightly litter starting from the flux measurements. The daily foliage litterfall was predicted as a linear function of time between two fortnightly consecutive measurements and will be used to constrain the parameter foliage turnover rate, *s~f~*.          	
#' 	
#' ##### Data used:   	
#' 1. [WTC_TEMP_CM_LEAFLITTER_20130913-20140528_L1.csv](https://hiev.uws.edu.au/data_files/55646)   	
#' 	
#' <!-- Process leaf litter data for each treatment  -->	
#' 	
source('R/leaf_litter_data.R')	
#' 	
#' 	
#' <!-- Merge leaf litter with all biomass data and plot the leaf litterfall  -->	
#' **Figure 5**: Dummy leaf litterfall pool of WTC-3 Eucalyptus tereticornis trees for all 4 treatments. Note that, the points are jittered to see the grey standard error bars (n=3).        	
#' 	
#' 	
source('R/plot_leaf_litter_data.R')	
#' 	
#' 	
#' 	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' #### Step 6: Estimate the leaf TNC (storage) pool   	
#' 1. First we considered the foliage TNC measurements over time and averaged the data across the treatments (n=3) to get the mean storage with standard errors.    	
#' 2. We then merged the diurnal foliage TNC data by averaging all 2-hourly interval data on a sunny (20 February) and an overcast/rainy day (26 March). Similarly this data were averaged across the treatments (n=3) to get the mean with standard errors.       	
#' 3. We also took onto board the diurnal leaf TNC data when trees were girdled on 15 May 2014 (only the ambient watered treatments).     	
#' 4. The unit was converted to gC in TNC per gC of plant for consistency.    	
#' 	
#' ##### Data used:   	
#' 1. [WTC_TEMP_CM_LEAFCARB_20130515-20140402_L2.csv](https://hiev.uws.edu.au/data_files/95740)   	
#' 2. [WTC_TEMP_CM_LEAFCARB-DIURNAL_20140220-20140326_R.csv](https://hiev.uws.edu.au/data_files/95732)      	
#' 3. [WTC_TEMP_CM_PETIOLEGIRDLE-LEAFMASS-AREA-CARB_20140507-20140512_L2.csv](https://hiev.uws.edu.au/data_files/96526)   	
#' 	
#' <!-- Process leaf TNC data for each treatment  -->	
#' 	
source('R/leaf_tnc_data.R')	
#' 	
#' 	
#' <!-- Merge leaf TNC with all biomass data and plot the leaf tnc storage data  -->	
#' **Figure 6**: Foliage TNC pool, *C~n,f~* over time of WTC-3 Eucalyptus tereticornis trees for all 4 treatments. Note that, the points are jittered to see the grey standard error bars (n=3).        	
#' 	
#' 	
source('R/plot_leaf_tnc_data.R')	
#' 	
#' 	
#' <!-- write the file of Daily GPP, Ra, LA and all mass pools (rootmass, woodmass, foliagemass, litterfall, tnc) data  -->	
#' 	
write.csv(data.GPP.Ra.LA.mass, file = "processed_data/GPP_Ra_LA_mass.csv", row.names = FALSE)	
#' 	
#' 	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' #### Step 7: Estimate the partitioning of woodmass to branch/stem and rootmass to fine/intermediate/coarse/bole roots  	
#' ##### *Purpose: To predict total wood and root maintenance respiration $R_{m,w}C_{t,w}$ and $R_{m,r}C_{t,r}$, we need both daily respiration rate and corresponding daily C mass of all wood and root components.*    	
#' ##### Woodmass partitioning    	
#' 1. Calculate partitioning of WTC-3 woodmass: We did not find any significant difference across the treatments for wood partitioning, so we considered equivalent wood partitioning for all treatments (Figure 7, solid lines).      	
#' 2. We then considered a linear variation over time in between the fortnightly measurements to predict the daily woodmass partitioning.     	
#' 	
#' ##### Rootmass partitioning     	
#' 1. Calculate partitioning of WTC-3 harvest rootmass: We did not find any significant difference across the treatments for root partitioning, so we considered equivalent root partitioning for all treatments (Figure 7, dotted lines end points).      	
#' 2. Estimate the partitioning of initial rootmass based on Court's Pot experiment harvest rootmass (control free seedling): We assume there were only fine and intermediate roots, with no coarse and bole roots at the start of flux measurements (Figure 7, starting points).   	
#' 3. We then considered a linear variation over time to predict the daily rootmass partitioning, in agreement with above ground biomass partitioning that follows linear trend.     	
#' 	
#' ##### Data used:   	
#' 1. [WTC_TEMP_CM_LEAFCARB_20130515-20140402_L2.csv](https://hiev.uws.edu.au/data_files/95740)   	
#' 2. [Below-ground Sink limited container experiment (Court's exp) - Harvest rootmass data of free seedlings](https://hiev.uws.edu.au/data_files/174409)     	
#' 	
#' <!-- Calculate partitioning of WTC-3 initial and harvest rootmass, and then plot the linear time series -->  	
#' **Figure 7**: Rootmass partitioning of WTC-3 Eucalyptus tereticornis trees for all 4 treatments. The grey shade shows the standard error (n=12).        	
#' 	
#' 	
source('R/woodmass_rootmass_partitioning.R')	
#' 	
#' 	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' #### Step 8: Plot leaf-scale respiration vs. leaf temperature to find Q~10~ for WTC-3    	
#' 1. We plotted the leaf-scale respiration (R~leaf~) against leaf temperature (T~leaf~) to test the short-term temperature sensitivity of the respiration of individual leaves similar to Drake et al. (2016). However, instead of using Arrhenius function (Drake et al. 2016), we fitted the exponential relationship: $R_{leaf} = R_{leaf,25}Q_{10}^{(T_{leaf}-25)/10}$. The statistical test did not show any significant differences for drought treatment, however warming reduced the pre-exponential coefficient (*R~leaf,25~*) but did not alter *Q~10~*, reflecting a *Q~10~* value of **2.26** at 25°C for all treatments that agrees with previous study of Drake et al. (2016).        	
#' 	
#' ##### Data used:    	
#' 1. [WTC_TEMP_CM_GX-RdarkVsT_20140207-20140423_L1.csv](https://hiev.uws.edu.au/data_files/53862)   	
#' 	
#' **Figure 8**: The short-term temperature sensitivity of the respiration of individual leaves relative to leaf temperature. Error bars reflect the standard errors of Eucalyptus tereticornis trees exposed to each temperature treatments (n = 6).    	
#' 	
#' 	
source('R/q10_calculate.R')	
#' 	
#' 	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' #### Step 9: Relationship between mass-based leaf-scale respiration at 25°C and air temperature        	
#' 1. We plotted the mass-based leaf-scale respiration at 25°C, R~leaf,25~ (nmol CO~2~ g^-1^ s^-1^) over time following Aspinwall et al. (2016). **Figure 9** (first panel) shows that leaf R (mass-basis) measured at 25°C varied over time and decreased as mean air temperature (T~air~) increased (rest of the panels).    	
#' 2. Using the site weather data, we fitted a linear regression between R~leaf25~ and T~air~. We tried different options with mean air temperature of previous 7-days, previous 3-days, previous day, 3-day including the day of respiration measurement and only the day of respiration measurement. Treatment had no effect on the intercept or slope of the relationship between R~leaf25~ and T~air~. Based on the data fit and considering the finding of previous study (Aspinwall et al. 2016), we found the best relationship of 25°C mass-based leaf-scale respiration (nmol CO~2~ g^-1^ s^-1^) with 3-day mean air temperature, T~air,3day~ (°C) with P < 0.001, r^2^ = 0.75:   	
#' $$R_{leaf,25} = 21.5 - 0.565 T_{air,3day}$$	
#' 	
#' ##### Data used:    	
#' 1. [WTC_TEMP_CM_GX-Rdark25_20130617-20140402_L2.csv](https://hiev.uws.edu.au/data_files/55626)    	
#' 2. [WTC_TEMP_CM_WTCMET_20130601-20130630_L1_v1.csv](https://hiev.uws.edu.au/data_files/52258) - [WTC_TEMP_CM_WTCMET_20140501-20140531_L1_v1.csv](https://hiev.uws.edu.au/data_files/52269)     	
#' 	
#' **Figure 9**: Mass-based night-time leaf-scale dark respiration measured in Eucalyptus tereticornis at a set temperature of 25°C through time (first panel) and in relation to prevailing mean air temperature, T~air~ (rest of the panels). Mean values are those of replicate whole-tree chambers (n = 6 for ambient/warmed treatment), determined based on three individual leaf measures per chamber. No statistically significant variation for drought/watered treatment on leaf respiration.           	
#' 	
#' 	
source('R/leafscaleR_vs_Tair_regression.R')	
#' 	
#' 	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' #### Step 10: Estimate daily mean respiration rates for foliage, wood and roots          	
#' 1. Foliage respiration: First we calculated the daily mean mass-based leaf-scale respiration, R~leaf,25~ (nmol CO~2~ g^-1^ s^-1^) at 25°C using the regression we fitted with 3-day mean air temperature, T~air,3day~ (°C) in Step 9: $R_{leaf,25} = 21.5 - 0.565 T_{air,3day}$.      	
#' 2. Bole wood respiration: There was no statistically significant difference in bole/stem wood respiration across the treatments, and the rate was 0.656 nmol CO~2~ g^-1^ s^-1^ at 15°C.                   	
#' 3. Branch wood respiration: We found statistically significant difference in branch wood respiration across the temperature treatments (ambient vs. warmed). The rates of branch wood respiration in WTC-3 were about 1.33 nmol CO~2~ g^-1^ s^-1^ in the ambient treatment and about 1.04 nmol CO~2~ g^-1^ s^-1^ in the warmed (+3°C) treatment when measured at 15°C. This indicates temperature acclimation in branch wood respiration.     	
#' 4. Fine root respiration: We assumed the specific rate of fine root respiration at 25°C was about 10 nmol CO~2~ g^-1^ s^-1^ (Drake et al. 2017). We estimated fine root respiration at 15°C using the Q~10~ value of 2.26 (as found in Step 8).                    	
#' 5. Coarse root respiration: For coarse roots (Diameter > 10 mm), we suggested using the rate we determined for branch wood in E. tereticornis. The rationale would be that the woody coarse roots were made up of more or less well developed xylem (and phloem) tissues, were thus more or less functioning as underground woody branches. Therefore, we would expect similar temperature acclimation in coarse root as in branch wood respiration.     	
#' 6. Bole root respiration: The bole and big tap roots were assigned the bole wood respiration, collected at the end of the WTC-3 study at 15°C.      	
#' 7. Intermediate root respiration: For the intermediate size class of 2 to 10 mm, We interpolated respiration rates between the size classes to estimate the middle class using a simple log-log plot.     	
#' 	
#' 8. **15-mins interval respiration rates**:     	
#'   - *We estimated 15-mins interval foliage respiration rates, R~leaf~ at 15-mins air temperatures (T~air~) using the site weather data, mass-based daily mean leaf-scale respiration at 25°C, R~leaf,25~ and Q~10~ of 2.26: $R_{leaf} = R_{leaf,25}Q_{10}^{(T_{leaf}-25)/10}$. We assume $T_{leaf} = T_{air}$. We then consider a 30% reduction in foliage respiration during day time.*    	
#'   - *Wood respiration rates in 15-mins interval were calculated at air temperature (T~air~) using the site weather data, mass-based wood respiration rates at 15°C and Q~10~ of 2.26.*    	
#'   - *Similarly root respiration rates in 15-mins interval were calculated using soil temperature at 10 cm depth, mass-based root respiration rates at 15°C and Q~10~ of 2.26.*    	
#' 9. **Daily mean respiration rates**: *Daily mean respiration rates for all tree components were calculated by summing all 15-mins data for each day.*    	
#' 10. *The units were converted to gC g^-1^C d^-1^ (from nmol CO~2~ g^-1^ s^-1^) for consistency.*    	
#' 	
#' ##### Data used:    	
#' 1. [WTC_TEMP_CM_GX-RBRANCH_20140513-20140522_L1_v1.csv](https://hiev.uws.edu.au/data_files/55727)     	
#' 2. [WTC_TEMP_CM_WTCFLUX-STEM_20140528_L1_v1.csv](https://hiev.uws.edu.au/data_files/179783)     	
#' 3. [WTC_TEMP_CM_WTCMET_20130601-20130630_L1_v1.csv](https://hiev.uws.edu.au/data_files/52258) - [WTC_TEMP_CM_WTCMET_20140501-20140531_L1_v1.csv](https://hiev.uws.edu.au/data_files/52269)	
#' 	
#' <!-- Calculate partitioning of WTC-3 initial and harvest rootmass, and then plot the linear time series -->  	
#' **Figure 10**: Daily mean respiration rates for all tree components (foliage, wood and roots). There were no statistically significant difference in respiration rates across drought/watered treatments. Both branch wood and coarse root respiration showed temperature acclimation having higher rates in ambient than warmed condition.               	
#' 	
#' 	
source('R/daily_R_rates.R')	
#' 	
#' 	
#' <!-- write csv file of all daily data (with woodmass and rootmass partitioning, respiration rates for roots, wood and foliage -->	
#' 	
write.csv(data.all, "processed_data/data_all.csv", row.names=FALSE) # unit of respiration rates: gC per gC plant per day	
#' 	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' #### Step 11: Apply Data Assimilation (DA) with the estimates of WTC-3 carbon pools and fluxes to predict all seven parameters (*k*, *Y*, *a~f~*, *a~w~*, *a~r~*, *s~f~*, *s~r~*), with the aim to quantify the main C balance processes (respiration, carbohydrate utilisation, allocation and turnover) in response to elevated temperature and drought on Eucalyptus tereticornis trees grown in WTCs.     	
#' 	
#' Hypothesis 1: NSC vs. Temperature	
#'   - As *Paul, Driscoll & Lawlor (1991)* have shown, low T generally results in accumulation of carbohydrates, indicating that growth or storage are limiting. Conversely, at high temperatures, sink demand is large and assimilate is depleted so there should be a marked higher utilization rate (*k*).    	
#'   - The carbohydrate concentration decreased in stem and root tissues for Citrus plants, while it increased in leaf tissues under moderate warm conditions (30/20ºC than at 25/20ºC, *Ribeiro et al. 2012*).     	
#' Hypothesis 2: NSC vs. Drought	
#' - We hypothesize that drought is carbon limiting and negatively impacts plant carbon balance and that plants will rely on stored carbon to survive carbon limitation and therefore deplete their stored carbon reserves over time.     	
#' 	
#' 	
#' 	
#' ##### *... in progress ...*       	
#' 	
#' 	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' <!-- #---------------------------------------------------------------------------------------------------------------- -->	
#' 	
#' 	
#' 	
#' 	
#' 	
