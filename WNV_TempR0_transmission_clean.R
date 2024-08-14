#clear memory:
rm(list=ls())

require('glmmTMB')
require('MuMIn')
require('tidyverse')
require('MASS')
require('performance')
require('DHARMa')
library('sjPlot')
library('sjmisc')
library('sjlabelled')
library('ggstats')
library('visreg')

#########################################################
### Pull in WNV data: 
#########################################################
WNV_dat <- read.csv("WNV_data.csv", header=T)
head(WNV_dat)
lapply(WNV_dat, class)
WNV_dat$NLCD_majority2 <- as.factor(WNV_dat$NLCD_majority2)

#########################################################
# base model:
#########################################################
fit_tweedie_base <- glmmTMB(
  MIR_Tar.mean ~
    trans_mean_std 
  + offset(NumObservations.sum_std)
  ,
  data=WNV_dat,
  ziformula= ~0,
  family=tweedie(link="log"),
  na.action=na.exclude
)
summary(fit_tweedie_base)
check_collinearity(fit_tweedie_base)
#DHARMa diagnostics:
testDispersion(fit_tweedie_base)
simulationOutput <- simulateResiduals(fittedModel = fit_tweedie_base, plot = F)
residuals(simulationOutput)
plot(simulationOutput)
plotResiduals(simulationOutput, form = WNV_dat$trans_mean_std)

#residual spatial autocorrelation: yes, still residual autocorrelation
testSpatialAutocorrelation(simulationOutput, x =  WNV_dat$xcoord.x, y = WNV_dat$ycoord.y)

#########################################################
# full model:
# for backward stepwise selection, include individual bird species variables;
# for dredge approach, replace these with the composite competent bird abundance variable
#########################################################
fit_tweedie_full <- glmmTMB(
  MIR_Tar.mean ~
    trans_mean_std 
  + amecro_abundance.mean_std # individual species in place of "competent_bird_abundance.mean_std" for stepwise selection
  + amerob_abundance.mean_std
  + houfin_abundance.mean_std
  + houspa_abundance.mean_std
  + cowscj_abundance.mean_std
  #+ competent_bird_abundance.mean_std # for dredge models, replacing the individual species above
  + Mean_Bird_Diversity.mean_std
  + GRIDMET_Precip_mm.mean_std
  + irrWater.mean_std
  + pdsi.mean_std
  + EVI.mean_std
  + VPD_kPa.mean_std
  + bite_mean_std
  + abund_mean_std
  + jrcStandingWater.mean_std
  + CApop_mean_std
  + Cattl_mean_std
  + Chick_mean_std
  + NLCD_majority2
  #+ xcoord_scale
  #+ ycoord_scale
  + offset(NumObservations.sum_std)
  ,
  data=WNV_dat,
  ziformula= ~0,
  family=tweedie(link="log"),
  na.action="na.fail"
)
summary(fit_tweedie_full)
check_collinearity(fit_tweedie_full)
#DHARMa diagnostics:
testDispersion(fit_tweedie_full)
simulationOutput <- simulateResiduals(fittedModel = fit_tweedie_full, plot = F)
residuals(simulationOutput)
plot(simulationOutput)
plotResiduals(simulationOutput, form = WNV_dat$trans_mean_std)

step_tweedie <- stepAIC(fit_tweedie_full, trace = TRUE, direction= "backward")
step_tweedie

dredge_results <- dredge(fit_tweedie_full, beta=c("sd"), fixed=c("cond(trans_mean_std)", "cond(offset(NumObservations.sum_std))"))
summary(get.models(dredge_results, 1)[[1]])

#residual spatial autocorrelation: no residual autocorrelation without x,y
testSpatialAutocorrelation(simulationOutput, x =  WNV_dat$xcoord.x, y = WNV_dat$ycoord.y)

#########################################################
# Best model:
# there IS residual spatial autocorrelation in the best fit model;
# so x, y are included here
#########################################################
fit_tweedie_best <- glmmTMB(
  MIR_Tar.mean ~
    trans_mean_std 
  + Mean_Bird_Diversity.mean_std
  + EVI.mean_std
  + VPD_kPa.mean_std
  + bite_mean_std
  + jrcStandingWater.mean_std
  + CApop_mean_std
  + Chick_mean_std
  + NLCD_majority2
  + xcoord_scale
  + ycoord_scale
  + offset(NumObservations.sum_std)
  ,
  data=WNV_dat,
  ziformula= ~0,
  family=tweedie(link="log"),
  na.action="na.fail"
)
summary(fit_tweedie_best)
check_collinearity(fit_tweedie_best)
#DHARMa diagnostics:
testDispersion(fit_tweedie_best)
simulationOutput <- simulateResiduals(fittedModel = fit_tweedie_best, plot = F)
residuals(simulationOutput)
plot(simulationOutput)
plotResiduals(simulationOutput, form = WNV_dat$trans_mean_std)

#residual spatial autocorrelation: yes, still residual autocorrelation without x,y, so add in x and y above
testSpatialAutocorrelation(simulationOutput, x =  WNV_dat$xcoord.x, y = WNV_dat$ycoord.y)

#for raw dependent: coef = 0.61 - for a 1 std increase in predicted temp-dependent transmission efficiency,
#we expect an increase of:
exp(0.61)-1 # 0.84, or 84%, in the minimum infection rate

#########################################################
#coefs:
#########################################################
ggcoef_model(fit_tweedie_best, exponentiate = T,
             show_p_values = FALSE,
             signif_stars = FALSE,
             significance = NULL,
             categorical_terms_pattern = "{level} (ref: {reference_level})",
             variable_labels = c(
               trans_mean_std = "Transmission Efficiency",
               Mean_Bird_Diversity.mean_std = "Bird Diversity",
               EVI.mean_std = "Enhanced Vegetation Index",
               VPD_kPa.mean_std = "Vapor Pressure Deficit",
               bite_mean_std = "Mos. Biting Rate",
               jrcStandingWater.mean_std = "Area Standing Water",
               CApop_mean_std = "Human Population Density",
               Chick_mean_std = "Density of Chickens",
               NLCD_majority2 = "Land Cover",
               xcoord_scale = "x coordinate",
               ycoord_scale = "y coordinate"
             )) +
  ggplot2::ggtitle("MIR ~ Transmission Efficiency") +
  theme(text=element_text(size=12)) +
  xlab(expression("exp(" ~ hat(beta) ~ ")"))

#########################################################
# Tables:
#########################################################
#SI models:
tab_model(fit_tweedie_base, fit_tweedie_full, file = "WNV_Transmission_GLM_base_and_full_R1.doc", pred.labels = c(
  "(Intercept)", "Transmission Efficiency", "American Crow Abundance", "American Robin Abundance", "House Finch Abundance",
  "House Sparrow Abundance", "Scrub Jay Abundance", "Bird Diversity", "Precip. (mm)", "Irrigated Ag. Area", "PDSI Drought", 
  "Enhanced Vegetation Index", "Vapor Pressure Deficit", "Mosquito Biting Rate", "Temp-Dependent Mosquito Abundance", 
  "Area of Standing Water", "Human Population Density", "Cattle Density", "Chicken Density", "Land Cover, Developed",
  "Land Cover, Natural", "Land Cover, Wetland"
), show.r2 = TRUE, show.p=FALSE, dv.labels = c("WNV MIR"), transform = "exp")
#best model
tab_model(fit_tweedie_best, file = "WNV_Transmission_GLM_best_R1.doc", pred.labels = c(
  "(Intercept)", "Transmission Efficiency", "Bird Diversity", "Enhanced Vegetation Index", "Vapor Pressure Deficit", 
  "Mosquito Biting Rate", "Area of Standing Water", "Human Population Density", "Chicken Density", "Land Cover, Developed", 
  "Land Cover, Natural", "Land Cover, Wetland", "x coordinate", "y coordinate"
), show.r2 = TRUE, show.p=FALSE, dv.labels = c("WNV MIR"), transform = "exp")

# full and best model from dredge:
tab_model(fit_tweedie_full, fit_tweedie_best, file = "WNV_Transmission_GLM_full_and_best_Dredge_R1.doc", pred.labels = c(
  "(Intercept)", "Transmission Efficiency", "Competent Bird Abundance", "Bird Diversity", "Precip. (mm)", "Irrigated Ag. Area", "PDSI Drought", 
  "Enhanced Vegetation Index", "Vapor Pressure Deficit", "Mosquito Biting Rate", "Temp-Dependent Mosquito Abundance", 
  "Area of Standing Water", "Human Population Density", "Cattle Density", "Chicken Density", "Land Cover, Developed",
  "Land Cover, Natural", "Land Cover, Wetland", "x coordinate", "y coordinate"
), show.r2 = TRUE, show.p=FALSE, dv.labels = c("WNV MIR"), transform = "exp")


#########################################################
# full Air Temp model: 
# dredge model selection was not feasibile here, 
# so model selection is based on backward stepwise selection
#########################################################
fit_tweedie_airT <- glmmTMB(
  MIR_Tar.mean ~
    airT_mean_std
  + I(airT_mean_std^2)
  + amecro_abundance.mean_std
  + amerob_abundance.mean_std
  + houfin_abundance.mean_std
  + houspa_abundance.mean_std
  + cowscj_abundance.mean_std
  + Mean_Bird_Diversity.mean_std
  + GRIDMET_Precip_mm.mean_std
  + irrWater.mean_std
  + pdsi.mean_std
  + EVI.mean_std
  + VPD_kPa.mean_std
  + jrcStandingWater.mean_std
  + CApop_mean_std
  + Cattl_mean_std
  + Chick_mean_std
  + NLCD_majority2
  + offset(NumObservations.sum_std)
  ,
  data=WNV_dat,
  ziformula= ~0,
  family=tweedie(link="log"),
  na.action="na.fail"
)
summary(fit_tweedie_airT)
check_collinearity(fit_tweedie_airT)
#DHARMa diagnostics:
testDispersion(fit_tweedie_airT)
simulationOutput <- simulateResiduals(fittedModel = fit_tweedie_airT, plot = F)
residuals(simulationOutput)
plot(simulationOutput)
plotResiduals(simulationOutput, form = WNV_dat$airT_mean_std)

#residual spatial autocorrelation: no residual autocorrelation without x,y
testSpatialAutocorrelation(simulationOutput, x =  WNV_dat$xcoord.x, y = WNV_dat$ycoord.y)

step_tweedie <- stepAIC(fit_tweedie_airT, trace = TRUE, direction= "backward")
step_tweedie

###############################################
# best model:
###############################################
fit_tweedie_airT_best <- glmmTMB(
  MIR_Tar.mean ~
  #MIR_Tar.mean_std ~ #standardize outcome, for figure
    airT_mean_std
  + I(airT_mean_std^2)
  + amerob_abundance.mean_std
  + houspa_abundance.mean_std
  + cowscj_abundance.mean_std
  + Mean_Bird_Diversity.mean_std
  + EVI.mean_std
  + CApop_mean_std
  + Chick_mean_std
  + offset(NumObservations.sum_std)
  ,
  data=WNV_dat,
  ziformula= ~0,
  family=tweedie(link="log"),
  na.action="na.fail"
)
summary(fit_tweedie_airT_best)
check_collinearity(fit_tweedie_airT_best)
#DHARMa diagnostics:
testDispersion(fit_tweedie_airT_best)
simulationOutput <- simulateResiduals(fittedModel = fit_tweedie_airT_best, plot = F)
residuals(simulationOutput)
plot(simulationOutput)
plotResiduals(simulationOutput, form = WNV_dat$airT_mean_std)

#residual spatial autocorrelation: no residual autocorrelation without x,y
testSpatialAutocorrelation(simulationOutput, x =  WNV_dat$xcoord.x, y = WNV_dat$ycoord.y)

#Find inflection point:
summary(fit_tweedie_airT_best)
fit <- fit_tweedie_airT_best$fit
coefs <- fit$par
coefs_lin <- coefs[2]
head(coefs_lin)
coefs_sq <- coefs[3]
coefs_mat <- cbind(coefs_lin, coefs_sq)
head(coefs_mat)
inf.point <- -coefs_mat[,1]/(2*coefs_mat[,2])
inf.point

#for raw MIR outcome
inf_WNV_dat <- WNV_dat[which(WNV_dat$airT_mean_std > 0.78 & WNV_dat$airT_mean_std < 0.88),]
head(inf_WNV_dat)
nrow(inf_WNV_dat)
min(inf_WNV_dat$airT_mean)
max(inf_WNV_dat$airT_mean)
mean(inf_WNV_dat$airT_mean)
# models with airT predict a peak of MIR at ~ 25.1 - 25.2 C (in the middle of the predicted optimal range in Shocket et al. 2020)
#full range in data:
min(WNV_dat$airT_mean)
max(WNV_dat$airT_mean)

#for scaled MIR outcome
inf_WNV_dat <- WNV_dat[which(WNV_dat$airT_mean_std > 0.25 & WNV_dat$airT_mean_std < 0.32),]
head(inf_WNV_dat)
nrow(inf_WNV_dat)
min(inf_WNV_dat$airT_mean)
max(inf_WNV_dat$airT_mean)
mean(inf_WNV_dat$airT_mean)
# models with airT predict a peak of MIR at ~ 24.6 - 24.7 C (in the middle of the predicted optimal range in Shocket et al. 2020)
#full range in data:
min(WNV_dat$airT_mean)
max(WNV_dat$airT_mean)

visreg(fit_tweedie_airT_best, "airT_mean_std", xlab="Air Temp, std", ylab="WNV MIR", main="WNV ~ Air Temperature", cex=1.5,
       #fill=list(col="maroon", alpha=0.05),
       line=list(lty=1, col="darkred"),
       points=list(col="red"))

#coefs:
ggcoef_model(fit_tweedie_airT_best, exponentiate = T,
             show_p_values = FALSE,
             signif_stars = FALSE,
             significance = NULL,
             variable_labels = c(
               airT_mean_std = "Air Temperature",
               "I(airT_mean_std^2)" = "Air Temperature squared",
               amerob_abundance.mean_std = "Avg. Robin Abundance",
               houspa_abundance.mean_std = "Avg. House Sparrow Abundance",
               cowscj_abundance.mean_std = "Avg. Jay Abundance",
               Mean_Bird_Diversity.mean_std = "Bird Diversity",
               EVI.mean_std = "Enhanced Vegetation Index",
               CApop_mean_std = "Human Population Density",
               Chick_mean_std = "Density of Chickens"
             )) +
  ggplot2::ggtitle("MIR ~ Air Temperature") +
  theme(text=element_text(size=12)) +
  xlab(expression("exp(" ~ hat(beta) ~ ")"))


# Tables:
#airT raw:
tab_model(fit_tweedie_airT_best, file = "WNV_airT_Transmission_GLM_best_R1.doc", pred.labels = c(
  "(Intercept)", "Air Temp (C)", "Air Temp sq. (C)", "American Robin Abundance", "House Sparrow Abundance", 
  "Scrub Jay Abundance", "Bird Diversity", "Enhanced Vegetation Index", "Human Population Density", "Chicken Density"
), show.r2 = TRUE, show.p=FALSE, dv.labels = c("WNV MIR"), transform = "exp")
#airT std:
tab_model(fit_tweedie_airT_best, file = "WNV_airT_Transmission_zscore_GLM_best_R1.doc", pred.labels = c(
  "(Intercept)", "Air Temp (C)", "Air Temp sq. (C)", "American Robin Abundance", "House Sparrow Abundance", 
  "Scrub Jay Abundance", "Bird Diversity", "Enhanced Vegetation Index", "Human Population Density", "Chicken Density"
), show.r2 = TRUE, show.p=FALSE, dv.labels = c("WNV MIR"), transform = "exp")

