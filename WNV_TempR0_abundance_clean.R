#clear memory:
rm(list=ls())

require('glmmTMB')
library('MuMIn')
require('tidyverse')
require('MASS')
require('performance')
require('DHARMa')
library('sjPlot')
library('sjmisc')
library('sjlabelled')
library('ggstats')

#########################################################
### Pull in mosquito data: *write new file with just vars used
#########################################################
abund_dat <- read.csv("Abundance_data.csv", header=T)
head(abund_dat)
lapply(abund_dat, class)
abund_dat$NLCD_majority2 <- as.factor(abund_dat$NLCD_majority2)

#########################################################
#Base model:
#########################################################
fit_tweedie_base <- glmmTMB(
  mosPerTrapNight.mean ~
    abund_mean_std 
  ,
  data=abund_dat,
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
plotResiduals(simulationOutput, form = abund_dat$abund_mean_std)

#residual spatial autocorrelation: yes, still residual autocorrelation
testSpatialAutocorrelation(simulationOutput, x =  abund_dat$xcoord.x, y = abund_dat$ycoord.y)

#########################################################
#Full model:
#here there IS significant spatial autocorrelation without x,y in full model, 
#so add them in for model selection:
#########################################################
fit_tweedie_full <- glmmTMB(
  mosPerTrapNight.mean ~
    abund_mean_std 
  + competent_bird_abundance.mean_std 
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
  + xcoord_scale #include x and y in model selection
  + ycoord_scale
  ,
  data=abund_dat,
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
plotResiduals(simulationOutput, form = abund_dat$abund_mean_std)

step_tweedie <- stepAIC(fit_tweedie_full, trace = TRUE, direction= "backward")
step_tweedie

dredge_results <- dredge(fit_tweedie_full, beta=c("sd"), fixed=c("cond(abund_mean_std)"))
summary(get.models(dredge_results, 1)[[1]])

#residual spatial autocorrelation: there is spatial autocorrelation without x,y, but not with x,y, so include in model selection
testSpatialAutocorrelation(simulationOutput, x =  abund_dat$xcoord.x, y = abund_dat$ycoord.y)

#########################################################
#Best model:
#########################################################
fit_tweedie_best <- glmmTMB(
  mosPerTrapNight.mean ~
    abund_mean_std 
  + competent_bird_abundance.mean_std
  + Mean_Bird_Diversity.mean_std
  + GRIDMET_Precip_mm.mean_std
  + irrWater.mean_std
  + EVI.mean_std
  + VPD_kPa.mean_std
  + jrcStandingWater.mean_std
  + Cattl_mean_std
  + NLCD_majority2
  + xcoord_scale
  + ycoord_scale
  ,
  data=abund_dat,
  ziformula= ~0,
  family=tweedie(link="log"),
  na.action=na.exclude
)
summary(fit_tweedie_best)
check_collinearity(fit_tweedie_best)
#DHARMa diagnostics:
testDispersion(fit_tweedie_best)
simulationOutput <- simulateResiduals(fittedModel = fit_tweedie_best, plot = F)
residuals(simulationOutput)
plot(simulationOutput)
plotResiduals(simulationOutput, form = abund_dat$abund_mean_std)

#residual spatial autocorrelation: no residual autocorrelation
testSpatialAutocorrelation(simulationOutput, x =  abund_dat$xcoord.x, y = abund_dat$ycoord.y)

#for raw dependent: coef = 0.35 - for a 1 std increase in predicted temp-dependent mosquito abundance,
#we expect an increase of:
exp(0.35)-1 # = 0.42, or 42%, in the number of mosquitos per trap night

#########################################################
#coef plots:
#########################################################
ggcoef_model(fit_tweedie_best, exponentiate = T, 
             show_p_values = FALSE,
             signif_stars = FALSE,
             significance = NULL,
             categorical_terms_pattern = "{level} (ref: {reference_level})",
             variable_labels = c(
               abund_mean_std = "Temp-Dependent Mos. Abund.",
               competent_bird_abundance.mean_std = "Competent Bird Host Abund.",
               Mean_Bird_Diversity.mean_std = "Bird Diversity",
               GRIDMET_Precip_mm.mean_std = "Precip. (mm)",
               irrWater.mean_std = "Area Irrigated Ag.",
               EVI.mean_std = "Enhanced Vegetation Index",
               VPD_kPa.mean_std = "Vapor Pressure Deficit",
               jrcStandingWater.mean_std = "Area Standing Water",
               Cattl_mean_std = "Density of Cattle",
               NLCD_majority2 = "Land Cover",
               xcoord_scale = "x coordinate",
               ycoord_scale = "y coordinate"
             )) +
  ggplot2::ggtitle("Mos. per Trap Night ~ \nTemp-Dependent Mos. Abund.") + 
  theme(text=element_text(size=12)) +
  xlab(expression("exp(" ~ hat(beta) ~ ")"))

#########################################################
# Tables:
#########################################################
#SI models:
tab_model(fit_tweedie_base, fit_tweedie_full, file = "WNV_Abundance_GLM_base_and_full_R1.doc", pred.labels = c(
  "(Intercept)", "Temp-Dependent Mos. Abundance", "Competent Bird Host Abundance", "Bird Diversity", "Precip. (mm)", "Irrigated Ag. Area", "PDSI Drought", 
  "Enhanced Vegetation Index", "Vapor Pressure Deficit", "Area of Standing Water", "Human Population Density", "Cattle Density", "Chicken Density",
  "Land Cover, Developed", "Land Cover, Natural", "Land Cover, Wetland", "x coordinate", "y coordinate"
), show.r2 = TRUE, show.p = FALSE, dv.labels = c("Mos. per Trap Night"), transform = "exp")
#best model
tab_model(fit_tweedie_best, file = "WNV_Abundance_GLM_best_R1.doc", pred.labels = c(
  "(Intercept)", "Temp-Dependent Mos. Abundance", "Competent Bird Host Abundance", "Bird Diversity", "Precip. (mm)", "Irrigated Ag. Area", 
  "Enhanced Vegetation Index", "Vapor Pressure Deficit", "Area of Standing Water", "Cattle Density", "Land Cover, Developed", "Land Cover, Natural", "Land Cover, Wetland",
  "x coordinate", "y coordinate"
), show.r2 = TRUE, show.p = FALSE, dv.labels = c("Mos. per Trap Night"), transform = "exp")

#dredge full and best table:
tab_model(fit_tweedie_full, fit_tweedie_best, file = "WNV_Abundance_GLM_full_and_best_Dredge_R1.doc", pred.labels = c(
  "(Intercept)", "Temp-Dependent Mos. Abundance", "Competent Bird Host Abundance", "Bird Diversity", "Precip. (mm)", "Irrigated Ag. Area", "PDSI Drought", 
  "Enhanced Vegetation Index", "Vapor Pressure Deficit", "Area of Standing Water", "Human Population Density", "Cattle Density", "Chicken Density",
  "Land Cover, Developed", "Land Cover, Natural", "Land Cover, Wetland", "x coordinate", "y coordinate"
), show.r2 = TRUE, show.p = FALSE, dv.labels = c("Mos. per Trap Night"), transform = "exp")

