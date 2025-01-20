#!/usr/bin/env Rscript

library(tidyverse)
library(lubridate)
library(rsofun)

# load script arguments
args <- commandArgs(trailingOnly = TRUE)

# load data
if (is.na(args[1])){
  load("data-raw/CH-LAE_forcing.rda")
} else {
  load(args[1])
}

# sitename
sitename <- "CH-Lae"

# Take only year 2004 to 2014, corresponding to subset of data for site CH-Lae
siteinfo <- tibble(
  lon = 8.365,
  lat = 47.47808,
  elv = 700,
  year_start = 2004,
  year_end = 2014,
  c4 = FALSE,
  igbp_land_use = "Mixed Forests",
  plant_functional_type = "Broadleaf trees")

siteinfo <- siteinfo %>% 
  dplyr::mutate(date_start = lubridate::ymd(paste0(year_start, "-01-01"))) %>%
  dplyr::mutate(date_end = lubridate::ymd(paste0(year_end, "-12-31")))

# load model parameters (valid ones)
params_siml_gs_leuning <- tibble(
  spinupyears = 250,
  recycle = 1,
  firstyeartrend = 2009,
  nyeartrend = 1,
  steps_per_day = 24,
  do_U_shaped_mortality = TRUE,
  do_closedN_run = TRUE,
  method_photosynth = "gs_leuning",
  method_mortality = "dbh"
)

params_siml_pmodel <- params_siml_gs_leuning
params_siml_pmodel$method_photosynth <- "pmodel"
params_siml_pmodel$steps_per_day <- 1

params_tile <- tibble(
  soiltype = 3,
  FLDCAP = 0.4,
  WILTPT = 0.05,
  K1 = 2.0,
  K2 = 0.05,
  K_nitrogen = 8.0,
  MLmixRatio = 0.8,
  etaN = 0.025,
  LMAmin = 0.02,
  fsc_fine = 1.0,
  fsc_wood = 0.0,
  GR_factor = 0.33,
  l_fract = 0.0,
  retransN = 0.0,
  f_initialBSW = 0.2,
  f_N_add = 0.02,
  tf_base = 1,
  par_mort = 1,
  par_mort_under = 1
)

params_species <- tibble(
  lifeform      = rep(1,4),
  phenotype     = c(0,1,1,1),
  pt            = rep(0,4),
  # Root parameters
  alpha_FR      = rep(1.2,4),
  rho_FR        = rep(200,4),
  root_r        = rep(2.9E-4,4),
  root_zeta     = rep(0.29,4),
  Kw_root       = rep(3.5e-09,4),
  leaf_size     = rep(0.04,4),
  # Photosynthesis parameters
  Vmax          = rep(35.0E-6,4),
  Vannual       = rep(1.2,4),
  wet_leaf_dreg = rep(0.3,4),
  m_cond        = rep(7.0,4),
  alpha_phot    = rep(0.06,4),
  gamma_L       = rep(0.02,4),
  gamma_LN      = rep(70.5,4),
  gamma_SW      = rep(0.08,4),
  gamma_FR      = rep(12.0,4),
  tc_crit       = rep(283.16,4),
  tc_crit_on    = rep(280.16,4),
  gdd_crit      = rep(280.0,4),
  betaON        = rep(0,2,4),
  betaOFF       = rep(0,1,4),
  # Allometry parameters
  alphaHT       = rep(36,4),
  thetaHT       = rep(0.5,4),
  alphaCA       = rep(150,4),
  thetaCA       = rep(1.5,4),
  alphaBM       = rep(5200,4),
  thetaBM       = c(2.36,2.30,2.54,2.30),
  # Reproduction parameters
  seedlingsize  = rep(0.05,4),
  maturalage    = rep(5,4),
  v_seed        = rep(0.1,4),
  # Mortality parameters
  mortrate_d_c  = rep(0.01,4),
  mortrate_d_u  = rep(0.075,4),
  # Leaf parameters
  LMA           = c(0.05,0.17,0.11,0.1),
  leafLS        = rep(1,4),
  LNbase        = rep(0.8E-3,4),
  CNleafsupport = rep(80,4),
  rho_wood      = c(590,370,350,300),
  taperfactor   = rep(0.75,4),
  lAImax        = rep(3.5,4),
  tauNSC        = rep(3,4),
  fNSNmax       = rep(5,4),
  phiCSA        = rep(0.25E-4,4),
  # C/N ratios for plant pools
  CNleaf0       = rep(25,4),
  CNsw0         = rep(350,4),
  CNwood0       = rep(350,4),
  CNroot0       = rep(40,4),
  CNseed0       = rep(20,4),
  Nfixrate0     = rep(0,4),
  NfixCost0     = rep(12,4),
  internal_gap_frac      = rep(0.1,4),
  # calibratable params
  kphio         = rep(0.05,4),
  phiRL         = rep(3.5,4),
  LAI_light     = rep(3.5,4)
)

init_cohort <- tibble(
  init_cohort_species = 2,    # indicates sps # 1 - Fagus sylvatica
  init_cohort_nindivs = 0.05, # initial individual density, individual/m2 ! 1 indiv/m2 = 10.000 indiv/ha
  init_cohort_bl      = 0.0,  # initial biomass of leaves, kg C/individual
  init_cohort_br      = 0.0,  # initial biomass of fine roots, kg C/individual
  init_cohort_bsw     = 0.05, # initial biomass of sapwood, kg C/individual
  init_cohort_bHW     = 0.0,  # initial biomass of heartwood, kg C/tree
  init_cohort_seedC   = 0.0,  # initial biomass of seeds, kg C/individual
  init_cohort_nsc     = 0.05  # initial non-structural biomass
)

init_soil <- tibble( #list
  init_fast_soil_C    = 0.01,
  init_slow_soil_C    = 0.001,
  init_Nmineral       = 0.015,
  N_input             = 0.0008
)

rh_to_vpd <- function(temp, # Air temperature (deg C)
                      rh    # Relative humidity (< 1)
) {
  esat <- 611.0 * exp( (17.27 * temp)/(temp + 237.3) )
  
  return(esat * (1.0 - rh)) # VPD (Pa)
}

rad_to_ppfd <- function(rad   # Downwelling radiation (W m-2)
) {
  kcFCE <- 2.04 # from flux to energy conversion, umol/J (Meek et al., 1984)
  
  return(rad * kcFCE * 1.0e-6)  # PPFD (mol m-2 s-1)
}

build_forcing <- function(forcing_data, hourly) {
  if (hourly)
    groups <- forcing_data %>% dplyr::group_by(
      lubridate::month(datehour),
      lubridate::day(datehour),
      lubridate::hour(datehour))
  else
    groups <- forcing_data %>% dplyr::group_by(
      lubridate::month(datehour),
      lubridate::day(datehour))
  forcing <- groups %>%
    summarise_at(vars(1:13), list(~mean(., na.rm = TRUE))) %>%
    rename(month=`lubridate::month(datehour)`,day=`lubridate::day(datehour)`) %>%
    ungroup() %>%
    rename(year=YEAR, hod=HOUR, rad=Swdown, temp=TEMP, rh=RH,
           rain=RAIN, wind=WIND, patm=PRESSURE, co2=aCO2_AW) %>%
    mutate(date = make_date(year,month,day),
           vpd = rh_to_vpd(temp, rh / 100.0),
           ppfd = rad_to_ppfd(rad)) %>%
    select(date,hod,temp,rain,vpd,ppfd,patm,wind,co2,)
  return(forcing)
}

build_driver <- function(params_siml, forcing) {
  drivers <- tibble(
    sitename,
    site_info = list(tibble(siteinfo)),
    params_siml = list(tibble(params_siml)),
    params_tile = list(tibble(params_tile)),
    params_species = list(tibble(params_species)),
    init_cohort = list(tibble(init_cohort)),
    init_soil = list(tibble(init_soil)),
    forcing = list(tibble(forcing))
  )
  return(drivers)
}

#---- gs leuning formatting -----
forcing_gs_leuning <- build_forcing(forcingLAE, TRUE)

biomee_gs_leuning_drivers <- build_driver(params_siml_gs_leuning, forcing_gs_leuning)

save(biomee_gs_leuning_drivers,
     file ="data/biomee_gs_leuning_drivers.rda",
     compress = "xz")

#---- p-model formatting -----
forcing_pmodel <- build_forcing(forcingLAE, FALSE)

biomee_p_model_drivers <- build_driver(params_siml_pmodel, forcing_pmodel)

save(biomee_p_model_drivers,
     file ="data/biomee_p_model_drivers.rda",
     compress = "xz")
