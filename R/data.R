#' rsofun P-model driver data
#'
#' Small dataset representing the driver to run the P-model at the FR-Pue site. 
#' It can also be used together with daily GPP flux time series data from CH-LAE 
#' (\code{\link{p_model_validation}}) to optimize model parameters.
#' To optimize model parameters to leaf traits data use the datasets \code{\link{p_model_drivers_vcmax25}} and \code{\link{p_model_validation_vcmax25}}.
#'
#' @format A tibble of driver data:
#' \describe{
#'   \item{sitename}{A character string containing the site name.}
#'   \item{forcing}{A tibble of a time series of forcing climate data, including
#'   the following data:
#'     \describe{
#'       \item{date}{Date of the observation in YYYY-MM-DD format.}
#'       \item{temp}{Daytime average air temperature in \eqn{^\circ}C.}
#'       \item{vpd}{Daytime average vapour pressure deficit in Pa.}
#'       \item{ppfd}{Photosynthetic photon flux density (PPFD) in
#'       mol m\eqn{^{-2}} s\eqn{^{-1}}. If all values are NA, it indicates that
#'       PPFD should be calculated by the SPLASH model.}
#'       \item{netrad}{Net radiation in W m\eqn{^{-2}}. This is currently
#'       ignored as a model forcing.}
#'       \item{patm}{Atmospheric pressure in Pa.}
#'       \item{snow}{Snow in water equivalents mm s\eqn{^{-1}}.}
#'       \item{rain}{Rain as precipitation in liquid form in mm s\eqn{^{-1}}.}
#'       \item{tmin}{Daily minimum air temperature in \eqn{^\circ}C.}
#'       \item{tmax}{Daily maximum air temperature in \eqn{^\circ}C.}
#'       \item{fapar}{Fraction of photosynthetic active radiation (fAPAR), taking
#'      values between 0 and 1.}
#'       \item{co2}{Atmospheric CO\eqn{_2} concentration.}
#'       \item{ccov}{Cloud coverage in \%. This is only used when either PPFD or
#'       net radiation are not prescribed.}
#'       }
#'   }
#'   \item{params_siml}{A tibble of simulation parameters, including
#'   the following data:
#'     \describe{
#'       \item{spinup}{A logical value indicating whether this simulation does spin-up.}
#'       \item{spinupyears}{Number of spin-up years.}
#'       \item{recycle}{Number of first N years of forcing data.frame that are recycled for spin-up.}
#'       \item{outdt}{An integer indicating the output periodicity.}
#'       \item{ltre}{A logical value, \code{TRUE} if evergreen tree.}
#'       \item{ltne}{A logical value, \code{TRUE} if evergreen tree and N-fixing.}
#'       \item{ltrd}{A logical value, \code{TRUE} if deciduous tree.}
#'       \item{ltnd}{A logical value, \code{TRUE} if deciduous tree and N-fixing.}
#'       \item{lgr3}{A logical value, \code{TRUE} if grass with C3 photosynthetic pathway.}
#'       \item{lgn3}{A logical value, \code{TRUE} if grass with C3 photosynthetic
#'       pathway and N-fixing.}
#'       \item{lgr4}{A logical value, \code{TRUE} if grass with C4 photosynthetic pathway.}
#'     }
#'   }
#'   \item{site_info}{A tibble containing site meta information.
#' This data structure can be freely used for documenting the dataset, but must include at least the following data:
#'     \describe{
#'       \item{lon}{Longitude of the site location in degrees east.}
#'       \item{lat}{Latitude of the site location in degrees north.}
#'       \item{elv}{Elevation of the site location, in meters above sea level.}
#'       \item{whc}{A numeric value for the rooting zone water holding capacity (in mm)}
#'     }
#'   }
#' }
#' 
#' @source Pastorello, G., Trotta, C., Canfora, E. et al. 
#' The FLUXNET2015 dataset and the ONEFlux processing pipeline for eddy covariance data. 
#' Sci Data 7, 225 (2020). https://doi.org/10.1038/s41597-020-0534-3
#' 
#' University of East Anglia Climatic Research Unit; Harris, I.C.; Jones, P.D.; Osborn, T. (2021): 
#' CRU TS4.05: Climatic Research Unit (CRU) Time-Series (TS) version 4.05 of high-resolution 
#' gridded data of month-by-month variation in climate (Jan. 1901- Dec. 2020). 
#' NERC EDS Centre for Environmental Data Analysis, date of citation. 
#' https://catalogue.ceda.ac.uk/uuid/c26a65020a5e4b80b20018f148556681
#' 
#' Weedon, G. P., G. Balsamo, N. Bellouin,S. Gomes, M. J. Best, and P. Viterbo(2014), 
#' The WFDEI meteorologicalforcing data set: WATCH Forcing Datamethodology applied 
#' to ERA-Interimreanalysis data,
#' Water Resour. Res.,50,7505–7514, doi:10.1002/2014WR015638.
#' 
#' Fick, S.E. and R.J. Hijmans, 2017. WorldClim 2: new 1km spatial resolution climate 
#' surfaces for global land areas. International Journal of Climatology 37 (12): 4302-4315.
"p_model_drivers"

#' rsofun P-model GPP validation data
#'
#' Small example dataset of target observations (daily GPP flux data) to optimize 
#' model parameters with the function \code{\link{calib_sofun}}
#'
#' @format A tibble of validation data:
#' \describe{
#'   \item{sitename}{A character string containing the site name (e.g. 'FR-Pue').}
#'   \item{data}{A tibble [ 2,920 x 3 ] with time series for the following variables:
#'     \describe{
#'       \item{date}{Date vector with format YYYY-MM-DD.}
#'       \item{gpp}{The observed Gross Primary Productivity (GPP) for each time stamp 
#'       (in gC m\eqn{^{-2}} d\eqn{^{-1}}).}
#'       \item{gpp_unc}{The uncertainty of the GPP (in gC m\eqn{^{-2}} d\eqn{^{-1}}).}
#'     }
#'   }
#' }
#' @examples require(ggplot2); require(tidyr)
#' p_model_validation %>% tidyr::unnest(data) 
#' 
#' @source Pastorello, G., Trotta, C., Canfora, E. et al. 
#' The FLUXNET2015 dataset and the ONEFlux processing pipeline for eddy covariance data. 
#' Sci Data 7, 225 (2020). https://doi.org/10.1038/s41597-020-0534-3
"p_model_validation"

#' rsofun P-model driver data (for leaf traits)
#'
#' Small dataset representing the driver to run the P-model at four separate sites. 
#' It can also be used together with leaf traits data from these four sites
#' (\code{\link{p_model_validation_vcmax25}}) to optimize model parameters.
#' To optimize model parameters to GPP flux data use the datasets \code{\link{p_model_drivers}} and \code{\link{p_model_validation}}.
#' 
#' @format See \code{\link{p_model_drivers}}
#' 
#' @source Atkin, O. K., Bloomfield, K. J., Reich, P. B., Tjoelker, M. G., Asner, G. P., Bonal, D., et al. (2015). 
#' Global variability in leaf respiration in relation to climate, plant functional types and leaf traits. 
#' New Phytol. 206 (2), 614–636. doi:10.1111/nph.13253
#' 
#' University of East Anglia Climatic Research Unit; Harris, I.C.; Jones, P.D.; Osborn, T. (2021): 
#' CRU TS4.05: Climatic Research Unit (CRU) Time-Series (TS) version 4.05 of high-resolution 
#' gridded data of month-by-month variation in climate (Jan. 1901- Dec. 2020). 
#' NERC EDS Centre for Environmental Data Analysis, date of citation. 
#' https://catalogue.ceda.ac.uk/uuid/c26a65020a5e4b80b20018f148556681
#' 
#' Weedon, G. P., G. Balsamo, N. Bellouin,S. Gomes, M. J. Best, and P. Viterbo(2014), 
#' The WFDEI meteorologicalforcing data set: WATCH Forcing Datamethodology applied 
#' to ERA-Interimreanalysis data,
#' Water Resour. Res.,50,7505–7514, doi:10.1002/2014WR015638.
#' 
#' Fick, S.E. and R.J. Hijmans, 2017. 
#' WorldClim 2: new 1km spatial resolution climate surfaces for global land areas. 
#' International Journal of Climatology 37 (12): 4302-4315.
#' 
#' C.D. Keeling, R.B. Bacastow, A.E. Bainbridge, C.A. Ekdahl, P.R. Guenther, and L.S. Waterman, (1976), 
#' Atmospheric carbon dioxide variations at Mauna Loa Observatory, Hawaii, Tellus, vol. 28, 538-551
"p_model_drivers_vcmax25"

#' rsofun P-model Vcmax25 validation data
#'
#' Small example dataset of target observations (leaf trait data) to optimize 
#' model parameters with the function \code{\link{calib_sofun}}
#'
#' @format A tibble of validation data:
#' \describe{
#'   \item{sitename}{A character string containing the site names (e.g. 'Reichetal_Colorado').}
#'   \item{data}{A tibble [ 1 x 2 ] with observations for the following variables:
#'     \describe{
#'       \item{vcmax25}{The observed maximum rate of carboxylation (Vcmax), normalised 
#'     to 25\eqn{^o} C (in mol C m\eqn{^{-2}} d\eqn{^{-1}}), aggregated over different plant species
#'     in each site.}
#'       \item{vcmax25_unc}{The uncertainty of the Vcmax25 (in mol C m\eqn{^{-2}} d\eqn{^{-1}}),
#'       calculated as the standard deviation among Vcmax25 observations for
#'       several species per site or as the total standard deviation across sites for
#'       single-plant-species sites.}
#'     }
#'   }
#' }
#' @examples require(ggplot2); require(tidyr)
#' p_model_validation_vcmax25 %>% tidyr::unnest(data) 
#' 
#' @source Atkin, O. K., Bloomfield, K. J., Reich, P. B., Tjoelker, M. G., Asner, G. P., Bonal, D., et al. (2015). 
#' Global variability in leaf respiration in relation to climate, plant functional types and leaf traits. 
#' New Phytol. 206 (2), 614–636. doi:10.1111/nph.13253
"p_model_validation_vcmax25"

#' rsofun P-model output data
#'
#' Example output dataset from a p-model run using \code{\link{p_model_drivers}}
#' See \code{\link{run_pmodel_f_bysite}} for a detailed 
#' description of the outputs.
"p_model_output"

#' rsofun P-model output data (using vcmax25 drivers)
#'
#' Example output dataset from a p-model run using \code{\link{p_model_drivers_vcmax25}}
#' See \code{\link{run_pmodel_f_bysite}} for a detailed 
#' description of the outputs.
"p_model_output_vcmax25"

#' rsofun BiomeE driver data (Leuning photosynthesis model)
#' 
#' Small dataset representing the driver to run the BiomeE-model at the CH-LAE site
#' using the Leuning photosynthesis specification (and half-hourly time step)
#' It can also be used together with leaf trait data from CH-LAE (\code{\link{biomee_validation}}) 
#' to optimize model parameters.
#'
#' @format A tibble of driver data.
#' \describe{
#'   \item{sitename}{Site name}
#'   \item{params_siml}{Simulation parameters as a data.frame, including
#'   the following data:
#'     \describe{
#'       \item{spinup}{Flag indicating whether this simulation does spin-up (deprecated).}
#'       \item{spinupyears}{Number of spin-up years. Set to 0 for no spinup.}
#'       \item{recycle}{Number of first N years of forcing data.frame that are recycled for spin-up.}
#'       \item{firstyeartrend}{Year of first transient year (AD) (optional). Is only used to set years in output data frames. Defaults to 0 if not provided.}
#'       \item{nyeartrend}{Number of transient years (optional). Determines the length of simulation output after spin-up. Defaults to number of years contained in the forcing data. (If longer than forcing data, last year of forcing is repeated until the end (spin-down).)}
#'       \item{steps_per_day}{Time resolution of the forcing (day-1).}
#'       \item{do_U_shaped_mortality}{Flag indicating whether U-shaped
#'         mortality is used.}
#'       \item{do_closedN_run}{Flag indicating whether doing N closed
#'         runs to recover N balance enforcing 0.2 kg N m-2 in the inorganic N pool.}
#'       \item{code_method_photosynth}{String specifying the method of photosynthesis
#'         used in the model, either "pmodel" or "gs_leuning".document()}
#'       \item{code_method_mortality}{String indicating the type of mortality in the
#'         model. One of the following: "dbh" is size-dependent mortality, "const_selfthin"
#'         is constant self thinning (in development), "cstarvation" is carbon starvation, and
#'         "growthrate" is growth rate dependent mortality.}
#'     }}
#'   \item{site_info}{Site meta info in a data.frame.
#' This data structure can be freely used for documenting the dataset, but must include at least the following data:
#'     \describe{
#'       \item{lon}{Longitude of the site location in degrees east.}
#'       \item{lat}{Latitude of the site location in degrees north.}
#'       \item{elv}{Elevation of the site location, in meters above sea level.}
#'     }}
#'   \item{forcing}{Forcing data.frame used as input
#'     \describe{
#'       \item{ppfd}{Photosynthetic photon flux density (mol s-1 m-2)}
#'       \item{tair}{Air temperature (deg C)}
#'       \item{vpd}{Vapor pressure deficit (Pa)}
#'       \item{rain}{Precipitation (kgH2O m-2 s-1 == mm s-1)}
#'       \item{wind}{Wind velocity (m s-1)}
#'       \item{pair}{Atmospheric pressure (pa)}
#'       \item{co2}{CO2 atmospheric concentration (ppm)}
#'     }}
#'   \item{params_tile}{Tile-level model parameters, into a single row data.frame, including
#'   the following data:
#'     \describe{
#'       \item{soiltype}{Integer indicating the type of soil: Sand = 1, LoamySand = 2,
#'         SandyLoam = 3, SiltLoam = 4, FrittedClay = 5, Loam = 6, Clay = 7.}
#'       \item{FLDCAP}{Field capacity (vol/vol). Water remaining in a soil after it
#'         has been thoroughly saturated and allowed to drain freely.}
#'       \item{WILTPT}{Wilting point (vol/vol). Water content of a soil at which
#'       plants wilt and fail to recover.}
#'       \item{K1}{Fast soil C decomposition rate (year\eqn{^{-1}}).}
#'       \item{K2}{Slow soil C decomposition rate (year\eqn{^{-1}}).}
#'       \item{K_nitrogen}{Mineral Nitrogen turnover rate (year\eqn{^{-1}}).}
#'       \item{MLmixRatio}{Ratio of C and N returned to litters from microbes.}
#'       \item{etaN}{N loss rate through runoff (organic and mineral) (year\eqn{^{-1}}).}
#'       \item{LMAmin}{Minimum LMA, leaf mass per unit area, kg C m\eqn{^{-2}}.}
#'       \item{fsc_fine}{Fraction of fast turnover carbon in fine biomass.}
#'       \item{fsc_wood}{Fraction of fast turnover carbon in wood biomass.}
#'       \item{GR_factor}{Growth respiration factor.}
#'       \item{l_fract}{Fraction of the carbon retained after leaf drop.}
#'       \item{retransN}{Retranslocation coefficient of nitrogen.}
#'       \item{f_initialBSW}{Coefficient for setting up initial sapwood.}
#'       \item{f_N_add}{Re-fill of N for sapwood.}
#'       \item{tf_base}{Calibratable scalar for respiration, used to increase LUE
#'         levels.}
#'       \item{par_mort}{Canopy mortality parameter.}
#'       \item{par_mort_under}{Parameter for understory mortality.}
#'     }}
#'   \item{params_species}{A data.frame containing species-specific model parameters,
#'       with one species per row, including the following data:
#'     \describe{
#'       \item{lifeform}{Integer set to 0 for grasses and 1 for trees.}
#'       \item{phenotype}{Integer set to 0 for deciduous and 1 for evergreen.}
#'       \item{pt}{Integer indicating the type of plant according to photosynthesis:
#'         0 for C3; 1 for C4}
#'       \item{alpha_FR}{Fine root turnover rate (year\eqn{^{-1}}).}
#'       \item{rho_FR}{Material density of fine roots (kg C m\eqn{^{-3}}).}
#'       \item{root_r}{Radius of the fine roots, in m.}
#'       \item{root_zeta}{e-folding parameter of root vertical distribution, in m.}
#'       \item{Kw_root}{Fine root water conductivity (mol m\eqn{^{-2}}
#'         s\eqn{^{-1}} MPa\eqn{^{-1}}).}
#'       \item{leaf_size}{Characteristic leaf size.}
#'       \item{Vmax}{Max RuBisCo rate, in mol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'       \item{Vannual}{Annual productivity per unit area at full sun (kg C
#'         m\eqn{^{-2}} year\eqn{^{-2}}).}
#'       \item{wet_leaf_dreg}{Wet leaf photosynthesis down-regulation.}
#'       \item{m_cond}{Factor of stomatal conductance.}
#'       \item{alpha_phot}{Photosynthesis efficiency.}
#'       \item{gamma_L}{Leaf respiration coefficient, in year\eqn{^{-1}}.}
#'       \item{gamma_LN}{Leaf respiration coefficient per unit N.}
#'       \item{gamma_SW}{Sapwood respiration rate, in kg C m\eqn{^{-2}} year\eqn{^{-1}}.}
#'       \item{gamma_FR}{Fine root respiration rate, kg C kg C\eqn{^{-1}}
#'         year\eqn{^{-1}}.}
#'       \item{tc_crit}{Critical temperature triggerng offset of phenology, in Kelvin.}
#'       \item{tc_crit_on}{Critical temperature triggerng onset of phenology, in Kelvin.}
#'       \item{gdd_crit}{Critical value of GDD5 for turning ON growth season.}
#'       \item{betaON}{Critical soil moisture for phenology onset.}
#'       \item{betaOFF}{Critical soil moisture for phenology offset.}
#'       \item{seedlingsize}{Initial size of seedlings, in kg C per individual.}
#'       \item{LNbase}{Basal leaf N per unit area, in kg N m\eqn{^{-2}}.}
#'       \item{lAImax}{Maximum crown LAI (leaf area index) (not used, see LAI_light).}
#'       \item{Nfixrate0}{Reference N fixation rate (kg N kg C\eqn{^{-1}} root).}
#'       \item{NfixCost0}{Carbon cost of N fixation (kg C kg N\eqn{^{-1}}).}
#'       \item{phiCSA}{Ratio of sapwood area to leaf area.}
#'       \item{mortrate_d_c}{Canopy tree mortality rate (year\eqn{^{-1}}).}
#'       \item{mortrate_d_u}{Understory tree mortality rate (year\eqn{^{-1}}).}
#'       \item{maturalage}{Age at which trees can reproduce (years).}
#'       \item{v_seed}{Fraction of G_SF to G_F.}
#'       \item{fNSmax}{Multiplier for NSNmax as sum of potential bl and br.}
#'       \item{LMA}{Leaf mass per unit area (kg C m\eqn{^{-2}}).}
#'       \item{rho_wood}{Wood density (kg C m\eqn{^{-3}}).}
#'       \item{alphaBM}{Coefficient for allometry (biomass = alphaBM * DBH ** thetaBM).}
#'       \item{thetaBM}{Coefficient for allometry (biomass = alphaBM * DBH ** thetaBM).}
#'       \item{kphio}{Quantum yield efficiency \eqn{\varphi_0},
#'        in mol mol\eqn{^{-1}}.}
#'       \item{phiRL}{Ratio of fine root to leaf area.}
#'       \item{LAI_light}{Maximum LAI limited by light.}
#'     }}
#'   \item{init_cohort}{A data.frame of initial cohort specifications, including
#'   the following data:
#'     \describe{
#'       \item{init_cohort_species}{Index of a species described in param_species.}
#'       \item{init_cohort_nindivs}{Initial individual density, in individuals per
#'         m\eqn{^{2}}.}
#'       \item{init_cohort_bl}{Initial biomass of leaf, in kg C per individual.}
#'       \item{init_cohort_br}{Initial biomass of fine root, in kg C per individual.}
#'       \item{init_cohort_bsw}{Initial biomass of sapwood, in kg C per individual.}
#'       \item{init_cohort_bHW}{Initial biomass of heartwood, in kg C per individual.}
#'       \item{init_cohort_seedC}{Initial biomass of seed, in kg C per individual.}
#'       \item{init_cohort_nsc}{Initial non-structural biomass, in kg C per individual.}
#'     }}
#'   \item{init_soil}{A data.frame of initial soil pools, including
#'   the following data:
#'     \describe{
#'       \item{init_fast_soil_C}{Initial fast soil carbon, in kg C m\eqn{^{-2}}.}
#'       \item{init_slow_soil_C}{Initial slow soil carbon, in kg C m\eqn{^{-2}}.}
#'       \item{init_Nmineral}{Mineral nitrogen pool, in kg N m\eqn{^{-2}}.}
#'       \item{N_input}{Annual nitrogen input to soil N pool, in kg N m\eqn{^{-2}}
#'         year\eqn{^{-1}}.}
#'     }}
#' }
"biomee_gs_leuning_drivers"

#' rsofun BiomeE driver data (P-model photosynthesis model)
#' 
#' Small dataset representing the driver to run the BiomeE-model at the CH-LAE site
#' using the P-model photosynthesis specification (and daily time step).
#' It can also be used together with leaf trait data from CH-LAE (\code{\link{biomee_validation}}) 
#' to optimize model parameters.
#'
#' @format See \code{\link{biomee_gs_leuning_drivers}}
#'
#' @inherit biomee_gs_leuning_drivers source
"biomee_p_model_drivers"

#' rsofun BiomeE targets validation data
#'
#' Small example dataset of target observations (leaf trait data) at the CH-LAE site
#' to optimize model parameters with the function \code{\link{calib_sofun}}
#'
#' @format A tibble of validation data:
#' \describe{
#'   \item{sitename}{site name}
#'   \item{data}{validation data}
#' }
#'
#' @source Lukas Hörtnagl,  Werner Eugster,  Nina Buchmann,  Eugenie Paul-Limoges,  Sophia Etzold,  Matthias Haeni,  Peter Pluess,  Thomas Baur  (2004-2014)
#' FLUXNET2015 CH-Lae Laegern,
#' Dataset. https://doi.org/10.18140/FLX/1440134
"biomee_validation"

#' rsofun BiomeE (gs_leuning) output data
#'
#' Example output dataset from a BiomeE-model run using divers \code{\link{biomee_gs_leuning_drivers}}
#' See \code{\link{run_biomee_f_bysite}} for a detailed description of the outputs.
"biomee_gs_leuning_output"

#' rsofun BiomeE (P-model) output data
#'
#' Example output dataset from a BiomeE-model run using divers \code{\link{biomee_p_model_drivers}}
#' See \code{\link{run_biomee_f_bysite}} for a detailed description of the outputs.
"biomee_p_model_output"
