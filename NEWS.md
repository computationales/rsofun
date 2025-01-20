# rsofun (development version v5.1.0)

* new BiomeE behavior to recycle last year of forcing if requested simulation time span (`nyeartrend`) is longer than 
available forcing data
* Breaking changes in biomee drivers
  * `init_cohort$init_n_cohorts` column has been phased out and must not be present in  
drivers to protect against data corruption.
  * `update_annualLAImax` flag was removed and the behavior is now slightly altered since `LAImax` and `underLAImax` are set and kept constant at the start of the simulation:
    * `LAImax` takes value `max(LAI_light, 0.5)`
    * `underLAImax` takes value `min(LAImax, 1.2)`

# rsofun v5.0.0

* new BiomeE forcing data matching that of P-model
  * `prec` is now called `rain`
  * `rh` is now provided as `vpd`
  * See `biomee_gs_leuning_drivers` for an example
* fix Fortran modules leading to segmentation faults using BiomeE model
* improved documentation

# rsofun v4.4.1

* bugfix Fortran modules and derived types
* compilation flags for different platforms

# rsofun v4.4

* bugfix on input data type P-model
* LM3-PPA to BiomeE renaming
* cost function rewrite
* update output format (consistency)
* add water balance variables to P-model output
* new p-model calibratable parameters
* rewrite of temperature and soil moisture stress functions
* update simulation parameters to take netrad, ppfd and ccov as input forcing
* documentation parameter sensitivity analysis
* documentation on the data format
* document calibration diagnostic and model uncertainty calculation
* more transparent licensing and copyright statements in COPYING file

# rsofun v4.3

* Consistent variable names p-model / biomee
* optimization stability BayesianTools
* documentation
* Bugfixes

# rsofun v4.2

* Canopy transfer solved for BiomeE with P-model
* Bugfixes

# rsofun v4.1

* Catching aborting in FORTRAN part to avoid R session crashes

# rsofun v4.0

* Public release after refactored code
