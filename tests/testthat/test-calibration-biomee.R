context("test BiomeE calibration framework and its parameters")
set.seed(10)

test_that("test calibration routine biomee (likelihood cost + Bayesiantools)", {
  skip_on_cran()
  df_drivers <- rsofun::biomee_gs_leuning_drivers
  ddf_obs <- rsofun::biomee_validation
  df_drivers$params_siml[[1]]$spinupyears <- 0
  
  settings <- list(
    method              = "bayesiantools",
    metric              = rsofun::cost_likelihood_biomee,
    control = list(
      sampler = "DEzs",
      settings = list(
        burnin = 1,
        iterations = 4,
        nrChains = 1
      )
    ),
    par = list(
      phiRL = list(lower = 0.5, upper = 5, init = 3.5),
      LAI_light = list(lower = 2, upper = 5, init = 3.5),
      tf_base = list(lower = 0.1, upper = 1, init = 0.5),
      par_mort = list(lower = 1, upper = 2, init = 1.1),
      
      # uncertainties
      err_GPP = list(lower = 0, upper = 30, init = 15),
      err_LAI = list(lower = 0, upper = 5, init = 3),
      err_Density = list(lower = 0, upper = 400, init = 280),
      err_Biomass = list(lower = 0, upper = 50, init = 45)
    )
  )
  
  pars <- rsofun::calib_sofun(
    drivers = df_drivers,
    obs = ddf_obs,
    settings = settings,
    # arguments for cost function
    targets = c("GPP","LAI","Density","Biomass")
  )
  
  # test for correctly returned values
  expect_type(pars, "list")
})

test_that("test calibration routine biomee (rmse cost + GenSA)", {
  skip_on_cran()
  df_drivers <- rsofun::biomee_gs_leuning_drivers
  ddf_obs <- rsofun::biomee_validation
  df_drivers$params_siml[[1]]$spinupyears <- 0
  
  settings <- list(
    method              = "gensa",
    # targets             = c("GPP","LAI","Density","Biomass"),
    metric              = rsofun::cost_rmse_biomee,
    control = list(
      maxit = 2
    ),
    par = list(
      phiRL = list(lower = 0.5, upper = 5, init = 3.5),
      LAI_light = list(lower = 2, upper = 5, init = 3.5),
      tf_base = list(lower = 0.1, upper = 1, init = 0.5),
      par_mort = list(lower = 1, upper = 2, init = 1.1)
    )
  )
  
  pars <- rsofun::calib_sofun(
    drivers = df_drivers,
    obs = ddf_obs,
    settings = settings
  )
  
  # test for correctly returned values
  expect_type(pars, "list")
})