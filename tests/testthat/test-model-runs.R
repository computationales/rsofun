context("test models and their parameters")
set.seed(10)

test_that("p-model run check GPP", {
  skip_on_cran()

  # load parameters (valid ones)
  params_modl <- list(
    kphio              = 0.04998, # setup ORG in Stocker et al. 2020 GMD
    kphio_par_a        = 0.01,  # set to zero to disable temperature-dependence of kphio, setup ORG in Stocker et al. 2020 GMD
    kphio_par_b        = 1.0,
    soilm_thetastar    = 0.6 * 240,  # to recover old setup with soil moisture stress
    soilm_betao        = 0.01,
    beta_unitcostratio = 146.0,
    rd_to_vcmax        = 0.014, # value from Atkin et al. 2015 for C3 herbaceous
    tau_acclim         = 30.0,
    kc_jmax            = 0.41
  )

  # read in demo data
  df_drivers <- p_model_drivers

  # run the SOFUN Fortran P-model
  mod <- run_pmodel_f_bysite(
    df_drivers$sitename[1],
    df_drivers$params_siml[[1]],
    df_drivers$site_info[[1]],
    df_drivers$forcing[[1]],
    params_modl = params_modl,
    makecheck = FALSE
  )

  # test if the returned values
  # are in a list (don't error / warning)
  expect_type(mod, "list")

  # test runread_pmodel_f
  df_output <- runread_pmodel_f(
    df_drivers,
    par = params_modl,
    makecheck = FALSE,
    parallel = FALSE
  )

  # test for correctly returned values
  expect_type(df_output, "list")

  # test runread_pmodel_f
  df_output_p <- runread_pmodel_f(
    df_drivers,
    par = params_modl,
    makecheck = TRUE,
    parallel = TRUE
  )

  # test for correctly returned values
  expect_type(df_output_p, "list")
})

test_that("p-model run check Vcmax25", {
  skip_on_cran()

  # load parameters (valid ones)
  params_modl <- list(
    kphio              = 0.04998, # setup ORG in Stocker et al. 2020 GMD
    kphio_par_a        = 0.01,  # set to zero to disable temperature-dependence of kphio, setup ORG in Stocker et al. 2020 GMD
    kphio_par_b        = 1.0,
    soilm_thetastar    = 0.6 * 240,  # to recover old setup with soil moisture stress
    soilm_betao        = 0.01,
    beta_unitcostratio = 146.0,
    rd_to_vcmax        = 0.014, # value from Atkin et al. 2015 for C3 herbaceous
    tau_acclim         = 30.0,
    kc_jmax            = 0.41
  )

  # read in demo data
  df_drivers <- p_model_drivers_vcmax25

  # run the SOFUN Fortran P-model
  mod <- run_pmodel_f_bysite(
    df_drivers$sitename[1],
    df_drivers$params_siml[[1]],
    df_drivers$site_info[[1]],
    df_drivers$forcing[[1]],
    params_modl = params_modl,
    makecheck = FALSE
  )

  # test if the returned values
  # are in a list (don't error / warning)
  expect_type(mod, "list")

  # test runread_pmodel_f
  df_output <- runread_pmodel_f(
    df_drivers,
    par = params_modl,
    makecheck = FALSE,
    parallel = FALSE
  )

  # test for correctly returned values
  expect_type(df_output, "list")

  # test runread_pmodel_f
  df_output_p <- runread_pmodel_f(
    df_drivers,
    par = params_modl,
    makecheck = TRUE,
    parallel = TRUE,
    ncores = 1
  )

  # test for correctly returned values
  expect_type(df_output_p, "list")
})

test_that("biomeE output check (gs leuning)", {
  skip_on_cran()

  out <- runread_biomee_f(
    biomee_gs_leuning_drivers,
    makecheck = TRUE,
    parallel = FALSE)

  expect_true(all.equal(colMeans(out$data[[1]]$output_daily_tile), colMeans(biomee_gs_leuning_output$data[[1]]$output_daily_tile), tolerance = 1e-4))
  expect_true(all.equal(colMeans(out$data[[1]]$output_annual_tile), colMeans(biomee_gs_leuning_output$data[[1]]$output_annual_tile), tolerance = 1e-4))
  expect_true(all.equal(colMeans(out$data[[1]]$output_annual_cohorts), colMeans(biomee_gs_leuning_output$data[[1]]$output_annual_cohorts), tolerance = 1e-4))
  expect_true(all.equal(colMeans(out$data[[1]]$output_annual_land_use), colMeans(biomee_gs_leuning_output$data[[1]]$output_annual_land_use), tolerance = 1e-4))

  # If this test fails it means that the output of the model is out of sync with the data in the data directory.
  # It could either mean that:
  # - the model was accidentally altered and should be fixed to deliver the expected output
  # - the model, drivers, or parameters was changed and the output data needs to be re-generetaed using the scripts in
  #   raw-data directory.
})


test_that("biomee output check (p-model)", {
  skip_on_cran()

  out <- runread_biomee_f(
    biomee_p_model_drivers,
    makecheck = TRUE,
    parallel = FALSE)

  expect_true(all.equal(colMeans(out$data[[1]]$output_daily_tile), colMeans(biomee_p_model_output$data[[1]]$output_daily_tile), tolerance = 1e-4))
  expect_true(all.equal(colMeans(out$data[[1]]$output_annual_tile), colMeans(biomee_p_model_output$data[[1]]$output_annual_tile), tolerance = 1e-4))
  expect_true(all.equal(colMeans(out$data[[1]]$output_annual_cohorts), colMeans(biomee_p_model_output$data[[1]]$output_annual_cohorts), tolerance = 1e-4))
  expect_true(all.equal(colMeans(out$data[[1]]$output_annual_land_use), colMeans(biomee_p_model_output$data[[1]]$output_annual_land_use), tolerance = 1e-4))

  # If this test fails it means that the output of the model is out of sync with the data in the data directory.
  # It could either mean that:
  # - the model was accidentally altered and should be fixed to deliver the expected output
  # - the model, drivers, or parameters was changed and the output data needs to be re-generetaed using the scripts in
  #   raw-data directory.
})

test_that("biomee parallel run check (gs leuning)", {
  skip_on_cran()
  
  df_drivers <- biomee_gs_leuning_drivers
  df_drivers$params_siml[[1]]$spinup <- FALSE
  
  df_output <- runread_biomee_f(
    df_drivers,
    makecheck = FALSE,
    parallel = TRUE,
    ncores = 2
  )
  
  # test for correctly returned values
  expect_type(df_output, "list")
  
})
