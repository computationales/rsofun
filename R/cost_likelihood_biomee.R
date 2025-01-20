#' Log-likelihood cost function for BiomeE with different targets
#' 
#' Cost function for parameter calibration, which
#' computes the log-likelihood for the biomee model fitting several target 
#' variables for a given set of parameters.
#' 
#' @param par A vector containing parameter values for \code{'phiRL',
#' 'LAI_light', 'tf_base', 'par_mort'} in that order, and for the error terms
#' corresponding to the target variables, e.g. \code{'err_GPP'} if GPP is a target. 
#' Make sure that
#' the order of the error terms in \code{par} coincides with the order provided in
#' the \code{targets} argument.
#' @param obs A nested data frame of observations, following the structure of \code{biomee_validation},
#' for example.
#' @param drivers A nested data frame of driver data, for example \code{biomee_gs_leuning_drivers}.
#' @param targets A character vector indicating the target variables for which the
#' optimization will be done. This should be a subset of \code{c("GPP", "LAI",
#' "Density", "Biomass")}.
#' 
#' @return The log-likelihood of the simulated 
#' targets by the biomee model versus the observed targets. 
#' 
#' @details The cost function performs a BiomeE model run for the value of
#' \code{par} given as argument. The likelihood is calculated assuming that the 
#' predicted targets are independent, normally distributed and centered on the observations. 
#' The optimization 
#' should be run using \code{BayesianTools}, so the likelihood is maximized.
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' # Compute the likelihood for a set of
#' # BiomeE model parameter values
#' # and the example data
#' cost_likelihood_biomee(
#'  par = c(3.5, 3.5, 1, 1,    # model params
#'          0.5),              # err_GPP
#'  obs = biomee_validation,
#'  drivers = biomee_gs_leuning_drivers,
#'  targets = c("GPP")
#' )
#' }

cost_likelihood_biomee <- function(
    par,
    obs,
    drivers,
    targets
){
  
  # predefine variables for CRAN check compliance
  GPP <- LAI <- Density12 <- plantC <- NULL
  
  # Add changed model parameters to drivers, overwriting where necessary.
  drivers$params_species[[1]]$phiRL[]  <- par[1]
  drivers$params_species[[1]]$LAI_light[]  <- par[2]
  drivers$params_tile[[1]]$tf_base <- par[3]
  drivers$params_tile[[1]]$par_mort <- par[4]
  
  # run model
  df <- runread_biomee_f(
    drivers,
    makecheck = TRUE,
    parallel = FALSE
  )
  
  # Aggregate variables from the model df taking the last 50 yrs (up to)
  df <- df$data[[1]]$output_annual_tile |>
    utils::tail(50) |>
    dplyr::summarise(
      GPP = mean(GPP),
      LAI = stats::quantile(LAI, probs = 0.95, na.rm=TRUE),
      Density = mean(Density12),
      Biomass = mean(plantC)
    )
  
  # reshuffle observed data
  col_names <- obs$data[[1]]$variables
  obs <- data.frame(t(obs$data[[1]]$targets_obs))
  colnames(obs) <- col_names
  
  # calculate the log likelihood, loop over targets
  ll <- lapply(seq(length(targets)), function(i){
    target <- targets[i]
    BayesianTools::likelihoodIidNormal(
      predicted = df[[target]],
      observed = obs[[target]],
      sd = par[4+i]
    )
  }) |>
    unlist() |>
    sum()     # sum log-likelihoods
  
  # trap boundary conditions
  if(is.nan(ll) || is.na(ll) | ll == 0){
    ll <- -Inf
  }
  
  return(ll)
}
