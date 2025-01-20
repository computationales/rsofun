#' Run BiomeE
#'
#' Runs BiomeE model for multiple sites.
#'
#' @param drivers A nested data frame with one row for each site and columns
#' named according to the arguments of function \code{\link{run_biomee_f_bysite}}.
#' Namely \code{sitename, params_siml, site_info} and \code{forcing}.
#' @param makecheck A logical specifying whether checks are performed 
#'  to verify forcings and model parameters. \code{TRUE} by default.
#' @param parallel Deprecated. Use ncores instead.
#' @param ncores An integer specifying the number of cores used for parallel 
#' computing (sites processed in parallel). Default: 1 (no parallel execution).
#'
#' @return A data frame (tibble) with one row for each site, site information 
#' stored in the nested column \code{site_info} and model outputs stored in the 
#' nested  column \code{data}. See \code{\link{run_biomee_f_bysite}} for a detailed 
#' description of the outputs.
#' Example outputs are provided as \code{\link{p_model_output}} and
#' \code{\link{p_model_output_vcmax25}}.
#' @export
#' 
#' @examples 
#' \donttest{
#' # Example BiomeE model run
#' 
#' runread_biomee_f(
#'   drivers = biomee_gs_leuning_drivers
#' )
#' runread_biomee_f(
#'   drivers = biomee_p_model_drivers
#' )
#' }

runread_biomee_f <- function(
    drivers,
    makecheck = TRUE,
    parallel = FALSE, # Not used
    ncores = 1
){
  
  # predefine variables for CRAN check compliance
  sitename <- params_siml <- site_info <- forcing <-
  params_tile <- params_species <- init_cohort <- init_soil <-
  init_lu <- luc_forcing <- . <- NULL

  if (parallel != (ncores > 1)) {
    warning("Warning: parallel flag is deprecated. Please set ncores to 1 to disable parallel execution.")
    parallel <- (ncores > 1)
  }

  parameters <- c(
    "sitename",
    "params_siml",
    "site_info",
    "forcing",
    "params_tile",
    "params_species",
    "init_cohort",
    "init_soil",
    "init_lu",
    "luc_forcing"
  )

  # Add potentially missing columns using NULL
  `%nin%` <- Negate(`%in%`)
  empty_col <- vector("list", nrow(drivers))
  if ('init_lu' %nin% colnames(drivers))
    drivers <- dplyr::mutate(drivers, 'init_lu'=empty_col)
  if ('luc_forcing' %nin% colnames(drivers))
    drivers <- dplyr::mutate(drivers, 'luc_forcing'=empty_col)

  df <- drivers %>%
    dplyr::select(any_of(parameters))
  
  if (ncores > 1){
    
    cl <- multidplyr::new_cluster(ncores) %>% 
      multidplyr::cluster_assign(makecheck = FALSE) %>% 
      multidplyr::cluster_library(
        packages = c("dplyr", "purrr", "rsofun")
      )

    df <- df %>%
      multidplyr::partition(cl) %>%
      dplyr::mutate('data' = purrr::pmap(list(
        sitename,
        params_siml,
        site_info,
        forcing,
        params_tile,
        params_species,
        init_cohort,
        init_soil,
        init_lu,
        luc_forcing
      ), run_biomee_f_bysite)) %>%
      dplyr::collect()
    
  } else {
    df <- df %>%
      dplyr::mutate('data' = purrr::pmap(
        .,
        run_biomee_f_bysite,
        makecheck = makecheck
      ))

  }

  df_out <- df |>
    dplyr::select('sitename', 'data') |>
    tidyr::unnest_wider('data')
  
  return(df_out)
}
