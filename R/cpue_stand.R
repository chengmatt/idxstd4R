
#' CPUE standardization  for all subsets of possible variables
#'
#' @param variables A vector of variables we want to consider
#' @param control_variables A formula in character form that specifies the response, as a function
#' of whatever control variables you want to include in all your models
#' @param data: A dataframe of the data you want to use
#' @param dir.models: The directory to your output folder
#' @param model_filename: The filename you want to give to your output object
#' @param error_dist: The error distribution you want to use for your standardization model. The default is to use a Tweedie distribution
#' @param stand_type: CPUE standardization model type. The default is to use mgcv to employ GAMs. Another option is to use GLMs
#'
#'
#'
#' @return A list object of models called model_subsets gets outputted into the environment and into the output directory
#' @export
#'
#' @examples
#' # Set working directory
#' dir.mod.output <- file.path(here(), "output")
#'
#' # Create an object with possible variables we want to consider
#' possible_variables <- c("s(doy)", "te(AdjLon,lat2, bs = c('tp', 'tp'))",
#' "s(bottom_depth)",  "te(AdjLon, lat2, by = gear_descrip, bs = c('tp', 'tp'))", "area")
#'
#' # Control variables we want to keep in every model
#' control_variables <- "weight ~ -1 + year + vessel_length_factor + offset(log(total_hooks_pots)) + type + gear_descrip +"
#'
#' # Run function here!
#' cpue_stand(variables = possible_variables, control_variables = control_variables,
#' data = data, dir.models = dir.mod.output, model_filename = "cpue_standardization_allgear.RData",
#' error_dist = tw(link = "log"))


cpue_stand <- function(variables, control_variables, data, dir.models,
                       model_filename, error_dist = tw(link = "log"), stand_type = "gam") {
  require(tidyverse)
  require(here)
  require(mgcv)
  require(beepr)
  require(doMC)

  # Set up parralellization process
  ncores <- detectCores()
  registerDoMC(cores = ((ncores/2))-1)

  # Set up model ------------------------------------------------------------

  # Vector to store variable combinations
  possible_variables <- variables

  # Get all possible subsets of variables
  variables <- unlist(lapply(1:length(possible_variables), combn,
                             x = possible_variables, simplify = FALSE),
                      recursive = FALSE)

  # Fit model here ----------------------------------------------------------
  log_obs_em <- data # Rename to fit through into models

  models_subsets <- foreach(i = 1:length(variables)) %dopar% {

    tryCatch({ # tryCatch keeps the loop going despite erros w/ deviance
      # NaN

      # Setting up loop for model formula
      # Model formula set up, and going through all combinations here
      formula <- as.formula(paste(control_variables, paste(variables[[i]], collapse = "+")))

      if(stand_type == "gam") {
        # Pipe this formula into our model formula for BAM
        model_all_dat <- bam(formula, data = log_obs_em, family = error_dist)
      } # if we are using gams for standardization

      if(stand_type == "glm") {
        # Pipe this formula into our model formula for BAM
        model_all_dat <- glm(formula, data = log_obs_em, family = error_dist)
      } # if we are using glms for standardization

      model_all_dat # Need to give the parrallelized code something to output - we want
      # our output to be a list of model objects

    }, error = function(error) {cat("ERROR :",conditionMessage(error), "\n")}) # end try catch statement to overide convergence errors
  } # End foreach loop

  # Save model objects
  save(models_subsets, file=file.path(dir.models, model_filename))

} # end fxn
