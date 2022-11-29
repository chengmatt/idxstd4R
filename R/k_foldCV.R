#' Conducts k-fold cross validation for cpue_standardization models
#'
#' @param variables A vector of variables we want to consider
#' @param control_variables A formula in character form that species the response, as a function
#' of whatever control variables you want to include in all your models
#' @param data: A dataframe of the data you want to use
#' @param dir.models: The directory to your output folder
#' @param model_filename: The filename you want to give to your output object
#' @param error_dist: The error distribution you want to use for your standardization model. The default is a Tweedie Distribution.
#' @param stand_type: CPUE standardization model type. The default is to use mgcv to employ GAMs. Another option is to use GLMs
#' @param seed: Set seed for k fold random splits
#' @param folds: Number of folds you want to conduct
#'
#' @return A list object of 5-fold cross validated standardized CPUE models. This object is called all_mods_cv_random
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
#' k_foldCV(variables = possible_variables, control_variables = control_variables,
#' data = data, dir.models = dir.mod.output, model_filename = "cpue_standardization_5foldCV.RData",
#' error_dist = tw(link = "log"), folds = 5)

k_foldCV <- function(data, variables, control_variables, dir.models,
                        model_filename, error_dist = tw(link = "log"), stand_type = "gam",
                        seed = 666, folds = 5) {

  require(tidyverse)
  require(here)
  require(mgcv)
  require(beepr)
  require(doMC)
  require(rsample)

  # Set up for success - parralellization!
  ncores <- detectCores()
  registerDoMC(cores = ((ncores/2))-1)

  # Renaming data object
  log_obs_em <- data

  # Create random splits here
  set.seed(seed)
  k_folds <- rsample::vfold_cv(log_obs_em, v = folds)

  # Next, put these samples in a list, so we can re-use our previous workflow...
  folds <- list()

  # Loop through to put k folds into the folds list
  for(i in 1:length(k_folds)) {
    folds[[i]] <- k_folds$splits[[i]]$in_id
  }

  # Vector to store variable combinations
  possible_variables <- variables

  # Get all possible subsets
  variables <- unlist(lapply(1:length(possible_variables), combn,
                             x = possible_variables, simplify = FALSE),
                      recursive = FALSE)

  start_time <- Sys.time()
  # Next, let's start the CV process!

  # Set up list of CV length variables
  all_mods_cv_random <- vector("list", length = length(variables))

  for (i in 1:length(variables)) {

    # Setting up loop for model formula
    formula <- as.formula(paste(control_variables, paste(variables[[i]], collapse = "+")))

    cross_vald_mods_random <- foreach(j = 1:length(folds)) %dopar%{

      tryCatch({ # tryCatch keeps the loop going despite errors w/ deviance

        if(stand_type == "gam") {
          model_all_cv <- bam(formula, data = log_obs_em[folds[[j]], ], family = error_dist)
        } # gam standardization

        if(stand_type == "glm") {
          model_all_cv <- glm(formula, data = log_obs_em[folds[[j]], ], family = error_dist)
        } # glm standardization

        # Make predictions
        predictions <- predict(model_all_cv, newdata = log_obs_em[-folds[[j]], ], type = "response")

        output <- list(model_all_cv, predictions)

      } , error = function(error) {cat("ERROR :",conditionMessage(error), "\n")}) # end try catch statement
    } # end foreach loop

    # Store for each output into our list of lists
    all_mods_cv_random[[i]] <- cross_vald_mods_random

    # Save model after it's done!
    save(all_mods_cv_random, file=file.path(dir.models, model_filename))
    print(paste0("done with model ",i,"/", paste(length(variables))))

  } # End for loop

} # end function
