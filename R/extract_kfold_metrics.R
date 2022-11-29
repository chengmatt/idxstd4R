#' Extracts metrics from 5-fold cross validated models
#'
#' @param variables A vector of variables we want to consider
#' @param model_object The model object we want to extract metrics from
#' @param data: A dataframe of the data you want to use
#' @param dir.models The directory to your output folder
#' @param metric_filename The filename you want to give to your output object (metrics of CV; i.e., RMSE, MAE, R2)
#' @param folds Number of fold we want to extract metrics out of
#' @param seed Seed used for k_fold_cv function for random splits
#'
#' @return CSVs of out-of-sample metrics derived from the 5-fold cross validation process
#' @export
#'
#' @examples
#' # Load in k fold CV model object
#' load(file.path(dir.mod.output, "cpue_stand_5foldcv.RData"))
#'
#' # Extract k fold CV metrics here
#' extract_metrics_5foldcv(variables = possible_variables,
#'                         model_object = cpue_stand_5foldcv,
#'                         data = data,
#'                         dir.models = dir.mod.output,
#'                         metric_filename =  "allgear_cv_metrics.csv",
#'                         folds = 5, seed = 666)
#'
extract_kfold_metrics <- function(variables, model_object, data,
                                    dir.models, metric_filename, folds, seed) {

  require(rsample)
  require(yardstick)
  require(tidyverse)

  # Vector to store variable combinations
  possible_variables <- variables

  # Get all possible subsets
  variables <- unlist(lapply(1:length(possible_variables), combn,
                             x = possible_variables, simplify = FALSE),
                      recursive = FALSE)

  # Rename model object here
  all_mods_cv_random <- model_object

  # Get folds
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

  # Create dataframe to store CV metrics in
  cv_results_random <- data.frame(matrix(data = NA, nrow = 15, length(variables)),
                                  metrics = c(rep("MAE", 5), rep("R2", 5), rep("RMSE", 5)))

  # Loop through this here
  for (i in 1:length(variables)) {

    tryCatch({
      # TESTING
      # i <- 1
      if (!is.null(all_mods_cv_random[[i]])){
        # Rename model runs as list object
        names(all_mods_cv_random)[i] <- paste("model", i) # Rename model run

        for (j in 1:length(folds)){
          # TESTING
          # j <- 1
          # Calculate CV metrics
          # Here, we are just subsetting the numbers we need ahead of time, before we put
          # the into the MAE, RMSE, and R2 functions.

          cv_res <- as.numeric(all_mods_cv_random[[i]][[j]][[2]])
          k_fold_obs <- as.numeric(log_obs_em[-folds[[j]],]$weight)

          # Put these in a dataframe
          cv_df <- data.frame(truth = k_fold_obs, estimate = cv_res)

          # MAE
          cv_results_random[j,i] <- yardstick::mae(data = cv_df, truth = truth , estimate, na.rm = TRUE)[[3]]

          # R2 error
          cv_results_random[j+5,i] <- yardstick::rsq(cv_df, truth, estimate, na.rm = TRUE)[[3]]

          # RMSE error
          cv_results_random[j+10,i] <- yardstick::rmse(cv_df, truth, estimate, na.rm = TRUE)[[3]]

          # Extract BIC, AIC from these as well
          cv_results_random[j+15, i] <- BIC(all_mods_cv_random[[i]][[j]][[1]])
          cv_results_random[j+20, i] <- AIC(all_mods_cv_random[[i]][[j]][[1]])
        }

        print(paste("done with model", i, "/", length(variables))) # Printing progress

        # Rename col names!
        names(cv_results_random)[i] <- paste("model", i)

      } # if statement here
    }, error = function(error) {cat("ERROR :",conditionMessage(error), "\n")}) # End tryCatch fxn

    if (i == length(variables)) { # If we are done with extracting all the CV metrics, let's
      # tidy up the data frame, and export this

      cv_results_random_long_random <- cv_results_random %>%
        pivot_longer(!"metrics", values_to = "metric_values", names_to = "models") %>%
        drop_na()

      # Write this out as a csv
      write.csv(cv_results_random_long_random, file.path(dir.models, metric_filename))
    }
  } # End loop


} # end function
