#' Extract insample metrics (BIC, AIC, R2, Year Effects) from cpue_stand function.
#'
#' @param variables: A vector of variables we want to consider
#' @param model_object: An RData file or a list object that contains your models
#' @param data: A dataframe of the data you want to use
#' @param dir.models: The directory to your output folder
#' @param metrics_ic_path: The filename and path you want your in sample metrics (BIC
# or AIC, Adjusted R2, Deviance Explained) to be called
#' @param year_coeffs_path: The filename and path you want your indices (year coefficients) and
#' corresponding standard errors to be called.
#'
#' @return CSVs of year effects and model metrics into directed folder path
#' @export
#'
#' @examples
#' # Load in model object output from the cpue_stand function
#' load(file.path(dir.mod.output, "models.RData"))
#'
#' # Extract in-sample metrics here from our models that were fit to all the data.
#' extract_insample_metrics(variables = possible_variables,
#'                        model_object = model_filename,
#'                        dir.models = dir.mod.output,
#'                        metrics_ic_path = "insample_metrics.csv",
#'                        year_coeffs_path = "year_coeffs,csv",
#'                        data = data)
#'
#'
extract_insample_metrics <- function(variables, model_object,
                                   dir.models, metrics_ic_path, year_coeffs_path)  {

  require(tidyverse)

  # Rename model_object to model_subsets
  model_subsets <- model_object

  # Vector to store variable combinations
  possible_variables <- variables

  # Get all possible subsets
  variables <- unlist(lapply(1:length(possible_variables), combn,
                             x = possible_variables, simplify = FALSE),
                      recursive = FALSE)

  # Set up dataframe to store metrics ---------------------------------------

  # Create year coefficient object to store se.fit and coeffs
  year_coeffs <- data.frame(matrix(data = NA, nrow = length(unique(data$year)),
                                   ncol = length(variables)),
                            year = sort(unique(data$year)))

  se_fit <- data.frame(matrix(data = NA, nrow = length(unique(data$year)),
                              ncol = length(variables)),
                       year = sort(unique(data$year)))

  metrics_ic <- data.frame(matrix(data = NA, nrow = 4,
                                  ncol = length(variables)),
                           type = c("AIC", "BIC", "dev", "adj-r2"))


  # Loop through to get metrics here ----------------------------------------

  for (i in 1:length(variables)) {
    tryCatch({
      # TESTING
      # i <- 203
      if (!is.null(models_subsets[[i]])){ # For models that are NULL, did not converge, etc.

        # Rename model runs as list object
        names(models_subsets)[i] <- paste("model", i) # Rename model run

        # Extract year coefficients and standard error fit!
        year_coeffs[,i] <- coef(models_subsets[[i]])[1:length(unique(data$year))]
        names(year_coeffs)[i] <- paste("model", i) # rename col names

        # Standard error fit
        se_fit[,i] <- summary(models_subsets[[i]])$se[1:length(unique(data$year))]
        names(se_fit)[i] <- paste("model", i) # rename col names

        # AIC and BIC calculation!
        metrics_ic[1,i] <- AIC(models_subsets[[i]]) # AIC
        metrics_ic[2,i] <- BIC(models_subsets[[i]]) # BIC
        names(metrics_ic)[i] <- paste("model", i) # rename col names

        # Extract deviance explained and adj-r2
        metrics_ic[3,i] <- summary(models_subsets[[i]])$dev.expl
        metrics_ic[4,i] <- summary(models_subsets[[i]])$r.sq

        print(paste("done with model", i, "/", length(variables))) # Printing progress

      }
    }, error = function(error) {cat("ERROR :",conditionMessage(error), "\n")}) # End tryCatch fxn

    if (i == length(variables)) { # When it's done running through these, then pivot, and clean the data

      # All we're doing here is pivoting from wide to long for all our metrics
      metrics_ic_alldat <- metrics_ic %>%
        pivot_longer(!"type", values_to = "metric_values", names_to = "values") %>%
        drop_na()

      # Binding year coefficients with se.fit
      year_coeffs_alldat <- year_coeffs %>%
        pivot_longer(!year, values_to = "year_coeffs", names_to = "model") %>%
        drop_na()

      se_fit_alldat <- se_fit %>%
        pivot_longer(!year, values_to = "se_fit", names_to = "model") %>%
        dplyr::select(-model) %>%
        drop_na()

      year_se_alldat <- cbind(year_coeffs_alldat, se_fit_alldat[-1])

      write.csv(metrics_ic_alldat, file.path(dir.models, metrics_ic_path))
      write.csv(year_se_alldat, file.path(dir.models, year_coeffs_path))
    }
  } # End loop

} # end function
