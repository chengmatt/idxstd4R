# Purpose: To present examples of how to implement tensor product smooths, offsets, and general introduction
# to the idxstd4R package
# NOTE: This is a rough example of how to use idxstd4R.
# Creator: Matthew LH. Cheng
# Date: 11/27/22


# Set up ------------------------------------------------------------------

library(here)
library(idxstd4R)

# Create directory to output
dir.out <- here()

# Create a random dataframe as an example
data <- data.frame(doy = rpois(1500, 150), AdjLon = rgamma(1500, 5, 3),
                   lat2 = rgamma(1500, 10, 5), bottom_depth = rgamma(1500, 1500, 10),
                   gear_descrip = sample(factor(c("LONGLINER", "POT")), replace = TRUE, size = 1500),
                   area = sample(c("a", "b", "c", replace = TRUE, size = 1500)),
                   type = sample(factor(c("observer", "logbook")), replace = TRUE, size = 1500),
                   weight = exp(rnorm(1500, 3, 1)), total_hooks_pots = rnorm(1500, 100, 10),
                   year = sample(factor(seq(1995:2020)), replace = TRUE, size = 1500))

# Variables that are considered
possible_variables <- c("s(doy, bs = 'cc')", "te(AdjLon,lat2, bs = c('tp', 'tp'))", "s(bottom_depth)",
                        "te(AdjLon, lat2, by = gear_descrip, bs = c('tp', 'tp'))", "area")

# s(doy, bs = 'cc') = smooth effect of day of year with a cyclic spline
# te(AdjLon,lat2, bs = c('tp', 'tp')) = tensor product smooth of lat and lon (avg spatial field)
# s(bottom_depth) = smooth effect of bottom depth as a thin plate spline
# te(AdjLon, lat2, by = gear_descrip, bs = c('tp', 'tp')) = tensor product smooth of lat and lon,
# estiamting independent spatial fields for each independent gear type.

# Variables that we know impact catchability
control_variables <- "weight ~ year + offset(log(total_hooks_pots)) + type + gear_descrip +"

# Above are our control variables that are meant to be included in all model variants


# CPUE standardization process --------------------------------------------

# Now, call our CPUE standardization function
cpue_stand(variables = possible_variables, control_variables = control_variables,
           data = data, dir.models = dir.out, model_filename = "test.RData",
           error_dist = tw(link = "log"), stand_type = "gam")

# Now, load in our object - called model_subsets
load(here(dir.out, "test.RData"))


# Extract insample metrics from the model above ------------------------------------------------

extract_insample_metrics(variables = possible_variables, model_object = models_subsets,
                         dir.models = dir.out, metrics_ic_path = "insample.csv",
                         year_coeffs_path = "year.csv")


# Five fold cross validation ----------------------------------------------

# Note that the seed has to be the same for the k_foldCV function and the extract metrics function
k_foldCV(data = data, variables = possible_variables, control_variables = control_variables,
            dir.models = dir.out, model_filename = "test_5fold.RData", error_dist = tw(link = "log"),
            stand_type = "gam", seed = 123, folds = 5)

# Load in 5 fold cross validated model objects - these are called all_mods_cv_random
load(here(dir.out, "test_5fold.RData"))


# Extract 5 fold cv metrics from the above --------------------------------

# Note that the seed has to be the same for the k_foldCV function and the extract metrics function
extract_kfold_metrics(variables = possible_variables, model_object = all_mods_cv_random,
                        data = data, metric_filename = "5foldcv.csv", dir.models = dir.out,
                        seed = 123, folds = 5)
