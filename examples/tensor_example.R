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
data <- data.frame(doy = rpois(5000, 150), AdjLon = rgamma(5000, 5, 3),
                   lat2 = rgamma(5000, 10, 5), bottom_depth = rgamma(5000, 5000, 10),
                   gear_descrip = sample(c("LONGLINER", "POT"), replace = TRUE, size = 5000),
                   area = sample(c("a", "b", "c", replace = TRUE, size = 5000)),
                   type = sample(c("observer", "logbook"), replace = TRUE, size = 5000),
                   weight = exp(rnorm(5000, 3, 1)), total_hooks_pots = rnorm(5000, 100, 10),
                   year = sample(factor(seq(1995:2020)), replace = TRUE, size = 5000))

# Variables that are considered
possible_variables <- c("s(doy)", "te(AdjLon,lat2, bs = c('tp', 'tp'))", "s(bottom_depth)",
                        "te(AdjLon, lat2, by = gear_descrip, bs = c('tp', 'tp'))", "area")

# Here, we are considering a day of year effect, a tensor product smooth (average spatial field),
# bottom depth, tensor product smooth with gear specific catchabilities, and fishery management areas

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

five_foldCV(data = data, variables = possible_variables, control_variables = control_variables,
            dir.models = dir.out, model_filename = "test_5fold.RData", error_dist = tw(link = "log"),
            stand_type = "gam", seed = 123)

# Load in 5 fold cross validated model objects - these are called all_mods_cv_random
load(here(dir.out, "test_5fold.RData"))

# Extract 5 fold cv metrics from the above --------------------------------

extract_5fold_metrics(variables = possible_variables, model_object = all_mods_cv_random,
                        data = data, metric_filename = "5foldcv.csv", dir.models = dir.out)
