% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_insample_metrics.R
\name{extract_insample_metrics}
\alias{extract_insample_metrics}
\title{Extract insample metrics (BIC, AIC, R2, Year Effects) from cpue_stand function.}
\usage{
extract_insample_metrics(
  variables,
  model_object,
  dir.models,
  metrics_ic_path,
  year_coeffs_path
)
}
\arguments{
\item{variables:}{A vector of variables we want to consider}

\item{model_object:}{An RData file or a list object that contains your models}

\item{data:}{A dataframe of the data you want to use}

\item{dir.models:}{The directory to your output folder}

\item{metrics_ic_path:}{The filename and path you want your in sample metrics (BIC}

\item{year_coeffs_path:}{The filename and path you want your indices (year coefficients) and
corresponding standard errors to be called.}
}
\value{
CSVs of year effects and model metrics into directed folder path
}
\description{
Extract insample metrics (BIC, AIC, R2, Year Effects) from cpue_stand function.
}
\examples{
# Load in model object output from the cpue_stand function
load(file.path(dir.mod.output, "models.RData"))

# Extract in-sample metrics here from our models that were fit to all the data.
extract_insample_metrics(variables = possible_variables,
                       model_object = model_filename,
                       dir.models = dir.mod.output,
                       metrics_ic_path = "insample_metrics.csv",
                       year_coeffs_path = "year_coeffs,csv",
                       data = data)


}
