% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frequency.R
\name{mz_frequency_interpolation}
\alias{mz_frequency_interpolation}
\title{convert mz to frequency using linear fit}
\usage{
mz_frequency_interpolation(
  mz_query,
  mz_values = NULL,
  frequency_values = NULL,
  model = NULL
)
}
\arguments{
\item{mz_query}{the M/Z value to fit}

\item{mz_values}{two M/Z values}

\item{frequency_values}{two frequency values}

\item{model}{a model to use instead of actual values}
}
\value{
data.frame with predicted_value, intercept, and slope
}
\description{
Given a query, and either two values of M/Z and two values of frequency or
a previously generated model, return a data.frame with the predicted value,
and the slope and the intercept so the model can be re-used later for
other points when needed.
}
