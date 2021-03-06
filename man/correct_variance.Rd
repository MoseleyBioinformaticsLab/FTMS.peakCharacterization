% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peak_sd_intensity_corrections.R
\name{correct_variance}
\alias{correct_variance}
\title{correct peak variance}
\usage{
correct_variance(observed_variance, fraction)
}
\arguments{
\item{observed_variance}{the observed \strong{variance}}

\item{fraction}{what fraction was it observed in}
}
\value{
corrected variance
}
\description{
Given a variance observed from a truncated normal distribution, correct it
assuming that it should have had 100\% observaionts
}
\references{
https://en.wikipedia.org/wiki/Truncated_normal_distribution
}
\seealso{
correct_mean correct_peak
}
