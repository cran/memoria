% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/memoria.R
\docType{data}
\name{palaeodataLagged}
\alias{palaeodataLagged}
\title{Lagged data generated by \code{\link{prepareLaggedData}}.}
\format{dataframe with 13 columns and 3988 rows.}
\usage{
data(palaeodataLagged)
}
\description{
A dataframe resulting from the application of \code{\link{prepareLaggedData}} to the dataset \code{\link{palaeodata}}. The dataframe columns are:
}
\details{
\itemize{
  \item \emph{Response_0} numeric, values of the response variable selected by the user in the argument \code{response} of the function \code{\link{prepareLaggedData}}. This column is used as response variable by the function \code{\link{computeMemory}}. In this case, Response represent pollen counts of Pinus.
  \item \emph{Response_0.2-1} numeric, time delayed values of the response for different lags (in ky). Considered together these columns represent the endogenous ecological memory.
  \item \emph{climate.temperatureAverage_0} numeric, values of the variable temperatureAverage for the lag 0 (no lag). This column represents the concurrent effect of the temperature over the response.
  #'   \item \emph{climate.rainfallAverage_0} numeric, values of the variable rainfallAverage for the lag 0 (no lag). This column represents the concurrent effect of rainfall over the response.
  \item \emph{climate.temperatureAverage_0.2-1} numeric, time delayed values of temperatureAverage for lags 0.2 to 1 (in ky).
  #'   \item \emph{climate.rainfallAverage_0.2-1} numeric, time delayed values of rainfallAverage for lags 0.2 to 1 (in ky).
}
}
\author{
Blas M. Benito  <blasbenito@gmail.com>
}
\keyword{datasets}
