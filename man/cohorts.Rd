\name{cohorts}
\alias{cohorts}
\docType{data}
\title{Age, Period and Cohort Effects for Vital Rates}
\description{Number of prostate cancer deaths and midperiod population for nonwhites in the USA by age and period.
  The cohort index \eqn{k} is related to age and period indices (\eqn{i} and \eqn{j}, respectively) by \eqn{k = j + I - i}, 
  where \eqn{I = max(i)} (Holford, 1983, Table 2).}
\usage{data(cohorts)}
\format{A data frame with 49 observations on the following 4 variables.
  \describe{
    \item{period}{A factor with levels \code{1935-}, \code{1940-}, \ldots, \code{1965-}.}
    \item{age}{A factor with levels \code{50-}, \code{55-}, \ldots, \code{80-}.}
    \item{y}{Numeric: the number of prostate cancer deaths.}
    \item{n}{Numeric: the midperiod population size.}
    }
  }
\references{
  Holford, T.R., 1983. \emph{The estimation of age, period and cohort effects for vital rates}. Biometrics 39, 311-324.
  }
\keyword{datasets}
