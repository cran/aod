\name{mice}

\alias{mice}

\docType{data}

\title{ Pregnant Female Mice Experiment }

\description{
  Unpublished laboratory data on the proportion of affected foetuses in two groups (control and treatment) of 10 
  pregnant female mice (Kupper and Haseman, 1978, p. 75).
  }

\usage{data(mice)}

\format{
  A data frame with 20 observations on the following 3 variables.
  \describe{
    \item{group}{a factor with levels \code{CTRL} and \code{TREAT}}
    \item{n}{a numeric vector: the total number of foetuses.}
    \item{y}{a numeric vector: the number of affected foetuses.}
    }
  }

\references{
  Kupper, L.L., Haseman, J.K., 1978. \emph{The use of a correlated binomial model for the analysis of a certain 
  toxicological experiments.} Biometrics 34, 69-76.
  }

\keyword{datasets}
