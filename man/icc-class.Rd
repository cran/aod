\name{icc-class}
\docType{class}
\alias{icc-class}
\title{Representation of Objects of Formal Class "icc"}
\description{Representation of the output of function \code{icc}.}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("icc", ...)} or, more commonly,
  via the function \code{icc}.}}
\section{Slots}{
  \describe{
    \item{\code{CALL}}{The call of the function.}
    \item{\code{rho}}{A numeric scalar giving the intra-cluster correlation.}
    \item{\code{f}}{A named numeric vector giving all the components of the \eqn{F} test for \eqn{\rho}.}
    \item{\code{rho.MC}}{A numeric scalar giving the Monte Carlo estimate for \eqn{\rho}.}
    \item{\code{varrho}}{A 2-component numeric vector giving the delta-method and anova estimates for the variance 
                         of \eqn{\rho}.}
    \item{\code{method}}{A character vector indicating which methods were used to estimate the \eqn{\rho}:
                         ML/REML and ANOVA.}
    \item{\code{features}}{A named numeric vector summarizing the design.}
    }
  }
\section{Methods}{
  \describe{
    \item{icc}{\code{signature(object = "icc")}: see \code{\link{icc}}.}
    }
  }
\keyword{classes}
