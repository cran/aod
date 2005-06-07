\name{varbin-class}
\docType{class}
\alias{varbin-class}
\alias{show,varbin-method}
\title{Representation of Objects of Formal Class "varbin"}
\description{Representation of the output of function \code{varbin} used to estimate proportions and their variance 
             under various distribution assumptions.}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("varbin", ...)} or, more commonly,
  via the function \code{varbin}.}}
\section{Slots}{
  \describe{
    \item{\code{CALL}}{The call of the function.}
    \item{\code{tab}}{A data frame containing the estimates, their variance and the confidence limits.}
    \item{\code{pboot}}{A numeric vector containing the bootstrap replicates.}
    \item{\code{alpha}}{The \eqn{\alpha} level to compute confidence intervals.}
    \item{\code{features}}{A named numeric vector summarizing the design.}
    }
  }
\section{Methods}{
  \describe{
    \item{varbin}{\code{signature(object = "varbin")}: see \code{\link{varbin}}.}
    }
  }
\keyword{classes}
