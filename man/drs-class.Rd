\name{drs-class}
\docType{class}
\alias{drs-class}
\alias{show,drs-method}
\title{Representation of Objects of Formal Class "drs"}
\description{Representation of the output of functions \code{donner} and \code{raoscott}.}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("drs", ...)} or, more commonly, 
  via the \code{donner} or \code{raoscott} functions.}
\section{Slots}{%
  \describe{%
    \item{\code{CALL}}{The call of the function.}
    \item{\code{tab}}{A data frame containing test information. The content of the data frame depends on the type 
      of the function which generated it.}
    \item{\code{rho}}{The ANOVA estimate of the intra-cluster correlation (function \code{donner}).}
    \item{\code{X2}}{The adjusted \eqn{\chi^2}{chi-squared} statistic.}
    }
  }
\section{Methods}{%
  \describe{
    \item{donner}{\code{signature(object = "drs")}: see \code{\link{donner}}.}
    \item{raoscott}{\code{signature(object = "drs")}: see \code{\link{raoscott}}.}
    }
  }
\keyword{classes}
