\name{iccbin-class}
\docType{class}
\alias{iccbin-class}
\title{Representation of Objects of Formal Class "iccbin"}
\description{Representation of the output of function \code{iccbin}.}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("iccbin", ...)} or, more commonly,
  via the function \code{iccbin}.}
\section{Slots}{
  \describe{
    \item{\code{CALL}}{The call of the function.}
    \item{\code{features}}{A character vector summarizing the main features of the method used.}
    \item{\code{rho}}{A numeric scalar giving the intra-cluster correlation.}
    }
  }
\section{Methods}{
  \describe{
    \item{icc}{\code{signature(object = "iccbin")}: see \code{\link{iccbin}}.}
    }
  }
\keyword{classes}
