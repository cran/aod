\name{logLik-methods}
\docType{methods}
%\alias{logLik-methods}
%\alias{logLik,ANY-method}
\alias{logLik,glimML-method}

\title{Methods for Functions "logLik" in Package "aod"}

\description{Extracts the maximized log-likelihood from fitted models of formal class \dQuote{glimML}.}

\usage{\S4method{logLik}{glimML}(object, \dots)}

\arguments{
  \item{object}{A fitted model of formal class \dQuote{glimML} (functions \code{betabin} or \code{negbin}).}
  \item{\dots}{Other arguments passed to methods.}
  }

\section{Methods}{\describe{
  \item{ANY}{Generic function: see \code{\link[stats]{logLik}}.}
  \item{glimML}{Extract the maximized log-likelihood from models of formal class \dQuote{glimML}, fitted by functions 
    \code{betabin} and \code{negbin}.}
  }}

\value{
  A numeric scalar with 2 attributes: \dQuote(df) (number of parameters in the model) and \dQuote{nobs} 
  (number of observations = degrees of freedom of the residuals + number of parameters in the model).
  }

\seealso{\code{\link[stats:logLik]{logLik}} in package \pkg{stats}.}

\keyword{methods}
