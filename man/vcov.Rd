\name{vcov-methods}
\docType{methods}
%\alias{vcov-methods}
%\alias{vcov,ANY-method}
\alias{vcov,glimML-method}
\alias{vcov,glimQL-method}
\alias{geese-class}
\alias{geeglm-class}
\alias{vcov,geese-method}
\alias{vcov,geeglm-method}
\title{Methods for Function "vcov" in Package "aod"}

\description{Extract the approximate var-cov matrix of estimated coefficients from fitted models.}

\section{Methods}{\describe{
  \item{ANY}{Generic function: see \code{\link[stats]{vcov}}.}
  \item{glimML}{Extract the var-cov matrix of estimated coefficients for fitted models of formal class \dQuote{glimML}.}
  \item{glimQL}{Extract the var-cov matrix of estimated coefficients for fitted models of formal class \dQuote{glimQL}.}
  \item{geese}{Extract the var-cov matrix of estimated coefficients for fitted models of class \dQuote{geese} 
               (contributed package \pkg{geepack}).}
  \item{geeglm}{Extract the var-cov matrix of estimated coefficients for fitted objects of class \dQuote{geeglm} 
                (contributed package \pkg{geepack}).}
  }}

\keyword{methods}
