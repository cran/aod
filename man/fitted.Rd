\name{fitted-methods}
\docType{methods}
\alias{fitted-methods}
\alias{fitted,ANY-method}
\alias{fitted,glimML-method}
\alias{fitted,glimQL-method}
\title{Methods for Function "fitted" in Package "aod"}

\description{Extracts the fitted values from models.}

\section{Methods}{\describe{
  \item{ANY}{Generic function: see \code{\link[stats]{fitted}}.}
  \item{glimML}{Extract the fitted values from models of formal class \dQuote{glimML}, presently generated 
    by functions \code{betabin} (maximum-likelihood beta-binomial regression) and \code{negbin} (maximum-likelihood 
    negative-binomial regression).} 
  \item{glimQL}{Extract the fitted values from models of formal class \dQuote{glimQL}, presently generated 
    by functions \code{quasibin} (quasi-likelihood binomial regression) and \code{quasibin} (quasi-likelihood Poisson
    regression).} 
  }}

\keyword{methods}
