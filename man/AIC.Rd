\name{AIC-methods}
\docType{methods}
\alias{AIC-methods}
\alias{AIC,ANY-method}
\alias{AIC,logLik-method}
\alias{AIC,glimML-method}

\title{Methods for Function "AIC" in Package "aod"}

\description{Extracts the Akaike information criterion (AIC) from fitted models of formal class \dQuote{glimML}.}

\usage{\S4method{AIC}{glimML}(object, \dots, k = 2)}

\arguments{
  \item{object}{A fitted model of formal class \dQuote{glimML} (functions \code{betabin} or \code{negbin}).}
  \item{\dots}{An optional list of fitted models separated by commas.}
  \item{k}{A numeric scalar, with a default value set to 2, thus providing the regular AIC.}
  }

\section{Methods}{\describe{
  \item{ANY}{Generic function: see \code{\link[stats]{AIC}}.}
  \item{glimML}{Extract the AIC from models of formal class \dQuote{glimML}, fitted by functions 
    \code{betabin} and \code{negbin}.
    The AIC is defined as \eqn{-2 \mbox{log-likelihood} + 2*n_{par}}{-2*log-likelihood + 2*npar}, were \eqn{n_{par}}{npar} 
    represents the number of parameters in the fitted model.}
  }}

\seealso{Check the examples in \code{\link{betabin}} and see \code{\link[stats:AIC]{AIC}} in package \pkg{stats}.}

\references{
  Burnham, K.P., Anderson, D.R., 2002. \emph{Model selection and multimodel inference: a practical
  information-theoretic approach}. New-York, Springer-Verlag, 496 p.
  }

\keyword{methods}
