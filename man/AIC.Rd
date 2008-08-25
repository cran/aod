\name{AIC-methods}
\docType{methods}
\alias{AIC,glimML-method}
\title{Akaike Information Criteria}
\description{Extracts the Akaike information criterion (AIC) and the corrected AIC (AICc) from fitted models of formal 
  class \dQuote{glimML} and possibly computes derived statistics.}

\usage{\S4method{AIC}{glimML}(object, \dots, k = 2)}

\arguments{
  \item{object}{fitted model of formal class \dQuote{glimML} (functions \code{betabin} or \code{negbin}).}
  \item{\dots}{optional list of fitted models separated by commas.}
  \item{k}{numeric scalar, with a default value set to 2, thus providing the regular AIC.}
  }

\section{Methods}{
  \describe{
    \item{glimML}{Extracts the AIC and AICc from models of formal class \dQuote{glimML}, fitted by functions 
    \code{betabin} and \code{negbin}.}
    }
  }
  
\details{
  \eqn{AIC = -2~\mbox{log-likelihood} + 2*n_{par}}{-2 * log-likelihood + 2 * npar}, where \eqn{n_{par}}{npar}
  represents the number of parameters in the fitted model.\cr
  \eqn{AICc = AIC + 2 * n_{par} * (n_{par} + 1) / (n_{obs} - n_{par} + 1)}{AICc = AIC + 2 * npar * (npar + 1) / (nobs - npar + 1)},
  where \eqn{n_{obs}}{nobs} is the number of observations used to compute the log-likelihood. It should be used when the number 
  of fitted parameters is large compared to sample size, i.e., when  \eqn{n_{obs} / n_{par} < 40}{nobs / npar < 40} (Hurvich and Tsai, 1995).
  }
  
\seealso{Examples in \code{\link{betabin}} and see \code{\link[stats:AIC]{AIC}} in package \pkg{stats}.}

\references{
  Burnham, K.P., Anderson, D.R., 2002. \emph{Model selection and multimodel inference: a practical
  information-theoretic approach}. New-York, Springer-Verlag, 496 p.\cr
  Hurvich, C.M., Tsai, C.-L., 1995. \emph{Model selection for extended quasi-likelihood models in small samples}. 
  Biometrics, 51 (3): 1077-1084.}

\keyword{methods}
