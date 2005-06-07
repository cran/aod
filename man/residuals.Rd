\name{residuals-methods}
\docType{methods}
\alias{residuals-methods}
\alias{residuals,ANY-method}
\alias{residuals,glimML-method}
\alias{residuals,glimQL-method}

\title{Residuals for Maximum-Likelihood and Quasi-Likelihood Models}

\description{
  Residuals of models fitted with functions \code{betabin} and \code{negbin} (formal class \dQuote{glimML}), or 
  \code{quasibin} and \code{quasipois} (formal class \dQuote{glimQL}).
  }

\usage{
  \S4method{residuals}{glimML}(object, type = c("pearson", "response", "link"), \dots)
  \S4method{residuals}{glimQL}(object, type = c("pearson", "response", "link"), \dots)
  }
\arguments{
  \item{object}{Fitted model of formal class \dQuote{glimML} or \dQuote{glimQL}.}
  \item{type}{Character string for the type of residual: \dQuote{pearson} (default) or \dQuote{response}.}
  \item{\dots}{Further arguments to be passed to the function, such as \code{na.action}.}
  }

\details{
  For models fitted with \code{betabin} or \code{quasibin}, Pearson's residuals are computed as:
  \deqn{\frac{y - n * \hat{p}}{\sqrt{n * \hat{p} * (1 - \hat{p}) * (1 + (n - 1) * \hat{\phi})}}}{%
    (y - n * p.fit) / (n * p.fit * (1 - p.fit) * (1 + (n - 1) * \phi))^{0.5}}
  where \eqn{y} and \eqn{n} are respectively the numerator and the denominator of the response, \eqn{\hat{p}}{p.fit} 
  is the fitted probability and \eqn{\hat{\phi}}{\phi} is the fitted overdispersion parameter. When \eqn{n = 0}, the 
  residual is set to 0. Response residuals are computed as \eqn{y/n - \hat{p}}{y/n - p.fit}.\cr
  For models fitted with \code{negbin} or \code{quasipois}, Pearson's residuals are computed as:
  \deqn{\frac{y - \hat{y}}{\sqrt{\hat{y} + \hat{\phi} * \hat{y}^2}}}{%
    (y - y.fit) / (y.fit + \phi * y.fit^2)^{0.5}}
  where \eqn{y} and \eqn{\hat{y}}{y.fit} are the observed and fitted counts, respectively. Response residuals are 
  computed as \eqn{y - \hat{y}}{y - y.fit}.
  }

\value{A numeric vector of residuals.}

\author{Matthieu Lesnoff \email{matthieu.lesnoff@cirad.fr}, Renaud Lancelot \email{renaud.lancelot@cirad.fr}}

\seealso{ \code{\link[stats:residuals.glm]{residuals.glm}} }

\examples{
  data(orob2)
  fm <- betabin(cbind(y, n - y) ~ seed, ~ 1,
                 link = "logit", data = orob2)
  #Pearson's chi-squared goodness-of-fit statistic
  sum(residuals(fm, type = "pearson")^2)
  }

\keyword{regression}
