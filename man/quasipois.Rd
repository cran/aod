\name{quasipois}
\alias{quasipois}
\title{Quasi-Likelihood Model for Counts}
\description{The function fits the log linear model (\dQuote{Procedure II}) proposed by Breslow (1984) accounting for 
  overdispersion in counts \eqn{y}.}
\usage{quasipois(formula, data, phi = NULL, tol = 0.001)}
\arguments{
  \item{formula}{A formula for the fixed effects. The left-hand side of the formula must be the counts \code{y} i.e.,
    positive integers (\code{y >= 0}). The right-hand side can involve an offset term.}
  \item{data}{A data frame containing the response (\code{y}) and explanatory variable(s).}
  \item{phi}{When \code{phi} is NULL (the default), the overdispersion parameter \eqn{\phi} is estimated from the data.
    Otherwise, its value is considered as fixed.}
  \item{tol}{A positive scalar (default to 0.001). The algorithm stops at iteration \eqn{r + 1} when the condition
  \eqn{\chi{^2}[r+1] - \chi{^2}[r] <= tol} is met by the \eqn{\chi^2}{chi-squared} statistics .}
  }
\details{
  For a given count \eqn{y}, the model is:
  \deqn{y~|~\lambda \sim Poisson(~\lambda)}{y | \lambda ~ Poisson(\lambda)}
  with \eqn{\lambda} a random variable of mean \eqn{E[\lambda] = \mu}
  and variance \eqn{Var[\lambda] = \phi * \mu^2}.\cr
  The marginal mean and variance are:
  \deqn{E[y] = \mu}
  \deqn{Var[y] = \mu + \phi * \mu^2}
  The function uses the function \code{glm} and the recrystallization: \eqn{\mu = exp(X b) = exp(\eta)}, where \eqn{X} 
  is a design-matrix, \eqn{b} is a vector of fixed effects and \eqn{\eta = X b} is the linear predictor. \cr
  The estimate of \eqn{b} maximizes the quasi log-likelihood of the marginal model.
  The parameter \eqn{\phi} is estimated with the moment method or can be set to a constant
  (a regular \emph{glim} is fitted when \eqn{\phi} is set to 0). The literature recommends to estimate \eqn{\phi}
  with the saturated model. Several explanatory variables are allowed in \eqn{b}. None is allowed in \eqn{\phi}.\cr
  An offset can be specified in the argument \code{formula} to model rates \eqn{y/T} (see examples). The offset and the
  marginal mean are \eqn{log(T)} and \eqn{\mu = exp(log(T) + \eta)}, respectively.
  }
\value{An object of formal class \dQuote{glimQL}: see \code{\link{glimQL-class}} for details.}
\references{
  Breslow, N.E., 1984. \emph{Extra-Poisson variation in log-linear models}. Appl. Statist. 33, 38-44.\cr
  Moore, D.F., Tsiatis, A., 1991. \emph{Robust estimation of the variance in moment methods for extra-binomial
  and extra-poisson variation}. Biometrics 47, 383-401.
  }
\author{Matthieu Lesnoff \email{matthieu.lesnoff@cirad.fr}, Renaud Lancelot \email{renaud.lancelot@cirad.fr}}
\seealso{\code{\link{glm}}, \code{\link[MASS]{negative.binomial}} in the recommended package \pkg{MASS}, 
  \code{\link[geepack]{geese}} in the contributed package \pkg{geepack}, 
  \code{\link[dispmod]{glm.poisson.disp}} in the contributed package \pkg{dispmod}.}
\examples{
  # without offset
  data(salmonella)
  quasipois(y ~ log(dose + 10) + dose,
            data = salmonella)
  quasipois(y ~ log(dose + 10) + dose, 
            data = salmonella, phi = 0.07180449)
  summary(glm(y ~ log(dose + 10) + dose,
          family = poisson, data = salmonella))
  quasipois(y ~ log(dose + 10) + dose,
          data = salmonella, phi = 0)
  # with offset
  data(cohorts)
  i <- cohorts$age ; levels(i) <- 1:7
  j <- cohorts$period ; levels(j) <- 1:7
  i <- as.numeric(i); j <- as.numeric(j)
  cohorts$cohort <- j + max(i) - i
  cohorts$cohort <- as.factor(1850 + 5 * cohorts$cohort)
  fm1 <- quasipois(y ~ age + period + cohort + offset(log(n)),
                   data = cohorts)
  fm1
  quasipois(y ~ age + cohort + offset(log(n)),
            data = cohorts, phi = fm1@phi)
  }
\keyword{regression}
