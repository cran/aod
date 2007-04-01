\name{negbin}
\alias{negbin}
\title{Negative-Binomial Model for Counts}
\description{The function fits a negative-binomial log linear model accounting for overdispersion in counts \eqn{y}.
  }

\usage{
  negbin(formula, random, data, phi.ini = NULL, warnings = FALSE, 
         na.action = na.omit, fixpar = list(),
         hessian = TRUE, control = list(maxit = 2000), \dots)
  }

\arguments{
  \item{formula}{A formula for the fixed effects. The left-hand side of the formula must be the counts \code{y} i.e., 
    positive integers (\code{y >= 0}). The right-hand side can involve an offset term.}
  \item{random}{A right-hand formula for the overdispersion parameter(s) \eqn{\phi}.}
  \item{data}{A data frame containing the response (\code{y}) and explanatory variable(s).}
  \item{phi.ini}{Initial values for the overdispersion parameter(s) \eqn{\phi}. Default to 0.1.}
  \item{warnings}{Logical to control printing of warnings occurring during log-likelihood maximization. Default to 
    FALSE (no printing).}
  \item{na.action}{A function name. Indicates which action should be taken in the case of missing value(s).}
  \item{fixpar}{A list with 2 components (scalars or vectors) of the same size, indicating which parameters are 
    fixed (i.e., not optimized) in the global parameter vector \eqn{(b, \phi)} and the corresponding fixed values.\cr
    For example, \code{fixpar = list(c(4, 5), c(0, 0))} means that 4th and 5th parameters of the model are set to 0.}
  \item{hessian}{A logical. When set to \code{FALSE}, the hessian and the variances-covariances matrices of the 
    parameters are not computed.}
  \item{control}{A list to control the optimization parameters. See \code{\link{optim}}. By default, set the maximum number of iterations to 2000.}
  \item{\dots}{Further arguments passed to \code{\link{optim}}.}
    }

\details{For a given count \eqn{y}, the model is:
  \deqn{y~|~\lambda \sim Poisson(~\lambda)}{y | \lambda ~ Poisson(\lambda)}
  with \eqn{\lambda} following a Gamma distribution \eqn{Gamma(r,~\theta)}{Gamma(r, \theta)}.\cr
  If \eqn{G} denote the gamma function, then:
  \deqn{P(\lambda) = r^{-\theta} * \lambda^{\theta - 1} * \frac{exp(-\frac{\lambda}{r})}{G(\theta)}}{%
    P(\lambda) = r^{-\theta} * \lambda^{\theta - 1} * exp(-\lambda / r) / G(\theta)}
  \deqn{E[\lambda] = \theta * r}
  \deqn{Var[\lambda] = \theta * r^2}
  The marginal negative-binomial distribution is:
  \deqn{P(y) = G(y + \theta) * \left(\frac{1}{1 + r}\right)^\theta * \frac{(\frac{r}{1 + r})^y}{y! * G(\theta)} }{%
    P(y) = G(y + \theta) * (1 / (1 + r))^\theta * (r / (1 + r))^y / (y! * G(\theta)) }
  The function uses the parameterization \eqn{\mu = \theta * r = exp(X b) = exp(\eta)} and \eqn{\phi = 1 / \theta}, 
  where \eqn{X} is a design-matrix, \eqn{b} is a vector of fixed effects, \eqn{\eta = X b} is the linear predictor and 
  \eqn{\phi} the overdispersion parameter.\cr
  The marginal mean and variance are:
  \deqn{E[y] = \mu}
  \deqn{Var[y] = \mu + \phi * \mu^2}
  The parameters \eqn{b} and \eqn{\phi} are estimated by maximizing the log-likelihood of the marginal model (using the 
  function \code{optim()}). Several explanatory variables are allowed in \eqn{b}. Only one is allowed in \eqn{\phi}.\cr
  An offset can be specified in the \code{formula} argument to model rates \eqn{y/T}. The offset and the marginal mean
  are \eqn{log(T)} and \eqn{\mu = exp(log(T) + \eta)}, respectively.
  }
  
\value{An object of formal class \dQuote{glimML}: see \code{\link{glimML-class}} for details.
  }
  
\references{
  Lawless, J.F., 1987. \emph{Negative binomial and mixed Poisson regression}.
  The Canadian Journal of Statistics, 15(3): 209-225.
  }
  
\author{Matthieu Lesnoff \email{matthieu.lesnoff@cirad.fr}, Renaud Lancelot \email{renaud.lancelot@cirad.fr}
  }
\seealso{\code{\link{glimML-class}}, \code{\link{glm}} and \code{\link{optim}},\cr
         \code{\link[MASS:glm.nb]{glm.nb}} in the recommended package \pkg{MASS},\cr
         \code{gnlr} in package \pkg{gnlm} available at \url{www.luc.ac.be/~jlindsey/rcode.html}.
  }
  
\examples{
  # without offset
  data(salmonella)
  negbin(y ~ log(dose + 10) + dose, ~ 1, salmonella)
  library(MASS) # function glm.nb in MASS
  fm.nb <- glm.nb(y ~ log(dose + 10) + dose,
                  link = log, data = salmonella)
  coef(fm.nb)
  1 / fm.nb$theta # theta = 1 / phi
  c(logLik(fm.nb), AIC(fm.nb))
  # with offset
  data(dja)
  negbin(y ~ group + offset(log(trisk)), ~ group, dja)
  # phi fixed to zero in group TREAT
  negbin(y ~ group + offset(log(trisk)), ~ group, dja,
    fixpar = list(4, 0))
  # glim without overdispersion
  summary(glm(y ~ group + offset(log(trisk)),
    family = poisson, data = dja))
  # phi fixed to zero in both groups
  negbin(y ~ group + offset(log(trisk)), ~ group, dja,
    fixpar = list(c(3, 4), c(0, 0))) 
  }

\keyword{regression}
