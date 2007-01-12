\name{betabin}
\alias{betabin}
\title{Beta-Binomial Model for Proportions}

\description{
  Fits a beta-binomial generalized linear model accounting for overdispersion in clustered binomial data \eqn{(n, y)}.
  }

\usage{
  betabin(formula, random, data, link = c("logit", "cloglog"),
          phi.ini = NULL, warnings = FALSE, na.action = na.omit, 
          fixpar = list(), hessian = TRUE, \dots)
  }

\arguments{
  \item{formula}{A formula for the fixed effects \code{b}. The left-hand side of the formula must be of the form
    \code{cbind(y, n - y)} where the modelled probability is \code{y/n}.}
  \item{random}{A right-hand formula for the overdispersion parameter(s) \eqn{\phi}.}
  \item{link}{The link function for the mean \eqn{p}: \dQuote{logit} or \dQuote{cloglog}.}
  \item{data}{A data frame containing the response (\code{n} and \code{y}) and explanatory variable(s).}
  \item{phi.ini}{Initial values for the overdispersion parameter(s) \eqn{\phi}. Default to 0.1.}
  \item{warnings}{Logical to control the printing of warnings occurring during log-likelihood maximization. 
    Default to \code{FALSE} (no printing).}
  \item{na.action}{A function name: which action should be taken in the case of missing value(s).}
  \item{fixpar}{A list with 2 components (scalars or vectors) of the same size, indicating which parameters are 
    fixed (i.e., not optimized) in the global parameter vector \eqn{(b, \phi)} and the corresponding fixed values.\cr
    For example, \code{fixpar = list(c(4, 5), c(0, 0))} means that 4th and 5th parameters of the model are set to 0.}
  \item{hessian}{A logical. When set to \code{FALSE}, the hessian and the variances-covariances matrices of the 
    parameters are not computed.}
  \item{\dots}{Further arguments passed to \code{\link{optim}}.}
  }

\details{For a given cluster \eqn{(n, y)}, the model is:
  \deqn{y~|~\lambda \sim Binomial(n,~\lambda)}{%
    y | \lambda ~ Binomial(n, \lambda)}
  with \eqn{\lambda} following a Beta distribution \eqn{Beta(a1,~a2)}{Beta(a1, a2)}.\cr
  If \eqn{B} denotes the beta function, then:
   \deqn{P(\lambda) = \frac{\lambda^{a1~-~1} * (1~-~\lambda)^{a2 - 1}}{B(a1,~a2)}}{%
     P(\lambda) = \lambda^{a1 - 1} * (1 - \lambda)^{a2 - 1} / B(a1, a2)}
   \deqn{E[\lambda] = \frac{a1}{a1 + a2}}{E[\lambda] = a1 / (a1 + a2)}
   \deqn{Var[\lambda] = \frac{a1 * a2}{(a1 + a2 + 1) * (a1 + a2)^2}}{%
     Var[\lambda] = a1 * a2 / [(a1 + a2 + 1) * (a1 + a2)^2]}
  The marginal beta-binomial distribution is:
  \deqn{P(y) = \frac{C(n,~y) * B(a1 + y, a2 + n - y)}{B(a1,~a2)}}{%
    P(y) = C(n, y) * B(a1 + y, a2 + n - y) / B(a1, a2)}
  The function uses the parameterization \eqn{p = \frac{a1}{a1 + a2} = h(X b) = h(\eta)}{%
    p = a1 / (a1 + a2) = h(X b) = h(\eta)} and \eqn{\phi = \frac{1}{a1 + a2 + 1}}{\phi = 1 / (a1 + a2 + 1)}, 
  where \eqn{h} is the inverse of the link function (logit or complementary log-log), \eqn{X} is a design-matrix, \eqn{b} 
  is a vector of fixed effects, \eqn{\eta = X b} is the linear predictor and \eqn{\phi} is the overdispersion
  parameter (i.e., the intracluster correlation coefficient, which is here restricted to be positive).\cr
  The marginal mean and variance are:
  \deqn{E[y] = n * p}
  \deqn{Var[y] = n * p * (1 - p) * [1 + (n - 1) * \phi]}
  The parameters \eqn{b} and \eqn{\phi} are estimated by maximizing the log-likelihood of the marginal model (using the
  function \code{optim}). Several explanatory variables are allowed in \eqn{b}, only one in \eqn{\phi}.
  }

\value{An object of formal class \dQuote{glimML}: see \code{\link{glimML-class}} for details.
  }

\references{
  Crowder, M.J., 1978. \emph{Beta-binomial anova for proportions}. Appl. Statist. 27, 34-37.\cr
  Griffiths, D.A., 1973. \emph{Maximum likelihood estimation for the beta-binomial distribution and an application
  to the household distribution of the total number of cases of disease}. Biometrics 29, 637-648.\cr
  Prentice, R.L., 1986. \emph{Binary regression using an extended beta-binomial distribution, with discussion of
  correlation induced by covariate measurement errors}. J.A.S.A. 81, 321-327.\cr
  Williams, D.A., 1975. \emph{The analysis of binary responses from toxicological experiments involving
  reproduction and teratogenicity}. Biometrics 31, 949-952.
  }

\author{Matthieu Lesnoff \email{matthieu.lesnoff@cirad.fr}, Renaud Lancelot \email{renaud.lancelot@cirad.fr}
  }
  
\seealso{ \code{\link{glimML-class}}, \code{\link{glm}} and \code{\link{optim}}
  }
  
\examples{
  data(orob2)
  fm1 <- betabin(cbind(y, n - y) ~ seed, ~ 1, data = orob2)
  fm2 <- betabin(cbind(y, n - y) ~ seed + root, ~ 1, data = orob2)
  fm3 <- betabin(cbind(y, n - y) ~ seed * root, ~ 1, data = orob2)
  # show the model
  fm1; fm2; fm3
  # AIC
  AIC(fm1, fm2, fm3)
  summary(AIC(fm1, fm2, fm3), which = "AICc")
  # Wald test for root effect
  wald.test(b = coef(fm3), Sigma = vcov(fm3), Terms = 3:4)
  # likelihood ratio test for root effect
  anova(fm1, fm3)
  # model predictions
  New <- expand.grid(seed = levels(orob2$seed),
                     root = levels(orob2$root))
  data.frame(New, predict(fm3, New, se = TRUE, type = "response"))
  # Djallonke sheep data
  data(dja)
  betabin(cbind(y, n - y) ~ group, ~ 1, dja)
  # heterogeneous phi
  betabin(cbind(y, n - y) ~ group, ~ group, dja,
          control = list(maxit = 1000))
  # phi fixed to zero in group TREAT
   betabin(cbind(y, n - y) ~ group, ~ group, dja,
    fixpar = list(4, 0))
  # glim without overdispersion
  summary(glm(cbind(y, n - y) ~ group,
    family = binomial, data = dja))
  # phi fixed to zero in both groups
  betabin(cbind(y, n - y) ~ group, ~ group, dja,
    fixpar = list(c(3, 4), c(0, 0))) 
  }
\keyword{regression}
