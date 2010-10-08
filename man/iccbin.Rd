\name{iccbin}
\alias{iccbin}
\alias{show,iccbin-method}
\title{Intra-Cluster Correlation for Binomial Data}

\description{This function calculates point estimates of the intraclass correlation \eqn{\rho}
from clustered binomial data \eqn{{(n_1, y_1), (n_2, y_2), ..., (n_K, y_K)}} (with \eqn{K} the number of clusters),
using a 1-way random effect model. Three estimates, following methods referred to as \dQuote{A}, \dQuote{B} and \dQuote{C}
in Goldstein et al. (2002), can be obtained.
}

\usage{iccbin(n, y, data, method = c("A", "B", "C"), nAGQ = 1, M = 1000)}
\arguments{
  \item{n}{Vector of the denominators of the proportions.}
  \item{y}{Vector of the numerators of the proportions.}
  \item{data}{A data frame containing the variables n and y.}
  \item{method}{A character (\dQuote{A}, \dQuote{B} or \dQuote{C}) defining the calculation method. See Details.}
  \item{nAGQ}{Same as in function \code{glmer} of package
    \pkg{lme4}. Only for methods \dQuote{A} and \dQuote{B}. Default to 1.}
  \item{M}{Number of Monte Carlo (MC) replicates used in method \dQuote{B}. Default to 1000.}
  }

\details{
  Before computations, the clustered data are split into binary (0/1) observations \eqn{y_{ij}} (obs. \eqn{j} in cluster
  \eqn{i}).
  The calculation methods are described in Goldstein et al. (2002).
  Methods "A" and "B" assume a 1-way generalized linear mixed model,
  and method "C" a 1-way linear mixed model.\cr
  For "A" and "B", function \code{iccbin} uses the logistic binomial-Gaussian model:
  \deqn{y_{ij} | p_{ij} \sim Bernoulli(p_{ij}),}
  \deqn{logit(p_{ij}) = b_0 + u_i,}
  where \eqn{b_0} is a constant and \eqn{u_i} a cluster random effect with \eqn{u_i \sim N(0, s^2_u)}.
  The ML estimate of the variance component \eqn{s^2_u} is calculated with the function \code{glmer} of package \pkg{lme4}.
  The intra-class correlation \eqn{\rho = Corr[y_{ij}, y_{ij'}]} is then calculated with a first-order model linearization
  around \eqn{E[u_i]=0} in method \dQuote{A}, and with Monte Carlo simulations in method \dQuote{B}.\cr
  For \dQuote{C}, function \code{iccbin} provides the common ANOVA (moments) estimate of \eqn{\rho}.
  For details, see for instance Donner (1986), Searle et al. (1992) or Ukoumunne (2002).
  }

\value{
  An object of formal class \dQuote{iccbin}, with 3 slots:
  \item{CALL}{The call of the function.}
  \item{features}{A character vector summarizing the main features of the method used.}
  \item{rho}{The point estimate of the intraclass correlation \eqn{\rho}.}
  }

\references{
  Donner A., 1986, \emph{A review of inference procedures for the intraclass correlation coefficient in the one-way random
  effects model}. International Statistical Review 54, 67-82.\cr
  Searle, S.R., Casella, G., McCulloch, C.E., 1992. \emph{Variance components}. Wiley, New York.\cr
  Ukoumunne, O. C., 2002. \emph{A comparison of confidence interval methods for the intraclass
  correlation coefficient in cluster randomized trials}. Statistics in Medicine 21, 3757-3774.\cr
  Golstein, H., Browne, H., Rasbash, J., 2002. \emph{Partitioning variation in multilevel models}.
  Understanding Statistics 1(4), 223-231.
  }

\author{Matthieu Lesnoff \email{matthieu.lesnoff@cirad.fr}, Renaud Lancelot \email{renaud.lancelot@cirad.fr}}

\seealso{\code{\link{iccbin-class}}, \code{\link[lme4:glmer]{glmer}}\cr}

\examples{

data(rats)
tmp <- rats[rats$group == "TREAT", ]
# A: glmm (model linearization)
iccbin(n, y, data = tmp, method = "A")
iccbin(n, y, data = tmp, method = "A", nAGQ = 10)
# B: glmm (Monte Carlo)
iccbin(n, y, data = tmp, method = "B")
iccbin(n, y, data = tmp, method = "B", nAGQ = 10, M = 1500)
# C: lmm (ANOVA moments)
iccbin(n, y, data = tmp, method = "C")

  \dontrun{
  # Example of confidence interval calculation with nonparametric bootstrap
  require(boot)
  foo <- function(X, ind) {
    n <- X$n[ind]
    y <- X$y[ind]
    X <- data.frame(n = n, y = y)
    iccbin(n = n, y = y, data = X, method = "C")@rho[1]
    }
  res <- boot(data = tmp[, c("n", "y")], statistic = foo, R = 500, sim = "ordinary", stype = "i")
  res
  boot.ci(res, conf = 0.95, type = "basic")
  }

}

\keyword{htest}
