\name{donner}
\alias{donner}
\alias{show,donner-class}
\title{Test of Proportion Homogeneity using Donner's Adjustment}
\description{Tests the homogeneity of proportions between \eqn{I} groups (H0: \eqn{p_1 = p_2 = ... = p_I} ) from clustered
  binomial data \eqn{(n, y)} using the adjusted \eqn{\chi^2}{chi-squared} statistic proposed by Donner (1989).}
\usage{donner(formula = NULL, response = NULL,
       weights = NULL, group = NULL, data, C = NULL)}
\arguments{
  \item{formula}{An optional formula where the left-hand side is either a matrix of the form \code{cbind(y, n-y)}, where 
    the modelled probability is \code{y/n}, or a vector of proportions to be modelled (\code{y/n}). In both cases, the 
    right-hand side must specify a single grouping variable. When the left-hand side of the formula is a vector of 
    proportions, the argument \code{weight} must be used to indicate the denominators of the proportions.}
  \item{response}{An optional argument indicating either a matrix of the form \code{cbind(y, n-y)}, where the modelled 
    probability is \code{y/n}, or a vector of proportions to be modelled (\code{y/n}).}
  \item{weights}{An optional argument used when the left-hand side of \code{formula} or \code{response} is a vector 
    of proportions: \code{weight} is the denominator of the proportion.}
  \item{group}{An optional argument only used when \code{response} is used. In this case, this argument is a factor 
    indicating a grouping variable.}
  \item{data}{A data frame containing the response (\code{n} and \code{y}) and the grouping variable.}
  \item{C}{If not NULL, a numerical vector of \eqn{I} cluster correction factors.}}
\details{
  The \eqn{\chi^2}{chi-squared} statistic is adjusted with the correction factor \eqn{C_i}
  computed in each group \eqn{i}. The test statistic is given by:
  \deqn{X^2 = \sum_{i}\frac{(y_i - n_i * p)^2}{C_i * n_i * p * (1 - p)}}{%
    X^2 = sum( (y_i - n_i * p)^2 / (C_i * n_i * p * (1 - p)) )}
  where \eqn{C_i = 1 + (nA_i - 1) * \rho}, \eqn{nA_i} is a scalar depending on the cluster sizes,
  and \eqn{rho} is the ANOVA estimate of the intra-cluster correlation, assumed common
  across groups (see Donner, 1989 or Donner et al., 1994). The statistic is compared to a \eqn{\chi^2}{chi-squared}
  distribution with \eqn{I - 1} degrees of freedom. Fixed correction factors can be specified
  with the argument \code{C}.
  }
\value{An object of formal class \dQuote{drs}: see \code{\link{drs-class}} for details. The slot \code{tab}
  provides the proportion of successes and the correction factor for each group.
  }
\references{
  Donner, A., 1989. \emph{Statistical methods in ophthalmology: an adjusted chi-squared approach}. Biometrics 45, 605-611.\cr
  Donner, A., 1993. \emph{The comparison of proportions in the presence of litter effects}. Prev. Vet. Med. 18, 17-26.\cr
  Donner, A., Eliasziw, M., Klar, N., 1994. \emph{A comparison of methods for testing homogeneity of proportions in 
  teratologic studies}. Stat. Med. 13, 1253-1264.
  }
\author{Matthieu Lesnoff \email{matthieu.lesnoff@cirad.fr}, Renaud Lancelot \email{renaud.lancelot@cirad.fr}}
\seealso{\code{\link{chisq.test}}, \code{\link{raoscott}}, \code{\link{drs-class}} }
\examples{
  data(rats)
  donner(formula = cbind(y, n - y) ~ group, data = rats)
  donner(formula = y/n ~ group, weights = n, data = rats)
  donner(response = cbind(y, n - y), group = group, data = rats)
  donner(response = y/n, weights = n, group = group, data = rats)
  # standard test
  donner(cbind(y, n - y) ~ group, data = rats, C = c(1, 1))
  data(antibio)
  donner(cbind(y, n - y) ~ treatment, data = antibio)
  }
\keyword{htest}
