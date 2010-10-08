\name{raoscott}
\alias{raoscott}
\alias{show,raoscott-class}
\title{Test of Proportion Homogeneity using Rao and Scott's Adjustment}
\description{
  Tests the homogeneity of proportions between \eqn{I} groups (H0: \eqn{p_1 = p_2 = ... = p_I} ) from clustered binomial 
  data \eqn{(n, y)} using the adjusted \eqn{\chi^2}{chi-squared} statistic proposed by Rao and Scott (1993).
  }
\usage{raoscott(formula = NULL, response = NULL, weights = NULL, 
              group = NULL, data, pooled = FALSE, deff = NULL)}
\arguments{
  \item{formula}{An optional formula where the left-hand side is either a matrix of the form \code{cbind(y, n-y)},
    where the modelled probability is \code{y/n}, or a vector of proportions to be modelled (\code{y/n}). 
    In both cases, the right-hand side must specify a single grouping variable. When the left-hand side of the formula 
    is a vector of proportions, the argument \code{weight} must be used to indicate the denominators of the 
    proportions.}
  \item{response}{An optional argument: either a matrix of the form \code{cbind(y, n-y)}, where the modelled probability 
    is \code{y/n}, or a vector of proportions to be modelled (\code{y/n}).}
  \item{weights}{An optional argument used when the left-hand side of \code{formula} or \code{response} is a vector
    of proportions: \code{weight} is the denominator of the proportions.}
  \item{group}{An optional argument only used when \code{response} is used. In this case, this argument is a factor 
    indicating a grouping variable.}
  \item{data}{A data frame containing the response (\code{n} and \code{y}) and the grouping variable.}
  \item{pooled}{Logical indicating if a pooled design effect is estimated over the \eqn{I} groups.}
  \item{deff}{A numerical vector of \eqn{I} design effects.}
  }
\details{
  The method is based on the concepts of design effect and effective sample size.\cr\cr
  The design effect in each group \eqn{i} is estimated by \eqn{deff_i = vratio_i / vbin_i}, where \eqn{vratio_i} is 
  the variance of the ratio estimate of the probability in group \eqn{i} (Cochran, 1999, p. 32 and p. 66)
  and \eqn{vbin_i} is the standard binomial variance. A pooled design effect (i.e., over the \eqn{I} groups)
  is estimated if argument \code{pooled = TRUE} (see Rao and Scott, 1993, Eq. 6). Fixed design effects can be specified
  with the argument \code{deff}.\cr
  The \eqn{deff_i} are used  to compute the effective sample sizes \eqn{nadj_i = n_i / deff_i}, the effective numbers
  of successes \eqn{yadj_i = y_i / deff_i} in each group \eqn{i}, and the overall effective proportion 
  \eqn{padj = \sum_{i} yadj_i / \sum_{i} deff_i}{padj = sum(yadj_i) / sum(deff_i)}.
  The test statistic is obtained by substituting these quantities in the usual \eqn{\chi^2}{chi-squared} statistic,
  yielding:
  \deqn{X^2 = \sum_{i}\frac{(yadj_i - nadj_i * padj)^2}{nadj_i * padj * (1 - padj)}}{%
    X^2 = sum( (yadj_i - nadj_i * padj)^2 / (nadj_i * padj * (1 - padj)) )}
  which is compared to a \eqn{\chi^2}{chi-squared} distribution with \eqn{I - 1} degrees of freedom.
  }
\value{An object of formal class \dQuote{drs}: see \code{\link{drs-class}} for details. The slot \code{tab}
  provides the proportion of successes, the variances of the proportion and the design effect for each group.
  }
\references{
  Cochran, W.G., 1999, 2nd ed. \emph{Sampling techniques}. John Wiley & Sons, New York.\cr
  Rao, J.N.K., Scott, A.J., 1992. \emph{A simple method for the analysis of clustered binary data}.
  Biometrics 48, 577-585.
  }
\author{Matthieu Lesnoff \email{matthieu.lesnoff@cirad.fr}, Renaud Lancelot \email{renaud.lancelot@cirad.fr}}
\seealso{ \code{\link{chisq.test}}, \code{\link{donner}}, \code{\link{iccbin}}, \code{\link{drs-class}} }
\examples{
  data(rats)
  # deff by group
  raoscott(cbind(y, n - y) ~ group, data = rats)
  raoscott(y/n ~ group, weights = n, data = rats)
  raoscott(response = cbind(y, n - y), group = group, data = rats)
  raoscott(response = y/n, weights = n, group = group, data = rats)
  # pooled deff
  raoscott(cbind(y, n - y) ~ group, data = rats, pooled = TRUE)
  # standard test
  raoscott(cbind(y, n - y) ~ group, data = rats, deff = c(1, 1))
  data(antibio)
  raoscott(cbind(y, n - y) ~ treatment, data = antibio)
  }
\keyword{htest}
