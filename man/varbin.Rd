\name{varbin}
\alias{varbin}
\alias{show,varbin-class}
\title{Mean, Variance and Confidence Interval of a Proportion}
\description{
  This function computes the mean and variance of a proportion from clustered binomial data \eqn{(n, y)}, using various 
  methods. Confidence intervals are computed using a normal approximation, which might be inappropriate when the 
  proportion is close to 0 or 1.}
\usage{varbin(n, y, data, alpha = 0.05, R = 5000)}
\arguments{
  \item{n}{The denominator of the proportion.}
  \item{y}{The numerator of the proportion.}
  \item{data}{A data frame containing the data.}
  \item{alpha}{The significance level for the confidence intervals. Default to 0.05, providing 95\% CI's.}
  \item{R}{The number of bootstrap replicates to compute the bootstrap mean and variance.}
  }
\details{
  Five methods are used for the estimations. Let us
  consider \eqn{N} clusters of sizes \eqn{n_1, \ldots, n_N} with observed responses (counts) \eqn{y_1, \ldots, y_N}.
  We note \eqn{p_i = y_i / n_i} the observed proportions \eqn{(i = 1, \ldots, N)}. An underlying assumption is that the
  theoretical proportion is homogeneous across the clusters.
  
  \bold{Binomial method:} the proportion and its variance are estimated as \eqn{p = \frac{\sum_{i} y_i}{\sum_{i} n_i}}{%
  p = sum(y_i) / sum(n_i)} and \eqn{\frac{p * (1 - p)}{\sum_{i} n_i - 1}}{p * (1 - p) / sum(n_i - 1)}, respectively.
  
  \bold{Ratio method:} the one-stage cluster sampling formula is used to estimate the variance of the ratio estimate
  (see Cochran, 1999, p. 32 and p. 66). The proportion is estimated as above (\eqn{p}).
  
  \bold{Arithmetic method:} the proportion is estimated as \eqn{p_A = \frac{1}{N}\sum_{i}\frac{y_i}{n_i}}{%
   p_A = sum(y_i / n_i) / N}, with estimated variance \eqn{\frac{\sum_{i}(p_i - p_A)^2}{N * (N - 1)}}{%
   [1/(N * (N - 1))] sum((p_i - p_A)^2)}.

   \bold{Jackknife method:} the proportion \eqn{p_J} is the arithmetic mean of the pseudovalues \eqn{pv_i}, with estimated
   variance \eqn{\frac{\sum_{i}(pv_i - p_J)^2}{N * (N - 1)}}{[1/(N * (N - 1))]sum((pv_i - p_J)^2)}
   (Gladen, 1977, Paul, 1982).
   
    \bold{Bootstrap method:} \eqn{R} samples of size \eqn{N} are drawn with equal probability from the initial sample
    \eqn{(p_1, \ldots , p_N)} (Efron and Tibshirani, 1993). The bootstrap estimate \eqn{p_B} and its estimated variance 
    are the arithmetic mean and the empirical variance (computed with denominator \eqn{R - 1}) of the \eqn{R} binomial 
    estimates, respectively.
  }
\value{
  An object of formal class \dQuote{varbin}, with 5 slots:
  \item{CALL}{The call of the function.}
  \item{tab}{A 4-column data frame giving for each estimation method the mean, variance, upper and lower limits
  of the \eqn{(1 - \alpha)} confidence interval.}
  \item{boot}{A numeric vector containing the R bootstrap replicates of the proportion. Might be used to compute other
    kinds of CI's for the proportion.}
  \item{alpha}{The significance level used to compute the \eqn{(1 - \alpha)} confidence intervals.}
  \item{features}{A numeric vector with 3 components summarizing the main features of the data: \code{N} = number 
    of clusters, \code{n} = number of subjects, \code{y} = number of cases.}

  The \dQuote{show} method displays the slot \code{tab} described above, substituting the standard error to the variance.
  }
\references{
  Cochran, W.G., 1999, 3th ed. \emph{Sampling techniques}. Wiley, New York.\cr
  Efron, B., Tibshirani, R., 1993. \emph{An introduction to the bootstrap}. Chapman and Hall, London.\cr
  Gladen, B., 1977. \emph{The use of the jackknife to estimate proportions from toxicological data in the presence 
  of litter effects}. JASA 74(366), 278-283.\cr
  Paul, S.R., 1982. \emph{Analysis of proportions of affected foetuses in teratological experiments}. 
  Biometrics 38, 361-370.
 }
\seealso{ \code{\link{varbin-class}}, \code{\link[boot:boot]{boot}} }
\author{Matthieu Lesnoff \email{matthieu.lesnoff@cirad.fr}, Renaud Lancelot \email{renaud.lancelot@cirad.fr}}
\examples{
  data(rabbits)
  varbin(n, y, rabbits[rabbits$group == "M", ])
  by(rabbits,
     list(group = rabbits$group),
     function(x) varbin(n = n, y = y, data = x, R = 1000))
  }
\keyword{htest}
