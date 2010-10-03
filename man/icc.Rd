\name{icc}
\alias{icc}
\alias{show,icc-method}
\title{Intra-Cluster Correlation for Binomial Data}

\description{This function computes approximate estimates of the intra-cluster correlation \eqn{\rho} 
from clustered, binomial data \eqn{{(n_1, y_1), (n_2, y_2), ..., (n_K, y_K)}} (with \eqn{K} the number of clusters),
using a 1-way linear random effect model for continuous variables. Both ML/REML and ANOVA estimates are provided.
}

\usage{icc(n, y, data, method = c("REML", "ML"), R = NULL)}
\arguments{
  \item{n}{Vector of the denominators of the proportions.}
  \item{y}{Vector of the numerators of the proportions.}
  \item{data}{A data frame containing the variables n and y.}
  \item{method}{A character string \dQuote{ML} (maximum likelihood) or \dQuote{REML} (restricted ML) used in the
    ML/REML estimation method for \eqn{\rho}. Default to \dQuote{REML}.}
  \item{R}{A scalar indicating the number of Monte Carlo (MC) replicates used to estimate the confidence interval of
    \eqn{\rho} for ML/REML estimates. When R is NULL (default), the MC confidence interval is not computed.}
  }

\details{
  The function assumes an homogeneous proportion \eqn{p} across the clusters.
  Clustered data are split into binary (0/1) observations \eqn{y_{ij}} (obs. \eqn{j} in cluster
  \eqn{i}). A linear 1-way random effect model representing these binary data is defined as:
  \deqn{y_{ij} = \mu + u_i + e_{ij}}
  where \eqn{\mu} is the general mean, \eqn{a_i} is the random cluster effect and \eqn{e_{ij}} is the residual error.
  Model assumptions are: \eqn{u \sim N(0, \sigma_u^2)}{u ~ N(0, \sigma_u^2)} and 
  \eqn{e_{ij} \sim N(0, \sigma_e^2)}{e_{ij} ~ N(0, \sigma_e^2)}.
  The intra-cluster correlation is computed as \deqn{\rho = \sigma_u^2 / (\sigma_u^2 + \sigma_e^2)}.
  Variances components \eqn{\sigma_u^2} and \eqn{\sigma_e^2} are estimated using ML/REML and (moments) ANOVA methods.\cr\cr
  
  For the ML/REML method, variance components (actually, vector \eqn{\nu = log(\sigma_u, \sigma_e)}) are 
  estimated with the function \code{lme} (package \pkg{nlme}): see Pinheiro and Bates, 2000. 
  An approximate variance of \eqn{\rho} is computed with the \emph{Delta} method from
  the estimated \eqn{Var[\nu]} (= matrix \code{apVar} provided in the \code{lme} output).
  If the argument \code{R} is not null, a MC confidence interval of \eqn{\rho} is computed assuming
  that \eqn{\nu \sim N(\nu, Var[\nu])}{\nu ~ N(\nu, Var[\nu])}.\cr\cr

  For the ANOVA method, moments estimates of the variance components are used.
  See Donner (1989), Searle et al. (1992) or Zou and Donner (2004).
  A \eqn{F} test is provided to assess whether \eqn{\rho = 0} (actually, whether
  \eqn{\sigma_u^2 = 0}: see Searle et al, 1992, p. 76).
  }

\value{
  An object of formal class \dQuote{icc}, with 7 slots:
  \item{CALL}{The call of the function.}
  \item{rho}{A numeric vector with 2 components: REML (or ML) and ANOVA, containing the estimated values of \eqn{\rho}
    according to these 2 methods.}
  \item{varrho}{A numeric scalar giving the approximate variance of \eqn{\rho} (ML/REML estimate)
   estimated with the \emph{Delta} method.}
  \item{f}{A numeric vector with the results of the \eqn{F} test.}
  \item{rho.MC}{A numeric vector with the MC replicates of \eqn{\rho}.}
  \item{method}{A character string taking values \dQuote{ML} or \dQuote{REML}.}
  \item{features}{A numeric vector with 3 components summarizing the main features of the data: \code{N} = number 
    of clusters, \code{n} = number of subjects, \code{y} = number of cases.}
  }

\references{
  Donner, A., 1989. \emph{Statistical methods in ophthalmology: an adjusted chi-squared approach}.
  Biometrics 45, 605-611.\cr
  Pinheiro, J.C., Bates, D.M., 2000. \emph{Mixed-effects models in S and S-PLUS}. Springer-Verlag, New York.\cr
  Searle, S.R., Casella, G., McCulloch, C.E., 1992. \emph{Variance components}. Wiley, New York.\cr
  Zou, G., Donner, A., 2004. \emph{Confidence interval estimation of the intraclass correlation coefficient for 
  binary outcome data}. Biometrics 60, 807-811.
  }

\author{Matthieu Lesnoff \email{matthieu.lesnoff@cirad.fr}, Renaud Lancelot \email{renaud.lancelot@cirad.fr}}

\seealso{\code{\link{icc-class}}, \code{\link[nlme:lme]{lme}}\cr
         \code{icc} in the contributed packages \pkg{irr} and \pkg{psy}.}

\examples{
  data(rats)
  icc(n, y, rats[rats$group == "CTRL", ])
  \dontrun{
    res <- icc(n, y, rats[rats$group == "TREAT", ], R = 5000)
    res
    hist(res@rho.MC)
    }
  by(rats,
     list(group = rats$group),
     function(x) icc(n, y, data = x))
  }

\keyword{htest}
