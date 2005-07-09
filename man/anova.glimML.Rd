\name{anova-methods}
\docType{methods}
\alias{anova-methods}
\alias{anova,ANY-method}
\alias{anova,glimML-method}
\alias{anova.glimML-class}
\alias{show,anova.glimML-method}

\title{Likelihood-Ratio Tests for Nested ML Models}

\description{
  Performs likelihood-ratio tests on nested models. Currently, one method was implemented
  for beta-binomial models (\code{betabin}) or negative-binomial models (\code{negbin}).}

\usage{
  \S4method{anova}{glimML}(object, \dots)
  }

\arguments{
  \item{object}{Fitted model of class \dQuote{glimML}.}
  \item{\dots}{Further models to be tested or arguments passed to the \code{print} function.}
  }
\details{The \code{anova} method for models of formal class \dQuote{glimML} needs at least 2 nested models of the 
  same type (either beta-binomial or negative-binomial models: they cannot be mixed). The quantity of interest is 
  the deviance difference between the compared models: it is a log-likelihood ratio statistic. Under the null 
  hypothesis that 2 nested models fit the data equally well, the deviance difference has an approximate 
  \eqn{\chi^2}{chi-squared} distribution with degrees of freedom = the difference in the number of parameters between 
  the compared models (Mc Cullagh and Nelder, 1989).}
\value{
  An object of formal class \dQuote{anova.glimML} with 3 slots:
  \item{models}{A vector of character strings with each component giving the name of the models and the formulas for 
    the fixed and random effects.}
  \item{anova.table}{A data frame containing the results. Row names correspond to the models.
    \tabular{lll}{
      logL        \tab numeric   \tab maximized log-likelihood\cr
      k           \tab numeric   \tab number of parameters in the model\cr
      AIC         \tab numeric   \tab Akaike information criterion for the model\cr
      BIC         \tab numeric   \tab Bayesian information criterion the model\cr
      Resid. dev. \tab numeric   \tab Residual deviance\cr
      Resid. Df   \tab numeric   \tab df of the residuals\cr
      Test        \tab character \tab Nested models which are tested\cr
      Deviance    \tab numeric   \tab Deviance difference between the 2 models\cr
      Df          \tab numeric   \tab df associated with deviance difference\cr
      P(> Chi2)   \tab numeric   \tab P value associated with H0.
      }
    }
  \item{type}{A character chain indicating the kind of fitted model: \dQuote{BB} for beta-binomial, or \dQuote{NB} for 
     negative-binomial model.}
  }
\section{Warning}{The comparison between 2 or more models will only be valid if they are fitted to the same data set.}

\references{McCullagh, P., Nelder, J.A., 1989. \emph{Generalized linear models}. London, Chapman & Hall, 511 p.\cr 
  See Appendix C. Likelihood ratio statistics, p. 476-478.}

\seealso{\code{\link{anova.glm}}, \code{\link[aod:AIC]{AIC}}}

\examples{
  data(orob2)
  # likelihood ratio test for the effect of root
  fm1 <- betabin(cbind(y, n - y) ~ seed, ~ 1, data = orob2)
  fm2 <- betabin(cbind(y, n - y) ~ seed + root, ~ 1, data = orob2)
  anova(fm1, fm2)
  }

\keyword{regression}
