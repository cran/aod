\name{summary,aic-method}
\docType{methods}
\alias{summary,aic-method}
\alias{show,aic-method}

\title{Akaike Information Statistics}

\description{Computes Akaike difference and Akaike weights from an object of formal class \dQuote{aic}.}

\usage{
  \S4method{summary}{aic}(object, which = c("AIC", "AICc"))
  }

\arguments{
  \item{object}{An object of formal class \dQuote{aic}.}
  \item{which}{A character string indicating which information criterion is selected to compute Akaike difference and
    Akaike weights: either \dQuote{AIC} or \dQuote{AICc}.}
  }

\section{Methods}{\describe{
  \item{summary}{The models are ordered according to AIC or AICc and 3 statistics are computed:
  
      - the \emph{Akaike difference} \eqn{\Delta}: the change in AIC (or AICc) between successive (ordered) models,
      
      - the \emph{Akaike weight} \eqn{W}: when \eqn{r} models are compared, 
      \eqn{W = e^{-0.5 * \Delta} / \sum_r{e^{-\frac{1}{2} * \Delta}}}{W = exp(-0.5 * \Delta) / sum(exp(-0.5 * \Delta))},
      
      - the \emph{cumulative Akaike weight} \eqn{cum.W}: the Akaike weights sum to 1 for the \eqn{r} models which 
      are compared.
      }
    }
  }
  
\seealso{Examples in \code{\link{betabin}} and \code{\link[stats:AIC]{AIC}} in package \pkg{stats}.}

\references{
  Burnham, K.P., Anderson, D.R., 2002. \emph{Model selection and multimodel inference: a practical
  information-theoretic approach}. New-York, Springer-Verlag, 496 p.\cr
  Hurvich, C.M., Tsai, C.-L., 1995. \emph{Model selection for extended quasi-likelihood models in small samples}.
  Biometrics, 51 (3): 1077-1084.
  }

\keyword{methods}
