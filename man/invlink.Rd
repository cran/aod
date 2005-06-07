\name{invlink}
\alias{invlink}
\title{Transformation from the Link Scale to the Observation Scale}
\description{
  The function transforms a variable from the link scale to the observation scale: probability or count.
  }
\usage{
  invlink(x, type = c("cloglog", "log", "logit"))
  }
\arguments{
  \item{x}{A vector of real numbers.}
  \item{type}{A character string. Legal values are \dQuote{cloglog}, \dQuote{log} and \dQuote{logit}.}
  }
\value{
  \eqn{anti-logit(x) = exp(x) / (1 + exp(x))}\cr
  \eqn{anti-cloglog(x) = 1 - exp(-exp(x))}
  }
\seealso{\code{\link{link}}}
\examples{
  x <- seq(-5, 5, length = 100)
  plot(x, invlink(x, type = "logit"),
       type = "l", lwd = 2, ylab = "Probability")
  lines(x, invlink(x, type = "cloglog"), lty = 2, lwd = 2)
  grid(col = "black")
  legend(-5, 1, legend = c("alogit(x)", "acloglog(x)"),
         lty = c(1, 2), bg = "white")
  }
\keyword{math}
