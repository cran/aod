\name{splitbin}
\alias{splitbin}
\title{Split Grouped Data Into Individual Data}
\encoding{latin1}

\description{
  The function splits grouped data and optional covariates into individual data.
  Two types of grouped data are managed by \code{splitbin}:
     \itemize{
        \item Grouped data with weights;
        \item Grouped data of form \code{cbind(success, failure)}.
        }
  When weights, successes or failures involve non-integer numbers, these numbers are rounded before splitting.
  }

\usage{
  splitbin(formula, data, id = "idbin")
  }

\arguments{
  \item{formula}{A formula. The left-hand side represents grouped data.
    The right-hand side defines the covariates. See examples for syntax.}
  \item{data}{A data frame where all the variables described in \code{formula} are found.}
  \item{id}{An optional character string naming the identifier (= grouping factor). Default to \dQuote{idbin}.}
  }

\value{
  A data frame built according to the formula and function used in the call.
  }

\examples{
# Grouped data with weights
mydata <- data.frame(
    success = c(0, 1, 0, 1),
    f1 = c("A", "A", "B", "B"),
    f2 = c("C", "D", "C", "D"),
    n = c(4, 2, 1, 3)
    )
mydata
splitbin(formula = n ~ f1, data = mydata)$tab
splitbin(formula = n ~ f1 + f2 + success , data = mydata)$tab

# Grouped data of form "cbind(success, failure)"
mydata <- data.frame(
    success = c(4, 1),
    failure = c(1, 2),
    f1 = c("A", "B"),
    f2 = c("C", "D")
    )
mydata$n <- mydata$success + mydata$failure
mydata
splitbin(formula = cbind(success, failure) ~ 1, data = mydata)$tab
splitbin(formula = cbind(success, failure) ~ f1 + f2, data = mydata)$tab
splitbin(formula = cbind(success, n - success) ~ f1 + f2, data = mydata)$tab
splitbin(formula = cbind(success, n - 0.5 * failure - success) ~ f1 + f2,
         data = mydata)$tab
}

\keyword{datagen}

