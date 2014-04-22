\name{methods}
\alias{print.mcrot}
\alias{print.mcroastci}
\alias{print.mixedrot}
\alias{print.glmrot}

\title{mcroast Methods}
\description{
Print method for the output of function mcroast
}
\usage{
\method{print}{mcrot}(x, digits = max(3, getOption("digits") - 3), ...)
\method{print}{mixedrot}(x, digits = max(3, getOption("digits") - 3), ...)
\method{print}{glmrot}(x, digits = max(3, getOption("digits") - 3), ...)
\method{print}{mcroastci}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{an object of class mcrot, mixedrot, or glmrot}
  \item{digits}{minimum number of significant digits}
  \item{\dots}{further arguments passed to or from other methods.}
}


\author{
Daniel Gerhard
}


\keyword{methods}
