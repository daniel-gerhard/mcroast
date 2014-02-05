\name{mcroast}
\alias{mcroast}
\title{
Multiple Contrast Rotation Tests
}
\description{
Adaption of the roast function in package limma to perform a maxT-test in the multivariate linear model, controlling FWER for multiple contrasts.
}
\usage{
mcroast(formula, data = NULL, K, nrot = 9999, adjusted = TRUE, moderated = FALSE)
}

\arguments{
  \item{formula}{A formula object, analogously to lm, with a matrix of response variables (independent subjects in rows) on the left, and design covariates on the right hand side.}
  \item{data}{A data.frame object, searching for column names in formula.}
  \item{K}{A contrast matrix (different contrasts in rows), with the same number of columns as number of parameters in the model}
  \item{nrot}{Number of Monte Carlo simulation steps}
  \item{adjusted}{logical, if TRUE single-step maxT FWER control is performed; if FALSE the comparison-wise error rate is controlled}
  \item{moderated}{logical, if TRUE, moderated statistics (Smyth 2004) are used.}
}
\details{
Most of the code was taken from function roast in package limma (Bioconductor).
}
\value{
A list with
\item{Y}{The matrix of response variables}
\item{X}{The design matrix}
\item{statistic}{A matrix of test statistics}
\item{pvalue}{A matrix of (adjusted) p-values}
\item{B}{A matrix with contrast parameter estimates}
\item{SD}{A matrix with estimated standard deviations}
\item{mrot}{A vector with simulated maxT statistics under H0}
}
\references{
Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader, JE, and Smyth, GK (2010). ROAST: rotation gene set tests for complex microarray experiments. Bioinformatics 26, 2176-2182. 

Langsrud, O (2005). Rotation tests. Statistics and Computing 15, 53-60.
}
\author{
Daniel Gerhard
}


\keyword{multivariate}
\keyword{htest}
