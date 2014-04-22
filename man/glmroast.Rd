\name{glmroast}
\alias{glmroast}
\title{
Multiple Contrast Rotation Tests based on Marginal Generalized Linear Models
}
\description{
Adaption of the roast function in package limma to perform a maxT-test for parameters in multiple generalized linear models, controlling FWER for multiple contrasts.
}
\usage{
glmroast(formula, data = NULL, K, family=gaussian(), nrot = 9999, adjusted = TRUE, responsenumber=NULL, ...)
}

\arguments{
  \item{formula}{A formula object, analogously to lmer, with a matrix of response variables (independent subjects in rows) on the left, and design covariates + random effects on the right hand side.}
  \item{data}{A data.frame object, searching for column names in formula.}
  \item{K}{A contrast matrix (different contrasts in rows), with the same number of columns as number of parameters in the model}
  \item{family}{An object of class family, or a list with family objects for each marginal model.}  
  \item{nrot}{Number of Monte Carlo simulation steps}
  \item{adjusted}{logical, if TRUE single-step maxT FWER control is performed; if FALSE the comparison-wise error rate is controlled}
  \item{responsenumber}{An integer vector, as a unique identifier of the response for each marginal model. Default is NULL, being equal to a sequence of integers from 1 to the number of response matrix columns.}
  \item{...}{Additional arguments passed to the glm function.}
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
