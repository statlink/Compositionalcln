\name{cln.regs}
\alias{cln.regs}
\title{
Many conditional logistic normal regressions
}

\description{
Many conditional logistic normal regressions.
}

\usage{
cln.regs(y, X, tol = 1e-6, maxit = 500)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). The number of observations
(vectors) with no zero values should be more than the columns of the predictor variables.
Otherwise, the initial values will not be calculated.
}
\item{X}{
A matrix with the predictor variable(s). A regression model will be fit to each
of them, separately.
}
\item{tol}{
The tolerance value to terminate the EM algorithm.
}
\item{maxit}{
The maximum number of iterations allowed for the EM algorithm.
}
}

\details{
A conditional logistic normal regression is being fitted for each of the columns
of X.
}

\value{
A list including:
\item{stat}{
A vector with the log-likelihood ratio test statistics.
}
\item{pvalue}{
A vector with the corresponding logged p-values.
}
}

\references{
Alzeley O. and Tsagris M. (2026). Modelling compositional data with structural zero values
https://arxiv.org/pdf/2608.29954
}


\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{cln.reg}
}
}

\examples{
x <- matrix( rnorm(150 * 10), ncol = 10 )
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
ind <- sample(150, 10)
for ( k in ind )  y[k, sample(3, 1)] <- 0
y <- y / rowSums(y)
mod <- cln.regs(y, x)
}

