\name{cln.mle}
\alias{cln.mle}
\alias{cln.mle2}

\title{
Maximum likelihood estimation of the conditional logistic normal model
}

\description{
Maximum likelihood estimation of the conditional logistic normal model.
}

\usage{
cln.mle(y, tol = 1e-6, maxit = 500)
cln.mle2(y, tol = 1e-6, maxit = 500)
}

\arguments{
\item{y}{
A numerical matrix with compositional data with zero values.
}
\item{tol}{
The tolerance value to terminate the EM algorithm.
}
\item{maxit}{
The maximum number of iterations allowed for the EM algorithm.
}
}

\details{
The function fits the conditional multivariate logistic normal distribution (Alzeley and Tsagris, 2026).

The function \code{cln.mle2} is a C++ implementation of \code{cln.mle} and returns the same output.
The observations with zero values are grouped by their pattern of zeros, so the conditional
moments of the E-step are computed once per pattern, which makes it considerably faster.
}

\value{
A list including:
\item{m}{
The mean vector in the Euclidean space.
}
\item{mesi}{
The mean vector in the simplex space.
}
\item{S}{
The estimated covariance matrix.
}
\item{loglik}{
The log-likelhiood value.
}
\item{patterns}{
A matrix with the patterns of zeros, where the value of 0 indicates the presence of a zero,
and the last column contains the percentage of occurrence each pattern. This is useful for
random values simulation.
}
\item{iters}{
The number of iterations required by the EM algorithm.
}
}

\references{
Alzeley O. and Tsagris M. (2026). Modelling compositional data with structural zero values
https://arxiv.org/pdf/2608.29954

Aitchison J. (1986). The statistical analysis of compositional data, pages 272--274.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{cln.reg}
}
}

\examples{
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
ind <- sample(150, 10)
for ( k in ind )  y[k, sample(3, 1)] <- 0
y <- y / rowSums(y)
cln.mle(y)
cln.mle2(y)
}


