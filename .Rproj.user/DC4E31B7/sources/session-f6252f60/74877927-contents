\name{boot.clnmle}
\alias{boot.clnmle}

\title{
Maximum likelihood estimation of the conditional logistic normal model
}

\description{
Maximum likelihood estimation of the conditional logistic normal model.
}

\usage{
boot.clnmle(y, tol = 1e-6, maxit = 500, R = 1000)
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
\item{R}{
The number of bootstrap resamples to perform.
}
}

\details{
The function performs non-parametric bootstrap for the conditional multivariate
logistic normal distribution (Alzeley and Tsagris, 2026).
}

\value{
A list including:
\item{m}{
The bootsrap mean vectors in the Euclidean space.
}
\item{mesi}{
The bootstrap mean vectors in the simplex space.
}
\item{S}{
The bootstrap estimated covariance matrices.
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
mod <- boot.clnmle(y, R = 300)
}


