\name{cln.james}
\alias{cln.james}

\title{
James test for equality of two mean vectors based on the conditional logistic normal model
}

\description{
James test for equality of two mean vectors based on the conditional logistic normal model.
}

\usage{
cln.james(y1, y2, a = 0.05, tol = 1e-6, maxit = 500)
}

\arguments{
\item{y1}{
A numerical matrix with compositional data with zero values.
}
\item{y2}{
A numerical matrix with compositional data with zero values.
}
\item{a}{
The significance level, set to 0.05 by default.
}
\item{tol}{
The tolerance value to terminate the EM algorithm.
}
\item{maxit}{
The maximum number of iterations allowed for the EM algorithm.
}
}

\details{
The function fits the conditional multivariate logistic normal distribution (Alzeley and Tsagris, 2026)
to each group, estimates the mean vectors and covariance matrices and uses them to perform the James
test for equality of two mean vectors without assuming equality of the covariance matrices.
}

\value{
A vector with 4 numbers, the value of the test statistic, the p-value, the correction
and the corrected critical value.
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
\code{ \link{cln.mle}, \link{cln.reg}
}
}

\examples{
y1 <- as.matrix(iris[1:75, 1:4])
ind <- sample(75, 10)
for ( k in ind )  y1[k, sample(3, 1)] <- 0
y1 <- y1 / rowSums(y1)

y2 <- as.matrix(iris[76:150, 1:4])
ind <- sample(75, 10)
for ( k in ind )  y2[k, sample(3, 1)] <- 0
y2 <- y2 / rowSums(y2)
cln.james(y1, y2)
}


