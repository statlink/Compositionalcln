\name{fbed.cln}
\alias{fbed.cln}
\title{
Forward-backward with early dropping variable selection
}

\description{
Forward-backward with early dropping variable selection.
}

\usage{
fbed.cln(y, X, prior = NULL, alpha = 0.05, univ = NULL, K = 0,
tol = 1e-6, maxit = 500)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable). The number of observations
(vectors) with no zero values should be more than the columns of the predictor variables.
Otherwise, the initial values will not be calculated.
}
\item{x}{
A matrix with the available predictor variables.
}
\item{prior}{
A matrix or a vector with some predictors that will be in the model.
}
\item{alpha}{
The significance level.
}
\item{univ}{
If you have the test statistics and the p-values (the output of \code{\link{cln.regs}}
or of \code{\link{cln.condregs}}) of the simple regressions pass it here.
}
\item{K}{
How many times should the FBED be repeated?
}
\item{tol}{
The tolerance value to terminate the EM algorithm.
}
\item{maxit}{
The maximum number of iterations allowed for the EM algorithm.
}
}

\details{
The function performs the FBED variable selection algorithm of
Borboudakis and Tsamardinos (2019).
}

\value{
A list including:
\item{univ}{
The result of the simple regressions, the same result as \code{\link{cln.regs}}.
This is a list that contains the log-likelihood ratio test statistics and their
corresponding logged p-values.
}
\item{res}{
A matrix with the selected variables, their test statistic and their logged p-value.
}
\item{info}{
A matrix with two columns, the number of variables selected at each step and the number of
tests performed at each step. Each step means at each value of K.
}
\item{runtime}{
The running time of the algorithm.
}
}

\references{
Alzeley O. and Tsagris M. (2026). Modelling compositional data with structural zero values
https://arxiv.org/pdf/2608.29954

Borboudakis G. and Tsamardinos I. (2019). Forward-backward selection with early dropping.
Journal of Machine Learning Research, 20(8): 1-39.
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
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
ind <- sample(150, 10)
for ( k in ind )  y[k, sample(3, 1)] <- 0
y <- y / rowSums(y)
X <- matrix(rnorm(150 * 30), ncol = 30)
mod <- fbed.cln(y, X)
}

