\name{cln.bcr}
\alias{cln.bcr}
\title{
Non-parametric confidence region of the simplicial mean on a ternary diagram
}

\description{
Non-parametric confidence region of the simplicial mean on a ternary diagram.
}

\usage{
cln.bcr(y, tol = 1e-6, maxit = 500, alpha = 0.05, R = 1000,
dg = TRUE, hg = TRUE, ploty = TRUE)
}

\arguments{
\item{y}{
A matrix with the compositional data (with zero values).
}
\item{tol}{
The tolerance value to terminate the EM algorithm.
}
\item{maxit}{
The maximum number of iterations allowed for the EM algorithm.
}
\item{alpha}{
The significance level to construct the (1-alpha)\% confidence region.
}
\item{dg}{
Do you want diagonal grid lines to appear? If yes, set this TRUE.
}
\item{hg}{
Do you want horizontal grid lines to appear? If yes, set this TRUE.
}
\item{ploty}{
Do you want to plot the data as well?
}
}

\details{
A non-parametric confidence region of the mean in \eqn{S^2} is constructed based
on the method of Fiksel, Zeger and Datta (2022).
}

\value{
A ternary diagram with the simplicial mean and the non-parametric confidence region.
}

\references{
Fiksel J., Zeger S. and Datta A. (2022). A transformation-free linear regression for
compositional outcomes and predictors. Biometrics, 78(3): 974--987.

Alzeley O. and Tsagris M. (2026). Modelling compositional data with structural zero values
https://arxiv.org/pdf/2608.29954

Aitchison J. (1986). The statistical analysis of compositional data, pages 272--274.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{ \link{cln.cr}, \link{cln.contour}, \link{cln.mle}
}
}

\examples{
y <- as.matrix(iris[1:50, 1:3])
y <- y / rowSums(y)
ind <- sample(50, 10)
for ( k in ind )  y[k, sample(3, 1)] <- 0
y <- y / rowSums(y)
mod <- cln.mle(y)
m <- mod$m
S <- mod$S
cln.bcr(y, m, S, R= 300, ploty = FALSE)
}

