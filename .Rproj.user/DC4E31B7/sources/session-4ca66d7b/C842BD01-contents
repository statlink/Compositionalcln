\name{cln.cr}
\alias{cln.cr}
\title{
Confidence region of the simplicial mean on a ternary diagram
}

\description{
Confidence region of the simplicial mean on a ternary diagram.
}

\usage{
cln.cr(y, m, S, alpha = 0.05, dg = TRUE, hg = TRUE, ploty = TRUE)
}

\arguments{
\item{y}{
A matrix with the compositional data (with zero values).
}
\item{m}{
A value with the mean vector in \eqn{S^2}.
}
\item{S}{
The covariance matrix in \eqn{S^2}.
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
The confidence region of the mean in \eqn{S^2} is constructed and then mapped back onto the simplex.
}

\value{
A ternary diagram with the simplicial mean and the confidence region.
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

\seealso{
\code{ \link{cln.contour}, \link{cln.mle}
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
cln.cr(y, m, S, ploty = FALSE)
}

