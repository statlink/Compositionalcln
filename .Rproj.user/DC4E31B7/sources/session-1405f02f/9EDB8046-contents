\name{cln.residplot}
\alias{cln.residplot}
\title{
Residual diagnostic plot for the CLN and CLNR models
}

\description{
Residual diagnostic plot for the CLN and CLNR models.
}

\usage{
cln.residplot(y, m = NULL, S, x = NULL, beta = NULL)
}

\arguments{
\item{y}{
A numerical matrix with the compositional data, possibly containing zero values. The last
column is the reference component of the additive log-ratio transformation.
}
\item{m}{
The estimated mean vector of the alr transformed data, of length \eqn{D-1}{D-1}, as returned
by \code{\link{cln.mle}}. It is ignored when \code{x} is supplied.
}
\item{S}{
The estimated \eqn{(D-1) \times (D-1)}{(D-1) x (D-1)} covariance matrix of the (residuals of
the) alr transformed data, as returned by \code{\link{cln.mle}} or \code{\link{cln.reg}}.
}
\item{x}{
The predictor variable(s), the same as those given to \code{\link{cln.reg}}. Do not supply a
design matrix, since the intercept is added internally. If \code{NULL} (default), the model
without covariates is assumed.
}
\item{beta}{
The \eqn{p \times (D-1)}{p x (D-1)} matrix of the estimated regression coefficients, as
returned by \code{\link{cln.reg}}, where the first row contains the intercepts. Required if
\code{x} is supplied.
}
}

\details{
For the \eqn{i}{i}-th composition let \eqn{b_i}{b_i} be the vector of the log-ratios among
its \eqn{C_i}{C_i} nonzero components, each taken against the last nonzero one. Under the
model, \eqn{b_i}{b_i} follows a multivariate normal distribution with mean
\eqn{Q_i \mu_i}{Q_i mu_i} and covariance matrix \eqn{Q_i \Sigma Q_i^\top}{Q_i Sigma Q_i'}, where
\eqn{Q_i}{Q_i} depends only on the zero pattern and \eqn{\mu_i}{mu_i} is either the common mean
or \eqn{B^\top x_i}{B'x_i}. The squared Mahalanobis distance
\deqn{m_i = (b_i - Q_i \mu_i)^\top (Q_i \Sigma Q_i^\top)^{-1} (b_i - Q_i \mu_i)}{
m_i = (b_i - Q_i mu_i)' (Q_i Sigma Q_i')^{-1} (b_i - Q_i mu_i)}
is compared to a chi-square distribution with \eqn{c_i = C_i - 1}{c_i = C_i - 1} degrees of
freedom. The chi-square probabilities are transformed to normal scores,
\eqn{z_i = \Phi^{-1}(G_{c_i}(m_i))}{z_i = qnorm(G_{c_i}(m_i))}, which have a common scale for
all rows, and are plotted against the standard normal quantiles. Only the observed part of each
vector is used, so no zero value is replaced.

The reference distribution ignores that the parameters were estimated from the same data, so
for small samples or many predictors the plot may look slightly better than the fit deserves.
Compositions with only one nonzero component contain no log-ratio and are omitted (\code{NA}).
}

\value{
The QQ-plot is produced and a list including:
\item{m2}{
The squared Mahalanobis distance \eqn{m_i}{m_i} of the observed log-ratios of each
row. \code{NA} for rows with a single nonzero component.
}
\item{dof}{
The degrees of freedom of each distance, that is, the number of nonzero components
minus one.
}
\item{cdf}{
The value of the chi-square cumulative distribution function with \code{dof} degrees
of freedom evaluated at \code{m2}. Under a correct model it is approximately uniform. It is not a
p-value, whose upper tail is \code{1 - cdf}.
}
\item{z}{
The normal scores of \code{cdf}, that is, the quantities plotted. Under a correct model
they are approximately standard normal, and large positive values indicate a poorly fitted
composition.
}
}

\references{
Alzeley O. and Tsagris M. (2026). Modelling compositional data with structural zero values
https://arxiv.org/pdf/2608.29954

Aitchison J. (2003). The statistical analysis of compositional data. New Jersey: Reprinted by
The Blackburn Press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{\link{cln.mle}}, \code{\link{cln.reg}}
}

\examples{
x <- as.vector(iris[, 4])
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
ind <- sample(150, 10)
for ( k in ind )  y[k, sample(3, 1)] <- 0
y <- y / rowSums(y)
## without covariates
mod <- cln.mle(y)
cln.residplot(y, m = mod$m, S = mod$S)

## with covariates
mod <- cln.reg(y, x)
cln.residplot(y, S = mod$S, x = x, beta = mod$beta)
}
