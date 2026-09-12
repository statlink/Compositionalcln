\name{Ternary diagram}
\alias{ternary}

\title{
Ternary diagram
}

\description{
Ternary diagram.
}

\usage{
ternary(y, dg = FALSE, hg = FALSE, m = NULL, colour = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{dg}{
Do you want diagonal grid lines to appear? If yes, set this TRUE.
}
\item{hg}{
Do you want horizontal grid lines to appear? If yes, set this TRUE.
}
\item{m}{
If you know the estimated mean pass it here, otherwise leave it NULL. It will appear
with the symbol +.
}
\item{colour}{
If you want the points to appear in different colour put a vector with the colour
numbers or colours.
}
}

\details{
There are two ways to create a ternary graph. We used here that one where each edge is equal to 1
and it is what Aitchison (1986) uses. For every given point, the sum of the distances from the edges
is equal to 1. Horizontal and or diagonal grid lines, and the mean can appear.
Zeros in the data appear with red x.
}

\value{
The ternary plot. Additionally, horizontal or diagonal grid lines can appear as well.
}

\references{
Aitchison, J. (1983). Principal component analysis of compositional data. Biometrika 70(1): 57--65.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{cln.contour}
}
}

\examples{
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
ind <- sample(150, 10)
for ( k in ind )  y[k, sample(3, 1)] <- 0
y <- y / rowSums(y)
mod <- cln.mle(y)
ternary(y, hg = TRUE, dg = TRUE, m = mod$mesi)
}

