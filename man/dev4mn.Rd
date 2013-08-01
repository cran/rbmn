\name{dev4mn}
\alias{dev4mn}
\title{
Computes the deviance for a sample of multinormal vector
}
\description{
From the \code{n} observed values of a vector of size \code{p} (Y),
their expectations (EY) and the variance matrix (VY) supposed
identical for all vectors, returns the deviance, i.e.
\code{-2*log(p(Y))}.
}
\usage{
dev4mn(Y, EY, VY)
}
\examples{
 dev4mn(matrix(runif(3), 1), t(rbmn0mn.01$mu), rbmn0mn.01$gamma);
}
\arguments{
\item{Y}{Matrix \code{nxp} of the \code{n} observed values of length
\code{p}.}
\item{EY}{Expectation of \code{Y} (matrix \code{nxp} or vector
\code{p}).}
\item{VY}{Matrix of the variance of each row of \code{Y} (matrix
\code{pxp}).}
}
\details{
When \code{EY} is a vector with length \code{ncol(Y)} this supposes
that all observations have the same expectation.
}
\value{
A scalar
}
