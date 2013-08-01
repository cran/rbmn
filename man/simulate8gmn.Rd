\name{simulate8gmn}
\alias{simulate8gmn}
\title{
simulates a multinormal vector with varying expectation
}
\description{
returns a matrix of simulated values with the variable in columns and
the simulations in rows.
}
\usage{
simulate8gmn(loi, cova, nbs, tol=1e-7)
}
\examples{
 loi <- list(mu=c(D=2, E=4), 
 rho=matrix(1:6, 2, dimnames=list(LETTERS[4:5], 
 LETTERS[1:3])), 
 gamma=matrix(c(1, 1, 1, 2), 2));
 cova <- matrix(runif(36), 12, dimnames=list(NULL, LETTERS[1:3]));
 print(simulate8gmn(loi, cova, 12));
}
\arguments{
\item{loi}{ list defining the distribution of the initial vector with
\code{$mu}, its expectation, \code{$gamma}, its variance matrix and
\code{$rho} a matrix of regression coefficients for the covariables
modifying the expectation.}
\item{cova}{ Values to give to the covariables. Must be a matrix with
\code{nbs} rows and \code{ncol(loi$rho)} columns or a vector with
\code{ncol(loi$rho)} values to be used for all simulations (i.e to
replace a matrix with identical rows..}
\item{nbs}{ number of simulations to return.}
\item{tol}{ tolerance value to be transmitted to \code{mvrnorm}.}
}
\details{
Just a call to the function \code{simulate8mn}, adding the terms to
the expectation due to the regression...
}
\value{
A matrix of size : \code{nbs x length(loi$mu)}
}
