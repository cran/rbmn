\name{simulate8mn}
\alias{simulate8mn}
\title{
simulates a multinormal vector
}
\description{
returns a matrix of simulated values with the variable in columns and
the simulations in rows.
}
\usage{
simulate8mn(mn, nbs, tol=1e-7)
}
\examples{
 print(simulate8mn(rbmn0mn.01, 12));
}
\arguments{
\item{mn}{ list defining the distribution of the initial vector with
\code{$mu}, its expectation, and \code{$gamma}, its variance matrix.}
\item{nbs}{ number of simulations to return.}
\item{tol}{ tolerance value to be transmitted to \code{mvrnorm}.}
}
\details{
Just a call to the basic function \code{mvrnorm}. Names of the
variables are taken from those of \code{mn$mu}, when these does not
exist, standard ones are provided.
}
\value{
A matrix/data frame of size : \code{nbs x length(mn$mu)}
}
