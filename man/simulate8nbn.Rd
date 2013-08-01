\name{simulate8nbn}
\alias{simulate8nbn}
\title{
simulates from a /nbn/ object
}
\description{
returns a matrix of simulated values with the variable in columns and
the simulations in rows.
}
\usage{
simulate8nbn(nbn, nbs)
}
\examples{
 simulate8nbn(rbmn0nbn.01, 10);
}
\arguments{
\item{nbn}{ The \code{nbn} object.}
\item{nbs}{ number of simulations to return.}
}
\details{
Just the sequential simulations of the nodes
}
\value{
A matrix of size : \code{nbs x length(nbn)}
}
