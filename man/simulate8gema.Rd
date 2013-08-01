\name{simulate8gema}
\alias{simulate8gema}
\title{
simulates from a /gema/ object
}
\description{
returns a matrix of simulated values with the variable in columns and
the simulations in rows.
}
\usage{
simulate8gema(gema, nbs)
}
\examples{
 simulate8gema(rbmn0gema.01, 10);
}
\arguments{
\item{gema}{ The \code{gema} object.}
\item{nbs}{ number of simulations to return.}
}
\details{
Just the application of the standard formula to a white noise.
Variables names are taken from those of \code{gema$mu}, when these
does not exist, standard ones are provided.
}
\value{
A matrix of size : \code{nbs x length(gema$mu)}
}
