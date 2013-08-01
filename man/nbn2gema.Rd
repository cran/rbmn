\name{nbn2gema}
\alias{nbn2gema}
\title{
computes a /gema/ from a /nbn/
}
\description{
from a /nbn/ object defining a normal Bayesian network, computes the
vector \code{mu} and the matrix \code{li} such that if the vector
\code{E} is a vector of i.i.d. centred and standardized normal, then
\code{mu + li %
*%
 E} has the same distribution as the input /nbn/.
}
\usage{
nbn2gema(nbn)
}
\examples{
 identical(nbn2gema(rbmn0nbn.02), rbmn0gema.02);
}
\arguments{
\item{nbn}{\code{nbn} object for which the generating matrices.}
}
\value{
a list with the two following components: \code{mu} and \code{li}.
}
