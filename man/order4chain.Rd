\name{order4chain}
\alias{order4chain}
\title{
returns a topological order of a /chain/ or checks a proposed order.
}
\description{
From a \code{chain} object returns one of the possible topological
orders, through a permutation when \code{is.null(ord)}. If not
\code{ord} must be a proposed order to be checked given as a
permutation if \code{is.numeric(ord)} or a vector of ordered names if
\code{is.character(ord)}.
}
\usage{
order4chain(chain, ord=NULL)
}
\examples{
 order4chain(rbmn0chain.02);
 order4chain(rbmn0chain.02, order4chain(rbmn0chain.02));
}
\arguments{
\item{chain}{ the \code{chain} object to be considered.}
\item{ord}{ Indicates what must be done. \code{NULL} to get a
topological order associated to the chain otherwise a permutation to
be checked as one of the possible topological orders of the chain.}
}
\details{
For the moment the \code{ord} option is bad and an error message is
returned when used.
}
\value{
a permutation vector of the nodes of the /nbn/ or a named character
with the nodes not having their parents before them; when it is of
length zero this means that the check was successful.
}
