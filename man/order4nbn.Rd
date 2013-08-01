\name{order4nbn}
\alias{order4nbn}
\title{
topological order of a /nbn/
}
\description{
returns one of the orders of the nodes such as the parents of any
node are less ranked than it when \code{is.null(ord)}. If not check
that the proposed order is either a right permutation
(\code{is.numeric(ord)}) or a vector of node names providing a
topological order (\code{is.character(ord)}).
}
\usage{
order4nbn(nbn, ord=NULL)
}
\examples{
 names(rbmn0nbn.04)[order4nbn(rbmn0nbn.04)];
}
\arguments{
\item{nbn}{\code{nbn} object for which the order must be computed.}
\item{ord}{\code{NULL} or an order to test as a permutation or a
vector of names.}
}
\details{
When \code{!is.null(ord)} the order must be an order, if not an error
is issued.
}
\value{
a permutation vector of the nodes of the /nbn/ or a named list with
the nodes not having their parents before them.
}
