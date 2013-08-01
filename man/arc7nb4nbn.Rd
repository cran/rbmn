\name{arc7nb4nbn}
\alias{arc7nb4nbn}
\title{
returns the number(s) of arcs of a /nbn/
}
\description{
returns the arc numbers of the node of /nbn/ object.
}
\usage{
arc7nb4nbn(nbn, each=FALSE)
}
\examples{
 arc7nb4nbn(rbmn0nbn.05);
}
\arguments{
\item{nbn}{ The \code{nbn} object to consider.}
\item{each}{ When \code{TRUE}, returns a named vector of the number
of parents of each node. If not the total number of arcs.}
}
\details{
Parents associated with a zero regression coefficient are not
excluded in the counting.
}
\value{
Either a number or a named vector of numbers (names being the node
names).
}
