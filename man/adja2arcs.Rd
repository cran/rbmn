\name{adja2arcs}
\alias{adja2arcs}
\title{
Arc matrix from an adjacency matrix
}
\description{
returns the arc matrix from an adjacency matrix.
}
\usage{
adja2arcs(adj)
}
\examples{
 adja2arcs(rbmn0adja.02)
}
\arguments{
\item{adj}{ The adjacency matrix.}
}
\value{
a matrix with two columns ("from","to")
}
