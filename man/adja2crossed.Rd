\name{adja2crossed}
\alias{adja2crossed}
\title{
creates a crossed-adjacency matrix from two ones
}
\description{
Like crossed4nbn1nbn but at the level of adjacency matrices. Must be
much efficient when regression coefficients are not needed.
}
\usage{
adja2crossed(adj1, adj2, nona=as.vector(outer(dimnames(adj1)[[1]],
	  dimnames(adj2)[[1]], paste, sep="_")))
}
\examples{
 print(adja2crossed(rbmn0adja.01, rbmn0adja.01));
}
\arguments{
\item{adj1}{ The first adjacency matrix.}
\item{adj2}{ The second adjacency matrix.}
\item{nona}{ The node names to give to the crossed /nbn/, the nodes
of the \code{nbn1} varying first.}
}
\details{
Just two Kronecker products of matrices.
}
\value{
The resulting crossed adjacency matrix.
}
