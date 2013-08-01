\name{rm8nd4adja}
\alias{rm8nd4adja}
\title{
removes somes nodes from an adjacency matrix
}
\description{
Eliminates from the adjacency matrix (\code{adja})all \code{nodes}
not breaking the existing links.\cr Important: the node order in
\code{adja} must be topological.
}
\usage{
rm8nd4adja(adja, nodes)
}
\examples{
 rm8nd4adja(rbmn0adja.04, "1.1");
}
\arguments{
\item{adja}{The relation matrix to be consider (same format as those
provided by the function \code{adja4nbn}. Must be in topological
order, roots first.}
\item{nodes}{Numeric or character vector providing the node numbers
to use for the generation of the subset.}
}
\details{
When a node is removed, all its parents become parent of its
children.
}
\value{
The reduced adjacency matrix.
}
\keyword{utilities}
\keyword{PKEYWORDS}
