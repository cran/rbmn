\name{adja4nbn}
\alias{adja4nbn}
\title{
adjacency matrix of a /nbn/
}
\description{
returns a dimnamed matrix indicating with 1 an arc from row to column
nodes (0 everywhere else); i.e. the adjacency matrix.
}
\usage{
adja4nbn(nbn)
}
\examples{
 adja4nbn(rbmn0nbn.04);
}
\arguments{
\item{nbn}{ The initial \code{nbn} object.}
}
\value{
A dimnamed matrix
}
