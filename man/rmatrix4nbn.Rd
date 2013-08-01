\name{rmatrix4nbn}
\alias{rmatrix4nbn}
\title{
regression matrix of a /nbn/
}
\description{
returns a dimnamed matrix indicating with \code{rho} an arc from row
to column nodes (0 everywhere else) where \code{rho} is the
regression coefficient. Also conditional standard deviations can be
introduced as diagonal elements but \code{mu} coefficient are lost...
It is advisable to normalize the /nbn/ first.
}
\usage{
rmatrix4nbn(nbn, stdev=TRUE)
}
\examples{
 rmatrix4nbn(rbmn0nbn.02);
 (rmatrix4nbn(rbmn0nbn.02, FALSE)>0)*1;
}
\arguments{
\item{nbn}{ The initial \code{nbn} object.}
\item{stdev}{Indicates if the standard deviations must placed in the
diagonal positions.}
}
\value{
A dimnamed matrix
}
