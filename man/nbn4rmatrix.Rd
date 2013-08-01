\name{nbn4rmatrix}
\alias{nbn4rmatrix}
\title{
a /nbn/ from a regression matrix
}
\description{
reverse of \code{rmatrix4nbn} but the standard deviations must be
included.
}
\usage{
nbn4rmatrix(rmatrix)
}
\examples{
 print8nbn(nbn4rmatrix(rmatrix4nbn(rbmn0nbn.02)));
}
\arguments{
\item{rmatrix}{ The regression coefficient matrix with the standard
deviations in the diagonal.}
}
\details{
\code{mu}s are put to nought
}
\value{
A /nbn/ object
}
