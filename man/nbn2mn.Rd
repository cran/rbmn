\name{nbn2mn}
\alias{nbn2mn}
\title{
computes the joint distribution of a /nbn/
}
\description{
Computes the joint distribution of a /nbn/ with three possible
algorithms according to \code{algo}.
}
\usage{
nbn2mn(nbn, algo=3)
}
\examples{
 print8mn(nbn2mn(rbmn0nbn.05));
}
\arguments{
\item{nbn}{The \code{nbn} object to be converted.}
\item{algo}{either \code{1}: transforming the \code{nbn} into a
\code{gema} first before getting the \code{mn} form; or \code{2}: one
variable after another is added to the joint distribution following a
topological order; or \code{3}: variances are computed through the
differents paths o}
}
\details{
To be explained if it works
}
\value{
the resulting /mn/ object
}
