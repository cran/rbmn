\name{normalize8nbn}
\alias{normalize8nbn}
\title{
normalizes a /nbn/
}
\description{
returns a \code{nbn} with a given expectation and variance through an
transformation leaving the correlation unchanged.
}
\usage{
normalize8nbn(nbn, mu=0, sigma=1)
}
\examples{
 print8nbn(normalize8nbn(rbmn0nbn.01));
}
\arguments{
\item{nbn}{ The \code{nbn} object to transform.}
\item{mu}{ Imposed expectations. When \code{NULL} nothing is changed.
When of length one, this value is given to all the node expectations.
If not the complete vector of expect}
\item{sigma}{ The same as \code{mu} but for the standard deviations.}
}
\value{
The transformed \code{nbn}.
}
