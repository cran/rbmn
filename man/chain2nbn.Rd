\name{chain2nbn}
\alias{chain2nbn}
\title{
transforms a /chain/ to a /nbn/
}
\description{
From a \code{chain} object returns the \code{nbn} translation.
}
\usage{
chain2nbn(chain)
}
\examples{
 print8nbn(chain2nbn(rbmn0chain.02), ordering=names(rbmn0nbn.02));
}
\arguments{
\item{chain}{ the \code{chain} object to be transformed.}
}
\value{
The corresponding \code{nbn} object.
}
