\name{bnfit2nbn}
\alias{bnfit2nbn}
\title{
transforms a /bn.fit/ of /bnlearn/ package to a /nbn/
}
\description{
returns a \code{nbn} object from a Gaussian \code{bn.fit} object of
/bnlearn/ package.
}
\usage{
bnfit2nbn(bn.fit)
}
\arguments{
\item{bn.fit}{ The object to be transformed.}
}
\details{
If \code{bn.fit} is not pertinent, a fatal error is issued.
}
\value{
A list following the \code{nbn} specification
}
