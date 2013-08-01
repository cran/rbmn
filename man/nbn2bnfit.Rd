\name{nbn2bnfit}
\alias{nbn2bnfit}
\title{
transforms a /nbn/ to a /bn.fit/ of /bnlearn/ package
}
\description{
returns a \code{bn.fit} object from a Gaussian \code{nbn} object of
/rbmn/ package.
}
\usage{
nbn2bnfit(nbn, onlydag=FALSE)
}
\arguments{
\item{nbn}{ The object to be transformed.}
\item{onlydag}{ Indicates if only the DAG must be computed. In that
case a /bn/ object of /bnlearn/}
}
\value{
The resulting \code{bn.fit} (or \code{bn}) object.
}
