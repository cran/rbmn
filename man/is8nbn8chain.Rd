\name{is8nbn8chain}
\alias{is8nbn8chain}
\title{
Checks if a given /nbn/ is a /chain/
}
\description{
returns \code{TRUE} [the order] or \code{FALSE} [NULL] according that
\code{nbn} is a chain of not [according to \code{order}].
}
\usage{
is8nbn8chain(nbn, order=FALSE)
}
\examples{
is8nbn8chain(rbmn0nbn.01);
is8nbn8chain(rbmn0nbn.04);
}
\arguments{
\item{nbn}{ The nbn object to consider.}
\item{order}{ When \code{FALSE} the answer to the question is
returned with \code{TRUE} or \code{FALSE}.\cr When \code{TRUE} the
chain order of the nodes is returned if it is a /chain/ else
\code{NULL}.}
}
\value{
A \code{logical(1)} when \code{order} si \code{TRUE} if not the
resulting chain order versus NULL.
}
