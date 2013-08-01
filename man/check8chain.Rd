\name{check8chain}
\alias{check8chain}
\title{
checks a /chain/ object
}
\description{
checks the consistency of \code{chain} as a /chain/ object issues a
fatal error with some clues if inconsistent.
}
\usage{
check8chain(chain)
}
\examples{
 check8chain(rbmn0chain.01);
 res <- check8chain(rbmn0adja.01);
 if (is.na(as.logical(res))) { print(res);}
}
\arguments{
\item{chain}{ The \code{chain} object to check.}
}
\details{
Looking a the code of this function provides a way to know which are
the requirements of a /chain/ object.
}
\value{
\code{TRUE} or a \code{character} containing some clue about the
discovered inconsistency.
}
