\name{diff8nbn}
\alias{diff8nbn}
\title{
returns a score of the difference between two /nbn/s
}
\description{
Returns a positive scalar value measuring, in some way, the
difference existing within two /nbn/s sharing the same structure.\cr
}
\usage{
diff8nbn(nbn1, nbn2, type=1, scalar=TRUE)
}
\examples{
 diff8nbn(rbmn0nbn.01, rbmn0nbn.01);
 diff8nbn(rbmn0nbn.01, rbmn0nbn.01, scalar=FALSE);
}
\arguments{
\item{nbn1}{First \code{nbn} object.}
\item{nbn2}{Second \code{nbn} object.}
\item{type}{When 1, the score includes the difference between the
sigmas. When -1, sigmas are not taken into account.}
\item{scalar}{When \code{TRUE} the squared norm is returned, if not
the vector of difference.}
}
\details{
For \code{type==1} it is the canonical euclidian difference between
all parameters, including the \code{sigma}. The score to use to
measure the differences between two successive estimations is not
well established (see the code).
}
\value{
Either a scalar or a named vector (according to \code{scalar}).
}
