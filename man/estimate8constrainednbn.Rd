\name{estimate8constrainednbn}
\alias{estimate8constrainednbn}
\title{
estimates the parameters of a nbn with equality constraints
}
\description{
Estimations of the parameters of a /nbn/ is done when there are some
equality constraints onto the regression coefficients.\cr Constant
terms (\code{mu}) and conditional standard deviations (\code{sigma})
are supposed independent (that is not constrained with
equalities).\cr Equality constraints are given by \code{sarc}, a list
of matrices with two columns, indicating each the series of arcs
having the same regression coefficient.
}
\usage{
estimate8constrainednbn(nbn, sarc, data, imp=0, nite=10, eps=10^-5)
}
\examples{
 data(boco);
 print8nbn(rbmn0nbn.05);
 print8nbn(estimate8nbn(rbmn0nbn.05, boco));
 print8nbn(estimate8constrainednbn(rbmn0nbn.05, rbmn0crarc.05, boco));
}
\arguments{
\item{nbn}{\code{nbn} object.}
\item{sarc}{List of Matrices with two columns indicating the tails
(1rst column) and the heads (2d column) of the arcs having a common
parameter. It is checked that these arcs are indeed included in
\code{nbn}. Nodes must be indicated by their names (not their
number).}
\item{data}{Data frame to be used for the estimation. It must
comprise all necessary nodes (not only those involved in \code{sarc}
but also the remaining parents of \code{sarc[,2]}. Usually, all used
variables are centred but this is not required.}
\item{imp}{When \code{0} nothing displayed. When \code{1} the number
of iterations is displayed. When \code{2} the successive values of
the criterion are also displayed. }
\item{nite}{Maximum number of iterations.}
\item{eps}{relative difference in successive scores needed to stop
the iterations.}
}
\details{
Not linked regression coefficients doesn't require to be included in
\code{sarc}, the function do it by itself.\cr The score to use to
measure the differences between two successive estimations is not
well established (see the code).
}
\value{
The resulting /nbn/ object with the estimated parameters.
}
