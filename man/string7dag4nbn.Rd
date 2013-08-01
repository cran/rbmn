\name{string7dag4nbn}
\alias{string7dag4nbn}
\title{
provides so-called string model of a /nbn/
}
\description{
returns a \code{character(1)} describing the dag of the nbn under the
string form.
}
\usage{
string7dag4nbn(nbn, sep=";")
}
\examples{
 string7dag4nbn(rbmn0nbn.01);
 string7dag4nbn(rbmn0nbn.04, sep=", ");
}
\arguments{
\item{nbn}{ The nbn.}
\item{sep}{Separation sign between parents after the conditioning
sign (\code{|}).}
}
\value{
A \code{character(1)}.
}
