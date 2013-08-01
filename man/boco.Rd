\name{body composition}
\docType{data}
\alias{boco}
\title{Body Composition Variables and Covariables }
\description{

  Real-world data set extracted from the Nhanes data base comprising
  nine variables describing the body composition and five easy
  measurable covariables.

}
\usage{
data(boco)
}
\format{

  The \code{boco} data set stored in variable \code{boco}
  comprises 100 individuals with the following variables:
  
  \itemize{
    \item \code{A} the age in years
    \item \code{H} the height in cm
    \item \code{W} the weight in kg
    \item \code{C} the waist circumference in cm
    \item \code{TF} the trunk fat in kg
    \item \code{LF} the leg fat in kg
    \item \code{AF} the arm fat in kg
    \item \code{TL} the trunk lean in kg
    \item \code{LL} the leg lean in kg
    \item \code{AL} the arm lean in kg
    \item \code{TB} the trunk bone in kg
    \item \code{LB} the leg bone in kg
    \item \code{AB} the arm bone in kg
  }
}

\source{

  Centers for Disease Control and Prevention. The 1999-2004 dual energy
  X-ray absorptiometry (DXA) multiple imputation data files and
  technical documentation.\cr Available from:
  http://www.cdc.gov/nchs/about/major/nhanes/dxx/dxa.html (accessed on
  13_07_03).
}

\examples{

# load the data and build the correct network from the model string.
data(boco);
print(head(boco));
boco7dag <- "[H][W|H][TF|W;H]";
# to be finished
}

\keyword{datasets}
