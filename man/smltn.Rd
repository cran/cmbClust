\name{smltn}
\alias{smltn}
\alias{smltn1}
\alias{smltn2}
\docType{data}
\encoding{UTF-8}
\title{
Simulated datasets
}
\description{
Datasets are simulated from conditional mixture models with different numbers of components.
}
\usage{data(smltn)}

\format{
Two datasets are stored in the data \code{smltn}. \code{smltn1} is a data matrix with 200 observations on two variables and one group membership; \code{smltn2} is a matrix with 300 observations on two variables and one group ID.
}


\examples{
data(smltn)
# view data matrices smltn1 and smltn2
print(smltn1)
print(smltn2)
}
\keyword{datasets}
