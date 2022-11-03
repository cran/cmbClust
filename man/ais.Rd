\name{ais}
\alias{ais}
\docType{data}
\encoding{UTF-8}
\title{Australian Institute of Sport data}
\description{
The data set considers body characteristics from 102 male and 100 female athletes at the Australian Institute of Sport. It is collected for a study of how data on various features varied with sport body size and sex of athlete.
}
\usage{data(ais)}
\format{ 
  A data frame with 202 observations on the following 13 variables.
  \describe{
   \item{sex}{Factor with levels: \code{female}, \code{male};}
   \item{sport}{Factor with levels: \code{B_Ball}, 
            \code{Field}, \code{Gym}, \code{Netball}, \code{Row}
            \code{Swim}, \code{T_400m}, \code{Tennis}, \code{T_Sprnt}, 
            \code{W_Polo};}
   \item{RCC}{Red cell count;}
   \item{WCC}{White cell count;}
   \item{Hc}{Hematocrit;}
   \item{Hg}{Hemoglobin;}
   \item{Fe}{Plasma ferritin concentration;}
   \item{BMI}{Body Mass Index;}
   \item{SSF}{Sum of skin folds;}
   \item{Bfat}{Body fat percentage;}
   \item{LBM}{Lean body mass;}
   \item{Ht}{Height, cm;}
   \item{Wt}{Weight, kg}
  }
}
\details{The data have been made publicly available in connection with the book by Cook, R.D. and Weisberg, S. (1994, ISBN-10:0471008397).}

\references{
Cook, R.D. and Weisberg, S. (1994). \emph{An introduction to regression graphics}. John Wiley & Sons.
}



\keyword{datasets}
