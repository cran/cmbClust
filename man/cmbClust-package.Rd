\name{cmbClust-package}
\alias{cmbClust-package}
\docType{package}
\title{
Conditional mixture modeling and model-based clustering
}
\description{
The utility of this package includes fitting a conditional mixture model with EM (Expectation Maximization) algorithm, model-based clustering based on conditional mixture modeling, conditional mixture modeling with parsimonious procedures, and optimal conditional order exploration by using either a full search or the proposed searching algorithm, and illustration of clustering results through pairwise plots.}
\details{

Function 'cmb.em' runs the parsimonious conditional mixture modeling for a user-specified conditioning order.

Function 'cmb.search' runs the 'cmb.em' procedure for all possible conditioning orders and determines the optimal order using BIC, or runs the proposed optimal order search algorithm and then the 'cmb.em' for the obtained optimal order.

Function 'cmb.plot' builds pairwise plots to present clustering results from functions 'cmb.em' and 'cmb.search".
}
\author{Yang Wang and Volodymyr Melnykov.

Maintainer: Yang Wang <wangy4@cofc.edu>
}

\references{
Melnykov, V., and Wang, Y. (2023). Conditional mixture modeling and model-based clustering. \emph{Pattern Recognition}, 133, p. 108994.
}


