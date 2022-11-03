\name{cmb.plot}
\alias{cmb.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphic display for the results of conditional mixture modeling and model-based clustering}
\description{
\code{cmb.plot} demonstrates the clustering results of functions \code{cmb.em} and \code{cmb.search}. A graph with a combination of pairwise scatter plot for data points, pairwise contour plot of estimated mixture density, and pairwise regression curves is produced.
}

\usage{
cmb.plot(obj, allcolors = NULL, allpch = NULL, lwd = 1, cex.text = 1, cex.point = 0.6, 
mar = c(0.6,0.6,0.6,0.6), oma = c(3.5,3.5,2.5,14), nlevels = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  	\item{obj }{output object of the function \code{cmb.em()} or \code{cmb.search()}}
		\item{allcolors}{colous of clusters (length K)}
		\item{allpch}{styles of data points in clusters (length K)}
		\item{lwd}{line width, a positive number, defaulting to 1}
		\item{cex.text}{magnification of labels and titles, defaulting to 1}		
		\item{cex.point}{magnification of plotting symbols, defaulting to 0.6}
		\item{mar}{margin sizes of plots in lines of text (length 4)}
		\item{oma}{outer margin sizes of a pairwise plot in lines of text (length 4)}
		\item{nlevels}{number of contour levels, defaulting to 30}
}

\value{
This function generates a graphic. 
}


\examples{
set.seed(4)
K <- 3
l <- 2
x <- as.matrix(iris[,-5])

# Run EM algorithm for fitting a conditioning mixture model 
\donttest{
obj <- cmb.em(x = x, order = c(1,2,3,4), l, K, method = "stepwise",
silent = TRUE, Parallel = FALSE)
cmb.plot(obj)
}	
}

\keyword{Pairwise plot}

