\name{cmb.em}
\alias{cmb.em}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conditional mixture modeling by EM algorithm}
\description{Runs conditional mixture modeling and model-based clustering by EM algorithm (Expectation Maximization) for a prespecified variables conditioning order. Runs variable selection procedure (forward, backward or stepwise) to achieve a parsimonious mixture model.}

\usage{
cmb.em(x, order = NULL, l, K, method = "stepwise", id0 = NULL, n.em = 200, em.iter = 5,
EM.iter = 200, nk.min = NULL, max.spur=5, tol = 1e-06, silent = FALSE, Parallel = FALSE,
n.cores = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x }{ dataset matrix (n x p)}
  \item{order }{ customized variables' conditioning order (length p)}
  \item{l }{ order of polynomial regression model }
  \item{K }{ number of clusters}
  \item{method }{variable selection method (options 'stepwise', 'forward', 'backward' and 'none')}
  \item{id0 }{initial membership vector (length n)}
  \item{n.em }{number of short EM in an emEM procedure}
  \item{em.iter }{maximum number of iterations of short EM in an emEM procedure}
  \item{EM.iter }{ maximum number of EM iterations}
  \item{nk.min}{spurious output control}
  \item{max.spur}{number of trials}
  \item{tol}{tolerance level}
  \item{silent}{output control (TRUE/FALSE)}
  \item{Parallel}{parallel computing (TRUE/FALSE)}
  \item{n.cores}{number of cores in parallel computing}
}
\details{
In conditional mixture modeling, each component is modeled by a product of conditional distributions with the means expressed by polynomial regression functions depending on other variables. Polynomial regression function order \code{l} and the number of clusters \code{K} are prespecified by user. The model's initialization can be determined by passing a group membership vector to the argument \code{id}, or obtained by the emEM algorithm (the default setting) in the function. There are two arguments related to the emEM procedure, the number of short EM \code{n.em} and maximum number of iterations for short EM \code{em.iter}. By default, the \code{n.em = 200} and \code{em.iter = 5}. The method of variable selection can be specified as \code{method = "stepwise", "forward", "backward", or "none"} where \code{method = none} means no parsimonious procedure conducted. During the model fitting and variable selection phases, EM algorithm will be applied multiple times, where options \code{EM.iter} and \code{tol} are stopping criteria of EM iteration. The spurious output control argument \code{nk.min}, by default \code{nk.min = (l x (p - 1) + 1) x 2}, can be set by user. When spurious output is obtained, \code{cmb.em} will be rerun. The maximum number of rerunning is \code{max.spur}. 

Notation: n - sample size, l - order of polynomial regression model, K - number of mixture components.
}
\value{
 \item{data}{input dataset}
  \item{model}{estimated regression models for each cluster (K x p matrix) }
  \item{id}{ vector of estimated membership (length n)}
  \item{loglik}{ estimated log likelihood}
  \item{BIC }{ Bayesian Information Criterion}
  \item{Pi   }{vector of estimated mixing proportions (length K)}
  \item{tau }{ matrix of estimated posterior probabilities (n x K)}
  \item{beta }{ matrix of estimated regression parameters (K x (p + p(p-1)l/2) )}
  \item{s2 }{ matrix of estimated variance (K x p)}
  \item{order}{applied conditioning order (length p)}
  \item{n_pars }{ number of model parameters}
}

\references{
Biernacki C., Celeux G., Govaert G. (2003). Choosing Starting Values for the EM Algorithm for Getting the Highest Likelihood in Multivariate Gaussian Mixture Models. \emph{Computational Statistics and Data Analysis}, 41(3-4), pp. 561-575.
}


\examples{
set.seed(1)
K <- 3
l <- 2
x <- as.matrix(iris[,-5])
id.true <- iris[,5]
\donttest{
# Run EM algorithm for fitting a conditioning mixture model 
obj <- cmb.em(x = x, order = c(1,3,2,4), l, K, method = "stepwise", silent = FALSE,
Parallel = FALSE)
id.cmb <- obj$id
table(id.true, id.cmb)
obj$BIC
}
}

\keyword{ Model-based clustering }
\keyword{conditional mixture modeling}
\keyword{ EM algorithm }
\keyword{ variable selection}
