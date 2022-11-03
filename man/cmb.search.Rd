\name{cmb.search}
\alias{cmb.search}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optimal conditioning order search}
\description{
Runs forward, backward, or stepwise variable selection procedure for obtaining the parsimonious conditional mixture models when all conditional orders are considered.
Alternatively, runs the optimal order search algorithm, and parsimonious conditional mixture modeling for the obtained order.}

\usage{
cmb.search(x, l, K, method = "stepwise", all.perms  = TRUE, id0 = NULL,  n.em = 200,
em.iter = 5, EM.iter = 200, nk.min = NULL, max.spur = 5, tol = 1e-06, silent = FALSE, 
Parallel = TRUE, n.cores = 4)}


\arguments{
  \item{x }{ dataset matrix (n x p)}
  \item{l }{ order of polynomial regression model }
  \item{K }{ number of clusters}
  \item{method }{  variable selection method (options 'stepwise', 'forward', 'backward' and 'none')}
  \item{all.perms}{conditioning order search algorithm (\code{TRUE}: full search; \code{FALSE} proposed search algorithm) }
  \item{id0 }{ initial group membership (length n)}
  \item{n.em }{ number of short EM in emEM procedure}
  \item{em.iter }{maximum number of short EM iterations in emEM}
  \item{EM.iter }{ maximum number of EM iterations}
  \item{nk.min}{spurious output control}
  \item{max.spur}{number of trials}
  \item{tol}{tolerance level}
  \item{silent}{output control}
  \item{Parallel}{Parallel computing}
  \item{n.cores}{number of cores in parallel computing}
}

\details{
Functions 'cmb.search' and 'cmb.em' have common arguments except 'all.perm'. With \code{all.perms = TRUE}, a full search is applied to data, that is running parsimonious conditional mixture modeling for all orders and recognizing the optimal order based on the BIC. Then two lists are returned: \code{best.model} stores the results for the conditional mixture model with the optimal order, and \code{models} has results for all orders. With the option \code{all.perms = FALSE}, the optimal conditional order search algorithm is applied, and then only the list \code{best.model} is returned. 
}

\value{
  The list \code{models} is returned when \code{all.perms = TRUE}.
  \item{best.model}{  membership assignments and estimated parameters of mixture model with the optimal contioning order.
  \describe{
  \item{data}{the input dataset}
  \item{model}{estimated regression models for each cluster (K x p) }
  \item{id}{ vector of estimated membership (length n)}
  \item{loglik}{ estimated log likelihood}
  \item{BIC }{Bayesian Information Criterion }
  \item{Pi   }{vector of estimated mixing proportions (length K)}
  \item{tau }{ matrix of estimated posterior probabilities  (n x K)}
  \item{beta }{matrix of estimated regression parameters (K x (p + p(p-1)l/2) )}
  \item{s2 }{ matrix of estimated variances (K x p)}
  \item{order}{applied conditioning order {length p}}
  \item{n_pars }{ number of parameters}
  }
  }

  \item{models}{ membership assignments and model parameters of mixture models with all conditioning orders.
  \describe{
  \item{model}{list of estimated regression models for all clusters (K x p x p!) }
  \item{id}{ p! vectors of estimated memberships (n x p!)}
  \item{loglik}{ estimated log likelihood values (length p!)}
  \item{BIC }{ Bayesian Information Criterion values (length p!)}
  \item{Pi   }{ p! vectors of estimated mixing proportions  (K x p!)}
  \item{tau }{ p! matrices of estimated posterior probabilities (K x p!)}
  \item{beta }{ p! matrices of estimated regression parameters  (K x (p + p(p-1)l/2) x p! )}
  \item{s2 }{ p! matrices of estimated variances (K x p x p!)}
  \item{order}{appled conditioning orders (p! x p )}
  \item{n_pars }{ number of parameters in p! models (length p!)}
  }}
}




\seealso{
cmb.em
}
\examples{
set.seed(1)
K = 3
l <- 2
x <- as.matrix(iris[,-5])

\donttest{
obj <- cmb.search(x = x, l, K, method = "stepwise", all.perms = FALSE,
Parallel = FALSE, silent = FALSE)
obj$best.model$BIC
}

}

\keyword{Optimal conditioning order search}

