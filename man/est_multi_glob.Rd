\name{est_multi_glob}
\alias{est_multi_glob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit marginal regression models for categorical responses}
\description{It estimates marginal regression models to datasets consisting of a categorical response and one or more covariates by a Fisher-scoring algorithm; this is an internal function.}

\usage{est_multi_glob(yv, X, model, ind = rep(1, length(yv)), w = rep(1, length(yv)), 
	be = NULL)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yv}{vector of responses}
  \item{X}{matrix of all distinct covariate configurations}
  \item{model}{type of logit (g = global, l = local, m = multinomial)}
  \item{ind}{vector to link responses to covariates}
  \item{w}{vector of weights}
  \item{be}{initial vector of regression coefficients}
}

%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}

\value{
  \item{be}{estimated vector of regression coefficients}
  \item{lk}{log-likelihood at convergence}
  \item{P}{matrix of the probabilities for each covariate configuration}
}

\references{Colombi, R. and Forcina, A. (2001), Marginal regression models for the analysis of positive
association of ordinal response variables, \emph{Biometrika}, \bold{88}, 1007-1019.

Glonek, G. F. V. and McCullagh, P. (1995), Multivariate logistic models, \emph{Journal of the Royal 
Statistical Society, Series B}, \bold{57}, 533-546.}

\author{Francesco Bartolucci - University of Perugia (IT)}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{maximum likelihood estimation}
\keyword{Fisher-scoring algorithm}% __ONLY ONE__ keyword per line