\name{est_multilogit}
\alias{est_multilogit}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Estimate multilogit model}

\description{The function performs maximum likelihood estimation of the multilogit model (internal function).}

\usage{est_multilogit(Y, Xdis, label = 1:n, be = NULL, Pdis = NULL, dis = F, fort = F)}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{Y}{matrix of all responses}
  \item{Xdis}{matrix of all covariates (two- or three-dimensional)}
  \item{label}{label associated to every covariate configuration}
  \item{be}{initial value of the parameter vector}
  \item{Pdis}{probability configurations as initial values}
  \item{dis}{to dispaly partial results}
  \item{fort}{to use fortran routines when possible}
}

%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}

\value{
\item{be}{estimated vector regression coefficients}
\item{P}{estimated matrix of probabilities}
\item{Pdis}{estimated matrix of distinct probabilities}
}

\author{Francesco Bartolucci - University of Perugia (IT)}

%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}