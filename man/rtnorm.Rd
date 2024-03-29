% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/truncated_normal.R
\name{rtnorm}
\alias{rtnorm}
\title{Univariate truncated normal sampler}
\usage{
rtnorm(n = 1, mean = 0, sd = 1, a = -Inf, b = Inf)
}
\arguments{
\item{n}{[integer] sample size}

\item{mean}{[numeric] scalar location parameter}

\item{sd}{[numeric] scalar scale parameter}

\item{a}{[numeric] lower bound}

\item{b}{[numeric] upper bound}
}
\description{
This implementations works for a < 37
or b > 37
}
