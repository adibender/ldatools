% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formula-parser.R
\name{vget_multappear}
\alias{vget_multappear}
\title{Check multiple strings for multiple appearances}
\usage{
vget_multappear(charvec, vec)
}
\arguments{
\item{charvec}{A Vector of strings. Each will be checked for multiple appearences
in \code{vec}}

\item{vec}{The character vector in which to check multiple appearences of \code{char}.}
}
\value{
\code{TRUE} if any element of \code{charvec} appears in \code{vec}
more than once, \code{FALSE} otherwise.
}
\description{
A vectorized version of \code{get_multappear}.
}
