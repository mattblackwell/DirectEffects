% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/telescope_match.R
\name{balance.tmatch}
\alias{balance.tmatch}
\title{Balance diagnostics for Telescope Match objects}
\usage{
balance.tmatch(object, vars, data, comparison = NULL)
}
\arguments{
\item{object}{an object of class \code{tmatch} -- results from a
call to \code{telescope_match}}

\item{vars}{a formula object containing either the treatment or the
mediator as the dependent variable (which denotes whether
first-stage or second-stage balance diagnostics are returned) and
the covariates for which balance diagnostics are requested as the
independent variables. Each covariate or function of covariates
(e.g. higher-order polynomials or interactions) should be
separated by a +.}

\item{data}{the data frame used in the call to
\code{telescope_match}}

\item{comparison}{a binary indicator for if the function should
return the balance for the treated group (`1`), for the control
group (`0`), or for overall combined balanced (`NULL`, the
default).}
}
\value{
Returns a data frame with the following columns.
\itemize{

\item variable: Name of covariate

\item before_0: Pre-matching average of the covariate in the
mediator == 0 (if first stage balance) or treatment == 0 (if second
stage balance) condition

\item before_1: Pre-matching average of the covariate in the
mediator == 1 (if first stage balance) or treatment == 1 (if second
stage balance) condition

\item after_0: Post-matching average of the covariate in the
mediator == 0 (if first stage balance) or treatment == 0 (if second
stage balance) condition

\item after_1: Post-matching average of the covariate in the
mediator == 1 (if first stage balance) or treatment == 1 (if second
stage balance) condition

\item before_sd: standard deviation of the outcome (pre-Matching)

\item before_diff: Pre-matching covariate difference between
mediator arms (if first stage balance) or treatment arms (if second
stage balance).

\item before_std_diff: Pre-matching standardized covariate
difference between mediator arms (if first stage balance) or
treatment arms (if second stage balance), Equal to Before_Diff/SD.

\item after_diff: Post--matching covariate difference between
mediator arms (if first stage balance) or treatment arms (if second
stage balance).

\item after_std_diff: Post-matching standardized covariate
difference between mediator arms (if first stage balance) or
treatment arms (if second stage balance), Equal to Before_Diff/SD.
}
}
\description{
Balance diagnostics for Telescope Match objects
}
\details{
Provides matching balance diagnostics for \code{tmatch}
objects returned by \code{telescope_match}
}
