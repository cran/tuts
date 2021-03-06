% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/b_tupoisBFS.R
\name{summary.tuts_poisBFS}
\alias{summary.tuts_poisBFS}
\title{Summary tables of tuts_poisBFS objects}
\usage{
\method{summary}{tuts_poisBFS}(object, ...)
}
\arguments{
\item{object}{A tuts_poisBFS object.}

\item{...}{A list of optional parameters: \cr
- burn: burn-in parameter ranging from 0 to 0.7, the default value is 0.\cr
- CI: confidence interval, the default value is set to 0.99.}
}
\description{
\code{summary.tuts_poisBFS} prints summary tables of tuts_poisBFS objects.
}
\examples{
# Note: Most of models included in tuts package are computationally intensive. In the example
# below I set parameters to meet CRAN’s testing requirement of maximum 5 sec per example.
# A more practical example would contain N=50 in the first line of the code and n.sim=10000.

 #1. Import or simulate the data (simulation is chosen for illustrative purposes):
DATA=simtuts(N=7,Harmonics=c(4,0,0), sin.ampl=c(10,0, 0), cos.ampl=c(0,0,0),
            trend=0,y.sd=2, ti.sd=0.2)
y=DATA$observed$y.obs
y=round(y-min(y))
ti.mu=DATA$observed$ti.obs.tnorm
ti.sd= rep(0.2, length(ti.mu))

#2. Fit the model:
n.sim=10
TUPOIS=tupoisbsf(y=y,ti.mu=ti.mu,ti.sd=ti.sd,freqs='internal',n.sim=n.sim,n.chains=2, CV=FALSE)

#3. Generate summary results (optional parameters are listed in brackets):
summary(TUPOIS)                               # Summary results (CI, burn).
summary(TUPOIS,burn=0.2)                      # Results after 20\% of burn-in (CI).

}
