\name{kumg}
\alias{dkumg} \alias{pkumg} \alias{qkumg} \alias{rkumg} \alias{mpskumg} \alias{qqkumg}
\title{Kumaraswamy G distribution}
\description{Computes the pdf, cdf, quantile, and random numbers, draws the q-q plot, and estimates the parameters of the Kumaraswamy \code{G} distribution. General form for the probability density function (pdf) of the Kumaraswamy \code{G} distribution due to Cordeiro and Castro (2011) is given by
\deqn{f(x,{\Theta}) = a\,b\,g(x-\mu,\theta ){\left( {G(x-\mu,\theta )} \right)^{a - 1}}{\left[ {1 - {{\left( {G(x-\mu,\theta )} \right)}^a}} \right]^{b - 1}},}
where \eqn{\theta} is the baseline family parameter vector. Also, a>0, b>0, and \eqn{\mu} are the extra parameters induced to the baseline cumulative distribution function (cdf) \code{G} whose pdf is \code{g}. The general form for the cumulative distribution function (cdf) of the Kumaraswamy \code{G} distribution is given by
\deqn{F(x,{\Theta}) = 1 - {\left[ {1 - {{\left( {G(x-\mu,\theta )} \right)}^a}} \right]^b}.}
Here, the baseline \code{G} refers to the cdf of famous families such as: Birnbaum-Saunders, Burr type XII, Exponential, Chen, Chisquare, F, Frechet, Gamma, Gompertz, Linear failure rate (lfr), Log-normal, Log-logistic, Lomax, Rayleigh, and Weibull. The parameter vector is \eqn{\Theta=(a,b,\theta,\mu)} where \eqn{\theta} is the baseline \code{G} family's parameter space. If \eqn{\theta} consists of the shape and scale parameters, the last component of \eqn{\theta} is the scale parameter (here, a and b are the first and second shape parameters). Always, the location parameter \eqn{\mu} is placed in the last component of \eqn{\Theta}.}
\usage{
dkumg(mydata, g, param, location = TRUE, log=FALSE)
pkumg(mydata, g, param, location = TRUE, log.p = FALSE, lower.tail = TRUE)
qkumg(p, g, param, location = TRUE, log.p = FALSE, lower.tail = TRUE)
rkumg(n, g, param, location = TRUE)
qqkumg(mydata, g, location = TRUE, method)
mpskumg(mydata, g, location = TRUE, method, sig.level)
}
\arguments{
  \item{g}{The name of family's pdf including: "\code{birnbaum-saunders}", "\code{burrxii}", "\code{chisq}", "\code{chen}", "\code{exp}", "\code{f}", "\code{frechet}", "\code{gamma}", "\code{gompetrz}", "\code{lfr}", "\code{log-normal}", "\code{log-logistic}", "\code{lomax}", "\code{rayleigh}", and "\code{weibull}".}
  \item{p}{a vector of value(s) between 0 and 1 at which the quantile needs to be computed.}
  \item{n}{number of realizations to be generated.}
  \item{mydata}{Vector of observations.}
  \item{param}{parameter vector \eqn{\Theta=(a,b,\theta,\mu)}}
  \item{location}{If \code{FALSE}, then the location parameter will be omitted.}
  \item{log}{If \code{TRUE}, then log(pdf) is returned.}
  \item{log.p}{If \code{TRUE}, then log(cdf) is returned and quantile is computed for \code{exp(-p)}.}
  \item{lower.tail}{If \code{FALSE}, then \code{1-cdf} is returned and quantile is computed for \code{1-p}.}
  \item{method}{The used method for maximizing the sum of log-spacing function. It will be  "\code{BFGS}", "\code{CG}", "\code{L-BFGS-B}", "\code{Nelder-Mead}", or "\code{SANN}".}
  \item{sig.level}{Significance level for the Chi-square goodness-of-fit test.}
}
\details{
It can be shown that the Moran's statistic follows a normal distribution. Also, a chi-square approximation exists for small samples whose mean and variance approximately are m(\code{log}(m)+0.57722)-0.5-1/(12m) and m(\eqn{\pi^2}/6-1)-0.5-1/(6m), respectively, with \code{m=n+1}, see Cheng and Stephens (1989). So, a hypothesis tesing can be constructed based on a sample of \code{n} independent realizations at the given significance level, indicated in above as \code{sig.level}.}
\value{
\enumerate{
\item A vector of the same length as \code{mydata}, giving the pdf values computed at \code{mydata}.
\item A vector of the same length as \code{mydata}, giving the cdf values computed at \code{mydata}.
\item A vector of the same length as \code{p}, giving the quantile values computed at \code{p}.
\item A vector of the same length as \code{n}, giving the random numbers realizations.
\item A sequence of goodness-of-fit statistics such as: Akaike Information Criterion (\code{AIC}), Consistent Akaike Information Criterion (\code{CAIC}), Bayesian Information Criterion (\code{BIC}), Hannan-Quinn information criterion (\code{HQIC}), Cramer-von Misses statistic (\code{CM}), Anderson Darling statistic (\code{AD}), log-likelihood statistic (\code{log}), and Moran's statistic (\code{M}). The Kolmogorov-Smirnov (\code{KS}) test statistic and corresponding \code{p-value}. The Chi-square test statistic, critical upper tail Chi-square distribution, related \code{p-value}, and the convergence status.
}
}
\references{
Cheng, R. C. H. and Stephens, M. A. (1989). A goodness-of-fit test using Moran's statistic with estimated parameters, \emph{Biometrika}, 76 (2), 385-392.

Cordeiro, G. M. and  Castro, M. (2011). A new family of generalized distributions, \emph{Journal of Statistical Computation and Simulation}, 81, 883-898.}
\author{Mahdi Teimouri}
\examples{
mydata<-rweibull(100,shape=2,scale=2)+3
dkumg(mydata, "weibull", c(1,1,2,2,3))
pkumg(mydata, "weibull", c(1,1,2,2,3))
qkumg(runif(100), "weibull", c(1,1,2,2,3))
rkumg(100, "weibull", c(1,1,2,2,3))
qqkumg(mydata, "weibull", TRUE, "Nelder-Mead")
mpskumg(mydata, "weibull", TRUE, "Nelder-Mead", 0.05)
}
