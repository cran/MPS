\name{gexppg}
\alias{dgexppg} \alias{pgexppg} \alias{qgexppg} \alias{rgexppg} \alias{mpsgexppg}
\title{geometric exponential Poisson G distribution}
\description{Computes the pdf, cdf, quantile, and random numbers of the geometric exponential Poisson \code{G} distribution. General form for the probability density function (pdf) of the geometric exponential Poisson \code{G} distribution due to Nadarajah et al. (2013) is given by
\deqn{f(x,{\Theta}) = \frac{{a(1 - b)\,g(x-\mu,\theta )\left( {1 - {e^{ - a}}} \right){e^{ - a + a\,G(x-\mu,\theta )}}}}{{{{\left( {1 - {e^{ - a}} - b + b{e^{ - a + a\,G(x-\mu,\theta )}}} \right)}^2}}},}
where \eqn{\theta} is the baseline family parameter vector. Also, a>0, 0<b<1, and \eqn{\mu} are the extra parameters induced to the baseline cumulative distribution function (cdf) \code{G} whose pdf is \code{g}. The general form for the cumulative distribution function (cdf) of the geometric exponential Poisson \code{G} distribution is given by
\deqn{F(x,{\Theta}) = \frac{{{e^{ - a + aG(x-\mu,\theta )}} - {e^{ - a}}}}{{{{{1 - {e^{ - a}} - b + b{e^{ - a + aG(x-\mu,\theta )}}}}}}}.}
Here, the baseline \code{G} refers to the cdf of famous families such as: Birnbaum-Saunders, Burr type XII, Exponential, Chen, Chisquare, F, Frechet, Gamma, Gompertz, Linear failure rate (lfr), Log-normal, Log-logistic, Lomax, Rayleigh, and Weibull. The parameter vector is \eqn{\Theta=(a,b,\theta,\mu)} where \eqn{\theta} is the baseline \code{G} family's parameter space. If \eqn{\theta} consists of the shape and scale parameters, the last component of \eqn{\theta} is the scale parameter (here, a and b are the first and second shape parameters). Always, the location parameter \eqn{\mu} is placed in the last component of \eqn{\Theta}.}
\usage{
dgexppg(mydata, g, param, location = TRUE, log=FALSE)
pgexppg(mydata, g, param, location = TRUE, log.p = FALSE, lower.tail = TRUE)
qgexppg(p, g, param, location = TRUE, log.p = FALSE, lower.tail = TRUE)
rgexppg(n, g, param, location = TRUE)
mpsgexppg(mydata, g, location = TRUE, method, sig.level)
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

Nadarajah, S., Cancho, V. G., and Ortega, E. M. M. (2013). The geometric exponential Poisson distribution, \emph{Statistical Methods & Applications}, 22, 355-380.}
\author{Mahdi Teimouri}
\examples{
x<-rweibull(100,shape=2,scale=2)+3
dgexppg(x, "weibull", c(1,0.5,2,2,3))
pgexppg(x, "weibull", c(1,0.5,2,2,3))
qgexppg(runif(100), "weibull", c(1,0.5,2,2,3))
rgexppg(100, "weibull", c(1,0.5,2,2,3))
mpsgexppg(x, "weibull", TRUE, "Nelder-Mead", 0.05)
}
