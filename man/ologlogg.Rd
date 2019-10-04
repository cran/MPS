\name{ologlogg}
\alias{dologlogg} \alias{pologlogg} \alias{qologlogg} \alias{rologlogg} \alias{mpsologlogg} \alias{qqologlogg}
\title{odd log-logistic G distribution}
\description{Computes the pdf, cdf, quantile, and random numbers, draws the q-q plot, and estimates the parameters of the odd log-logistic \code{G} distribution. General form for the probability density function (pdf) of the odd log-logistic \code{G} distribution due to Gauss et al. (2017) is given by
\deqn{f(x,{\Theta}) = \frac{{a\,b\,d\,g(x-\mu,\theta ){{\left( {G(x-\mu,\theta )} \right)}^{a\,d - 1}}{{\left[ {\bar G(x-\mu,\theta )} \right]}^{d - 1}}}}{{{{\left[ {{{\left( {G(x-\mu,\theta )} \right)}^d} - {{\left( {\bar G(x-\mu,\theta )} \right)}^d}} \right]}^{a + 1}}}}{\left\{ {1 - {{\left[ {\frac{{{{\left( {G(x-\mu,\theta )} \right)}^d}}}{{{{\left( {G(x-\mu,\theta )} \right)}^d} - {{\left( {\bar G(x-\mu,\theta )} \right)}^d}}}} \right]}^a}} \right\}^{b - 1}},}
with \eqn{\bar G(x-\mu,\theta ) = 1 - G(x-\mu,\theta )} where \eqn{\theta} is the baseline family parameter vector. Also, a>0, b>0, d>0, and \eqn{\mu} are the extra parameters induced to the baseline cumulative distribution function (cdf) \code{G} whose pdf is \code{g}. The general form for the cumulative distribution function (cdf) of the odd log-logistic \code{G} distribution is given by
\deqn{F(x,{\Theta}) = 1 - {\left\{ {1 - {{\left[ {\frac{{{{\left( {G(x-\mu,\theta )} \right)}^d}}}{{{{\left( {G(x-\mu,\theta )} \right)}^d} - {{\left( {\bar G(x-\mu,\theta )} \right)}^d}}}} \right]}^a}} \right\}^b}.}
Here, the baseline \code{G} refers to the cdf of famous families such as: Birnbaum-Saunders, Burr type XII, Exponential, Chen, Chisquare, F, Frechet, Gamma, Gompertz, Linear failure rate (lfr), Log-normal, Log-logistic, Lomax, Rayleigh, and Weibull. The parameter vector is \eqn{\Theta=(a,b,d,\theta,\mu)} where \eqn{\theta} is the baseline \code{G} family's parameter space. If \eqn{\theta} consists of the shape and scale parameters, the last component of \eqn{\theta} is the scale parameter (here, a, b, and d are the first, second, and the third shape parameters). Always, the location parameter \eqn{\mu} is placed in the last component of \eqn{\Theta}.}
\usage{
dologlogg(mydata, g, param, location = TRUE, log=FALSE)
pologlogg(mydata, g, param, location = TRUE, log.p = FALSE, lower.tail = TRUE)
qologlogg(p, g, param, location = TRUE, log.p = FALSE, lower.tail = TRUE)
rologlogg(n, g, param, location = TRUE)
qqologlogg(mydata, g, location = TRUE, method)
mpsologlogg(mydata, g, location = TRUE, method, sig.level)
}
\arguments{
  \item{g}{The name of family's pdf including: "\code{birnbaum-saunders}", "\code{burrxii}", "\code{chisq}", "\code{chen}", "\code{exp}", "\code{f}", "\code{frechet}", "\code{gamma}", "\code{gompetrz}", "\code{lfr}", "\code{log-normal}", "\code{log-logistic}", "\code{lomax}", "\code{rayleigh}", and "\code{weibull}".}
  \item{p}{a vector of value(s) between 0 and 1 at which the quantile needs to be computed.}
  \item{n}{number of realizations to be generated.}
  \item{mydata}{Vector of observations.}
  \item{param}{parameter vector \eqn{\Theta=(a,b,d,\theta,\mu)}}
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

Gauss, M. C., Alizadeh, M., Ozel, G., Hosseini, B. Ortega, E. M. M., and Altunc, E. (2017). The generalized odd log-logistic family of distributions: properties, regression models and applications, \emph{Journal of Statistical Computation and Simulation}, 87(5), 908-932.}
\author{Mahdi Teimouri}
\examples{
mydata<-rweibull(100,shape=2,scale=2)+3
dologlogg(mydata, "weibull", c(1,1,1,2,2,3))
pologlogg(mydata, "weibull", c(1,1,1,2,2,3))
qologlogg(runif(100), "weibull", c(1,1,1,2,2,3))
rologlogg(100, "weibull", c(1,1,1,2,2,3))
qqologlogg(mydata, "weibull", TRUE, "Nelder-Mead")
mpsologlogg(mydata, "weibull", TRUE, "Nelder-Mead", 0.05)
}

