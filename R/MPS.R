mpsbetaexpg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    starts<-c(1,1,1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
      starts<-c(1,1,1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
      starts<-c(1,1,1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    starts<-c(1,1,1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
      starts<-c(1,1,1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
      starts<-c(1,1,1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
      starts<-c(1,1,1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    d*f0*beta(a,b)^(-1)*(1-c0)^(d*a-1)*(1-(1-c0)^d)^(b-1)}

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    1-pbeta((1-c0)**d,shape1=a,shape2=b)}

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }

  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  if (location==TRUE){
    out$par[1:(n.p-1)]<-abs(out$par[1:(n.p-1)])}
  else{
    out$par[1:n.p]<-abs(out$par[1:n.p])
  }

  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }

  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}


mpsbetag<-function(mydata, g, location=TRUE, method, sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    starts<-c(1,1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
      starts<-c(1,1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
      starts<-c(1,1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    starts<-c(1,1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
      starts<-c(1,1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    f0*beta(a,b)^(-1)*c0^(a-1)*(1-c0)^(b-1)}

  cdf0<-function(par,x){
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    pbeta(c0,shape1=a,shape2=b)}

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log", "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsexpexppg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    starts<-c(1,1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
      starts<-c(1,1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
      starts<-c(1,1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    starts<-c(1,1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
      starts<-c(1,1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    a*b*(1-exp(-b))^(-1)*f0*c0^(a-1)*exp(-b*c0^a)}

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    (1-exp(-b))^(-1)*(1-exp(-b*c0^a))}

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsexpgg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    starts<-c(1,1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
      starts<-c(1,1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
      starts<-c(1,1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    starts<-c(1,1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
      starts<-c(1,1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    r=par[1]
    u=par[2]
    r*u*f0*((1-c0)^(r-1))*(1-(1-c0)^r)^(u-1)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    r=par[1]
    u=par[2]
    (1-(1-c0)^r)^u
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsexpg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    starts<-c(1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
      starts<-c(1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
      starts<-c(1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    starts<-c(1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
      starts<-c(1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
      starts<-c(1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
      starts<-c(1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
      starts<-c(1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    a*f0*c0^(a-1)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    c0^a
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsexpkumg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    starts<-c(1,1,1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
      starts<-c(1,1,1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
      starts<-c(1,1,1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    starts<-c(1,1,1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
      starts<-c(1,1,1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
      starts<-c(1,1,1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
      starts<-c(1,1,1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    a*b*d*f0*c0^(a-1)*(1-c0^a)^(b-1)*(1-(1-c0^a)^b)^(d-1)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    (1-(1-c0^a)^b)^d
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }

  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  if (location==TRUE){
    out$par[1:(n.p-1)]<-abs(out$par[1:(n.p-1)])}
  else{
    out$par[1:n.p]<-abs(out$par[1:n.p])
  }

  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }

  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsgammag<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    starts<-c(1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
      starts<-c(1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
      starts<-c(1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    starts<-c(1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
      starts<-c(1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
      starts<-c(1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
      starts<-c(1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
      starts<-c(1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    gamma(a)^(-1)*f0*(1-c0)^(-2)*(c0/(1-c0))^(a-1)*exp(-c0/(1-c0))
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    pgamma(c0/(1-c0),shape=a)
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsgammag1<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    starts<-c(1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
      starts<-c(1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
      starts<-c(1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    starts<-c(1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
      starts<-c(1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
      starts<-c(1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
      starts<-c(1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
      starts<-c(1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    gamma(a)^(-1)*f0*(-log(1-c0))^(a-1)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    pgamma(-log(1-c0),shape=a)
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsgammag2<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    starts<-c(1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
      starts<-c(1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
      starts<-c(1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    starts<-c(1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
      starts<-c(1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
      starts<-c(1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
      starts<-c(1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
      starts<-c(1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    gamma(a)^(-1)*f0*(-log(c0))^(a-1)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    pgamma(-log(c0),shape=a)
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsgbetag<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    starts<-c(1,1,1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
      starts<-c(1,1,1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
      starts<-c(1,1,1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    starts<-c(1,1,1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
      starts<-c(1,1,1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
      starts<-c(1,1,1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
      starts<-c(1,1,1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    d*beta(a,b)^(-1)*f0*c0^(a*d-1)*(1-c0^d)^(b-1)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    pbeta(c0^d,shape1=a,shape2=b)
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }

  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  if (location==TRUE){
    out$par[1:(n.p-1)]<-abs(out$par[1:(n.p-1)])}
  else{
    out$par[1:n.p]<-abs(out$par[1:n.p])
  }

  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }

  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}

mpsgtransg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,0.5,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,0.5,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    starts<-c(1,0.5,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
      starts<-c(1,0.5,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,0.5,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
      starts<-c(1,0.5,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,0.5,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,0.5,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,0.5,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,0.5,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
      starts<-c(1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    starts<-c(1,0.5,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
      starts<-c(1,0.5,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,0.5,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
      starts<-c(1,0.5,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,0.5,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
      starts<-c(1,0.5,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,0.5,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,0.5,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,0.5,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
      starts<-c(1,0.5,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,0.5,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
      starts<-c(1,0.5,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    a*f0*(c0)**(a-1)*(1+b-2*b*c0)*(1+b*(1-c0))**(a-1)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    (c0)**a*(1+b*(1-c0))**(a)
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpskumg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    starts<-c(1,1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
      starts<-c(1,1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
      starts<-c(1,1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    starts<-c(1,1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
      starts<-c(1,1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    r=par[1]
    u=par[2]
    r*u*f0*(c0^(r-1))*(1-c0^r)^(u-1)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    r=par[1]
    u=par[2]
    1-(1-c0^r)^u
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsologlogg<-function(mydata, g, location=TRUE, method, sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    starts<-c(1,1,1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
      starts<-c(1,1,1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
      starts<-c(1,1,1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    starts<-c(1,1,1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
      starts<-c(1,1,1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
      starts<-c(1,1,1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
      starts<-c(1,1,1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    a*b*d*f0*c0**(d*a-1)*(1-c0)**(d-1)*(c0**d+(1-c0)**d)**(-a-1)*(1-(c0**d/(c0**d+(1-c0)**d))**a)**(b-1)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    1-(1-(c0**d/(c0**d+(1-c0)**d))**a)**b
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }

  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)

  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }

  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsloggammag1<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    starts<-c(1,1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
      starts<-c(1,1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
      starts<-c(1,1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    starts<-c(1,1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
      starts<-c(1,1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    b^a*gamma(a)^(-1)*f0*(-log(1-c0))^(a-1)*(1-c0)^(b-1)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    pgamma(-b*log(1-c0),shape=a)
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsloggammag2<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    starts<-c(1,1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
      starts<-c(1,1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
      starts<-c(1,1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    starts<-c(1,1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
      starts<-c(1,1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    b^a*gamma(a)^(-1)*f0*(-log(c0))^(a-1)*(c0)^(b-1)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    1-pgamma(-b*log(c0),shape=a)
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsgxlogisticg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    starts<-c(1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
      starts<-c(1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
      starts<-c(1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    starts<-c(1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
      starts<-c(1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
      starts<-c(1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
      starts<-c(1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
      starts<-c(1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    a*f0/(1-c0)*(-log(1-c0))**(-a-1)*(1+(-log(1-c0))**(-a))**(-2)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    1/(1+(-log(1-c0))**(-a))
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsgmbetaexpg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    starts<-c(1,1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
      starts<-c(1,1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
      starts<-c(1,1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    starts<-c(1,1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
      starts<-c(1,1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    b*a*f0/(1-c0)**2*exp(-b*c0/(1-c0))*(1-exp(-b*c0/(1-c0)))**(a-1)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    (1-exp(-b*c0/(1-c0)))**a
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpstexpsg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    starts<-c(1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
      starts<-c(1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
      starts<-c(1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    starts<-c(1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
      starts<-c(1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
      starts<-c(1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
      starts<-c(1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
      starts<-c(1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    a*(1-exp(-a))^(-1)*f0*exp(-a*c0)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    (1-exp(-a*c0))*(1-exp(-a))^(-1)
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsweibullg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    starts<-c(1,1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
      starts<-c(1,1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
      starts<-c(1,1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    starts<-c(1,1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
      starts<-c(1,1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    a*b**(-a)*(f0/(1-c0))*(-log(1-c0))**(a-1)*exp(-b**(-a)*(-log(1-c0))**a)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    1-exp(-b**(-a)*(-log(1-c0))**a)
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}

mpsmbetag<-function(mydata, g, location=TRUE, method, sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,0.5,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,0.5,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    starts<-c(1,1,0.5,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
      starts<-c(1,1,0.5,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,0.5,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
      starts<-c(1,1,0.5,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,0.5,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,0.5,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,0.5,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,0.5,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    starts<-c(1,1,0.5,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
      starts<-c(1,1,0.5,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,0.5,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,0.5,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,0.5,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,0.5,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,0.5,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,0.5,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,0.5,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
      starts<-c(1,1,0.5,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,0.5,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
      starts<-c(1,1,0.5,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    d^a*f0*beta(a,b)^(-1)*c0^(a-1)*(1-c0)^(b-1)*(1-(1-d)*c0)^(-a-b)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    pbeta(d*c0/(1-(1-d)*c0),shape1=a,shape2=b)
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }

  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  if (location==TRUE){
    out$par[1:(n.p-1)]<-abs(out$par[1:(n.p-1)])}
  else{
    out$par[1:n.p]<-abs(out$par[1:n.p])
  }

  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }

  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}

mpsmog<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(0.5,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    starts<-c(1,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
      starts<-c(1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
      starts<-c(1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    starts<-c(1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
      starts<-c(1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
      starts<-c(1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
      starts<-c(1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
      starts<-c(1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
      starts<-c(0.5,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    r=par[1]
    r*f0/((1-(1-r)*(1-c0))^2)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    r=par[1]
    1-r*(1-c0)/(1-(1-r)*(1-c0))
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
mpsmokumg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,0.5,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,0.5,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    starts<-c(1,1,0.5,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
      starts<-c(1,1,0.5,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,0.5,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
      starts<-c(1,1,0.5,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,0.5,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,0.5,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-c(seq(1,(n-1)))
    yy<--log(1-(ii)/(n+1))
    res<-suppressWarnings(summary(lm(yy ~ -1+y+I(y^2)))$coefficient)
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,0.5,p.hat[1],p.hat[2],(min.mydata+1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,0.5,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,0.5,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    starts<-c(1,1,0.5,1/mean.mydata,(min.mydata+1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
      starts<-c(1,1,0.5,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,0.5,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,0.5,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,0.5,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,0.5,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,0.5,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,0.5,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,0.5,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
      starts<-c(1,1,0.5,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,0.5,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
      starts<-c(1,1,0.5,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    (d*a*b*c0^(a-1)*(1-c0^a)^(b-1))/((1-(1-a)*(1-c0^a)^b)^2)
  }
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    1-(d*(1-c0^a)^b)/(1-(1-d)*(1-c0^a)^b)
  }
  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=suppressWarnings((2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i]))
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}


mpsgexppg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,0.5,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,0.5,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    starts<-c(1,0.5,mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
      starts<-c(1,0.5,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,0.5,df.1,df.2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
      starts<-c(1,0.5,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,0.5,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,0.5,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,0.5,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,0.5,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
      starts<-c(1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    starts<-c(1,0.5,1/mean.mydata,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
      starts<-c(1,0.5,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,0.5,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
      starts<-c(1,0.5,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,0.5,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,0.5,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,0.5,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
      starts<-c(1,0.5,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,0.5,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,0.5,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,0.5,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
      starts<-c(1,0.5,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,0.5,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
      starts<-c(1,0.5,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    theta=par[1]
    eta=par[2]
    theta*(1-eta)*(1-exp(-theta))*f0*exp(-theta+theta*c0)/(1-exp(-theta)-eta+eta*exp(-theta+theta*c0))**2
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    theta=par[1]
    eta=par[2]
    (exp(-theta+theta*c0)-exp(-theta))/(1-exp(-theta)-eta+eta*exp(-theta+theta*c0))
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
pbetaexpg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    1-pbeta((1-c0)**d,shape1=a,shape2=b)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dbetaexpg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    d*f0*beta(a,b)^(-1)*(1-c0)^(d*a-1)*(1-(1-c0)^d)^(b-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pbetag<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    pbeta(c0,shape1=a,shape2=b)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dbetag<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    f0*beta(a,b)^(-1)*c0^(a-1)*(1-c0)^(b-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pexpexppg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    (1-exp(-b))^(-1)*(1-exp(-b*c0^a))
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dexpexppg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    a*b*(1-exp(-b))^(-1)*f0*c0^(a-1)*exp(-b*c0^a)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pexpgg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    r=par[1]
    u=par[2]
    (1-(1-c0)^r)^u
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dexpgg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    r=par[1]
    u=par[2]
    r*u*f0*((1-c0)^(r-1))*(1-(1-c0)^r)^(u-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pexpg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    c0^a
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dexpg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    a*f0*c0^(a-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pexpkumg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    (1-(1-c0^a)^b)^d
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dexpkumg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    a*b*d*f0*c0^(a-1)*(1-c0^a)^(b-1)*(1-(1-c0^a)^b)^(d-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pgammag<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    pgamma(c0/(1-c0),shape=a)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dgammag<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    gamma(a)^(-1)*f0*(1-c0)^(-2)*(c0/(1-c0))^(a-1)*exp(-c0/(1-c0))
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pgammag1<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    pgamma(-log(1-c0),shape=a)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dgammag1<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    gamma(a)^(-1)*f0*(-log(1-c0))^(a-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pgammag2<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    pgamma(-log(c0),shape=a)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dgammag2<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    gamma(a)^(-1)*f0*(-log(c0))^(a-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pgbetag<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    pbeta(c0^d,shape1=a,shape2=b)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dgbetag<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    d*beta(a,b)^(-1)*f0*c0^(a*d-1)*(1-c0^d)^(b-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pgtransg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    (c0)**a*(1+b*(1-c0))**(a)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dgtransg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    a*f0*(c0)**(a-1)*(1+b-2*b*c0)*(1+b*(1-c0))**(a-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pkumg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    r=par[1]
    u=par[2]
    1-(1-c0^r)^u
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dkumg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    r=par[1]
    u=par[2]
    r*u*f0*(c0^(r-1))*(1-c0^r)^(u-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pologlogg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    1-(1-(c0**d/(c0**d+(1-c0)**d))**a)**b
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dologlogg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    a*b*d*f0*c0**(d*a-1)*(1-c0)**(d-1)*(c0**d+(1-c0)**d)**(-a-1)*(1-(c0**d/(c0**d+(1-c0)**d))**a)**(b-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
ploggammag1<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    pgamma(-b*log(1-c0),shape=a)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dloggammag1<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    b^a*gamma(a)^(-1)*f0*(-log(1-c0))^(a-1)*(1-c0)^(b-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
ploggammag2<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    1-pgamma(-b*log(c0),shape=a)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dloggammag2<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    b^a*gamma(a)^(-1)*f0*(-log(c0))^(a-1)*(c0)^(b-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pgxlogisticg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    1/(1+(-log(1-c0))**(-a))
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dgxlogisticg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    a*f0/(1-c0)*(-log(1-c0))**(-a-1)*(1+(-log(1-c0))**(-a))**(-2)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pgmbetaexpg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    (1-exp(-b*c0/(1-c0)))**a
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dgmbetaexpg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    b*a*f0/(1-c0)**2*exp(-b*c0/(1-c0))*(1-exp(-b*c0/(1-c0)))**(a-1)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
ptexpsg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    (1-exp(-a*c0))*(1-exp(-a))^(-1)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dtexpsg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    a*(1-exp(-a))^(-1)*f0*exp(-a*c0)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pweibullg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    1-exp(-b**(-a)*(-log(1-c0))**a)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dweibullg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    a*b**(-a)*(f0/(1-c0))*(-log(1-c0))**(a-1)*exp(-b**(-a)*(-log(1-c0))**a)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pmbetag<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    pbeta(d*c0/(1-(1-d)*c0),shape1=a,shape2=b)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dmbetag<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    d^a*f0*beta(a,b)^(-1)*c0^(a-1)*(1-c0)^(b-1)*(1-(1-d)*c0)^(-a-b)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pmog<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    r=par[1]
    1-r*(1-c0)/(1-(1-r)*(1-c0))
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dmog<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[2]; loc=par[3]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[2]; loc=par[3]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[2]; dchisq(x,df)}
      cum=function(par,x){df=par[2]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
      cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[2]; b=par[3]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[2]; loc=par[3]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[2]; loc=par[3]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[2]; dexp(x,rate)}
      cum=function(par,x){rate=par[2]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[2]; loc=par[3]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[2]; loc=par[3]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[2]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[2]; d=par[3]; loc=par[4]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[2]; d=par[3]; loc=par[4]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[2]; b=par[3]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[2]; b=par[3]; loc=par[4]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[2]; b=par[3]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[2]; b=par[3]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[2]; scale=par[3]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[2]; scale=par[3]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    r=par[1]
    r*f0/((1-(1-r)*(1-c0))^2)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pgexppg<-function(mydata, g, param, location = TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    theta=par[1]
    eta=par[2]
    (exp(-theta+theta*c0)-exp(-theta))/(1-exp(-theta)-eta+eta*exp(-theta+theta*c0))
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}

dgexppg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    theta=par[1]
    eta=par[2]
    theta*(1-eta)*(1-exp(-theta))*f0*exp(-theta+theta*c0)/(1-exp(-theta)-eta+eta*exp(-theta+theta*c0))**2
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}
pmokumg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    1-(d*(1-c0^a)^b)/(1-(1-d)*(1-c0^a)^b)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dmokumg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; plnorm(x-loc,meanlog,sdlog)}
    if (location==FALSE){
      den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[4]; loc=par[5]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[4]; loc=par[5]; pchisq(x-loc,df)}
    if (location==FALSE){
      den=function(par,x){df=par[4]; dchisq(x,df)}
      cum=function(par,x){df=par[4]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; pf(x-loc,df1,df2)}
    if (location==FALSE){
      den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
      cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1/(((x-loc)/b)^(-a)+1)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[4]; b=par[5]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-(1+a*(x-loc))^(-b)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[4]; loc=par[5]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[4]; loc=par[5]; pexp(x-loc,rate)}
    if (location==FALSE){
      den=function(par,x){rate=par[4]; dexp(x,rate)}
      cum=function(par,x){rate=par[4]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[4]; loc=par[5]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[4]; loc=par[5]; 1-exp(-((x-loc)/a)^2)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[4]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[4]; d=par[5]; loc=par[6]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[4]; d=par[5]; loc=par[6]; 1-(1+(x-loc)^d)^(-a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; exp(-((x-loc)/b)^(-a))}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[4]; b=par[5]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[4]; b=par[5]; loc=par[6]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location==FALSE){
      den=function(par,x){a=par[4]; b=par[5]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[4]; b=par[5]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pweibull(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; pgamma(x-loc,shape,scale)}
    if (location==FALSE){
      den=function(par,x){shape=par[4]; scale=par[5]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[4]; scale=par[5]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    d=par[3]
    (d*a*b*c0^(a-1)*(1-c0^a)^(b-1))/((1-(1-a)*(1-c0^a)^b)^2)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}

mpsweibullextg<-function(mydata, g, location=TRUE, method , sig.level){
  min.mydata<-min(mydata)
  sort.mydata<-sort(mydata)
  n<-length(mydata)
  y<-(sort.mydata[2:n]-min.mydata)
  y<-(sort.mydata[2:n]-min.mydata)[is.finite(log(sort.mydata[2:n]-min.mydata))]
  n<-length(mydata)
  qp1<-y[floor(.25*n)]
  qp3<-y[floor(.75*n)]
  median.mydata<-y[floor(.5*n)]
  mean.mydata<-mean(y)
  inv.mydata<-1/y[is.finite(1/y)]
  inv.mean<-mean(inv.mydata)
  std.mydata<-sd(y)
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))),(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
      starts<-c(1,1,log(median.mydata),sqrt(2*abs(log(mean.mydata/median.mydata))))
    }}

  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    starts<-c(1,1,mean.mydata,(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
      starts<-c(1,1,mean.mydata)
    }}

  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    df.1<-ifelse(inv.mean>1,2*inv.mean/(inv.mean-1),1)
    df.2<-ifelse(mean.mydata>1,2*mean.mydata/(mean.mydata-1),1)
    starts<-c(1,1,df.1,df.2,(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
      starts<-c(1,1,df.1,df.2)
    }}

  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    standard.mydata<-(y-mean.mydata)/std.mydata
    s.hat<-ifelse(max(standard.mydata)<=1,-log(1-n/(n+.4))/(exp(1)-1),-
                    log(1-length(standard.mydata[standard.mydata<=1])/length(standard.mydata))/(exp(1)-1))
    r.hat<-abs(log(log(1-log(1-0.5)/s.hat))/log(median.mydata))
    starts<-c(1,1,r.hat,s.hat,(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
      starts<-c(1,1,r.hat,s.hat)
    }}

  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    ii<-seq(1,(n-1))
    yy<--log(1-(ii+0.3)/(n-1+0.4))
    res<-summary(lm(yy ~ -1+y+I(y^2)))$coefficient
    coeff<-c(abs(res[1,1]),abs(res[2,1]))
    z1=NULL;
    lfr.log<-function(p) {
      z1<--log(sum(p[1]+p[2]*y))+p[1]*sum(y)+p[2]/2*sum(y^2)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(abs(coeff[1]),2*abs(coeff[2])),lfr.log)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata,(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
      starts<-c(1,1,log(0.75/(1-0.75))/log(qp3/median.mydata),median.mydata)
    }}

  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    a.hat<-1
    b.hat<-(.5^(-1/a.hat)-1)/median.mydata
    z1=NULL;
    lomax.log<-function(p) {
      z1<--n*log(p[1])-n*log(p[2])+(p[2]+1)*sum(log(1+p[1]*y))
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(1,b.hat),lomax.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    a.hat<-log(log(1-.75)/log(1-.5))/(qp3-median.mydata)
    b.hat<--a.hat*log(0.5)/(exp(a.hat*median.mydata)-1)
    z1=NULL;
    gompertz.log<-function(p){
      z1<--n*log(p[2])-p[1]*sum(y)+p[2]/p[1]*sum(exp(p[1]*y)-1)
      z1[z1<1e-16]<-1e-16
    }
    p.hat<-suppressWarnings(optim(c(a.hat,b.hat),gompertz.log, method)$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    starts<-c(1,1,1/mean.mydata,(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
      starts<-c(1,1,1/mean.mydata)
    }}

  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    starts<-c(1,1,log(2)/(median.mydata)^2,(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
      starts<-c(1,1,log(2)/(median.mydata)^2)
    }}

  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    z1=NULL;
    burrxii.log<-function(p){
      z1<--n*log(p[1])-n*log(p[2])-(p[2]-1)*sum(log(y))+(p[1]+1)*sum(log(1+y^p[2]))
    }
    p.hat<-suppressWarnings(optim(c(1,1),burrxii.log, method="BFGS")$par)
    starts<-c(1,1,p.hat[1],p.hat[2],(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
      starts<-c(1,1,p.hat[1],p.hat[2])
    }}

  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    cons<-cor(inv.mydata,rank(inv.mydata))*(sd(inv.mydata)/mean(inv.mydata))*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median(inv.mydata)/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,1/b.hat,(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
      starts<-c(1,1,a.hat,1/b.hat)
    }}

  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    r<-(mean(inv.mydata))^(-1)
    b.hat<-sqrt(mean.mydata*r)
    a.hat<-sqrt(2*(sqrt(mean.mydata/r)))
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    cons<-cor(mydata,rank(mydata))*(std.mydata/mean.mydata)*sqrt((n+1)/(n-1))/sqrt(3)
    a.hat<-abs(ifelse(cons<1,-log(2)/log(1-cons),-log(2)/log(1-.98)))
    b.hat<-median.mydata/(log(2))^(1/a.hat)
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    a.hat<-uniroot(function(ss) trigamma(ss)-var(log(y)[is.finite(log(y))]),c(0,1000000))$root
    b.hat<-mean.mydata/a.hat
    starts<-c(1,1,a.hat,b.hat,(min.mydata-1/length(y)))
    if (location=="FALSE"){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
      starts<-c(1,1,a.hat,b.hat)
    }}

  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    a*b*f0*c0**(a-1)/((1-c0)**(a+1))*exp(-b*(c0/(1-c0))**a)
  }

  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    1-exp(-b*(c0/(1-c0))**a)
  }

  mpsw<-function(par,x){
    z=NULL
    z=diff(c(0,cdf0(par,x),1))
    for (j in 2:n){if((x[j]-x[j-1])==0){z[j]=pdf0(par,x[j])
    }}
    z[z<1e-16]=1e-16
    -sum(log(z))
  }
  out<-suppressWarnings(optim(starts,fn=mpsw,x=sort.mydata,method=method))
  n.p<-length(out$par)
  gammam<-(n+1)*(log(n+1)-digamma(1))-1/2-1/(12*(n+1))
  sigma2m<-(n+1)*(pi^2/6-1)-1/2-1/(6*(n+1))
  c1<-gammam-sqrt(n*sigma2m/2)
  c2<-sqrt(sigma2m/(2*n))
  stat.chisquare<-(mpsw(out$par,sort.mydata)+n.p/2-c1)/c2
  pvalue.chisquare<-pchisq(stat.chisquare,df=n,lower.tail=FALSE)
  Moran<-out$value
  log.likelihood=sum(log(pdf0(out$par,sort.mydata)))
  u=cdf0(out$par,sort.mydata)
  von<-c()
  anderson<-c()
  for(i in 1:n){
    u[i]<-ifelse(u[i]==1,0.999999999,u[i])
    von[i]=(u[i]-(2*i-1)/(2*n))^2
    anderson[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
  }
  anderson.stat=-n-mean(anderson)
  von.stat=sum(von)+1/(12*n)
  CAIC=-2*log.likelihood + 2*n.p + 2*(n.p*(n.p+1))/(n-n.p-1)
  AIC=-2*log.likelihood + 2*n.p
  BIC=-2*log.likelihood + n.p*log(n)
  HQIC=-2*log.likelihood + 2*log(log(n))*n.p
  ks.stat=suppressWarnings(ks.test(mydata, "cdf0", par=out$par))

  aux2=cbind(AIC, CAIC, BIC, HQIC, von.stat, anderson.stat,log.likelihood,Moran)
  colnames(aux2)=c("AIC","CAIC","BIC","HQIC","CM","AD", "log",
                   "Moran")
  rownames(aux2)=c("")

  aux3=cbind(ks.stat$statistic,ks.stat$p.value)
  colnames(aux3)=c("statistic","p-value")
  rownames(aux3)=c("")

  aux4=cbind(stat.chisquare,qchisq(sig.level,df=n,lower.tail=FALSE),1-pchisq(stat.chisquare,df=n))
  colnames(aux4)=c("statistic","chi-value","p-value")
  rownames(aux4)=c("")

  aux5=cbind(if(out$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  list("MPS"=out$par,"Measures"=aux2,"KS"=aux3,"chi-square"=aux4,"Convergence Status"=aux5)
}
pweibullextg<-function(mydata, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location=="FALSE"){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location=="FALSE"){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location=="FALSE"){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location=="FALSE"){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location=="FALSE"){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location=="FALSE"){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  cdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    1-exp(-b*(c0/(1-c0))**a)
  }
  cdf<-cdf0(param,mydata)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdf)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdf)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdf
  return(cdf)
}
dweibullextg<-function(mydata, g, param, location=TRUE, log=FALSE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    den=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; dlnorm(x-loc,meanlog,sdlog)}
    cum=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; plnorm(x-loc,meanlog,sdlog)}
    if (location=="FALSE"){
      den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
      cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    den=function(par,x){df=par[3]; loc=par[4]; dchisq(x-loc,df)}
    cum=function(par,x){df=par[3]; loc=par[4]; pchisq(x-loc,df)}
    if (location=="FALSE"){
      den=function(par,x){df=par[3]; dchisq(x,df)}
      cum=function(par,x){df=par[3]; pchisq(x,df)}
    }}
  if(g=="f"){
    den=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; df(x-loc,df1,df2)}
    cum=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; pf(x-loc,df1,df2)}
    if (location=="FALSE"){
      den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
      cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
    }}
  if(g=="chen"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; a*b*(x-loc)^(a-1)*exp((x-loc)^a)*exp(-b*(exp((x-loc)^a)-1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-b*(exp((x-loc)^a)-1))}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; a*b*(x)^(a-1)*exp((x)^a)*exp(-b*(exp((x)^a)-1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp((x)^a)-1))}
    }}
  if(g=="lfr"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a+b*(x-loc))*exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-a*(x-loc)-(b*(x-loc)^2)/2)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (a+b*(x))*exp(-a*(x)-(b*(x)^2)/2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*(x)-(b*(x)^2)/2)}
    }}
  if(g=="log-logistic"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b^(-a)*(x-loc)^(a-1))/((((x-loc)/b)^a +1)^2)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1/(((x-loc)/b)^(-a)+1)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*(x)^(a-1))/((((x)/b)^a +1)^2)}
      cum=function(par,x){a=par[3]; b=par[4]; 1/(((x)/b)^(-a)+1)}
    }}
  if(g=="lomax"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*b)/((1+a*(x-loc))^(b+1))}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-(1+a*(x-loc))^(-b)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*(x))^(b+1))}
      cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*(x))^(-b)}
    }}
  if(g=="gompertz"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*exp(a*(x-loc))*exp(-(exp(a*(x-loc))-1)*b/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; 1-exp(-(exp(a*(x-loc))-1)*b/a)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
      cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*(x))-1)*b/a)}
    }}
  if(g=="exp"){
    den=function(par,x){rate=par[3]; loc=par[4]; dexp(x-loc,rate)}
    cum=function(par,x){rate=par[3]; loc=par[4]; pexp(x-loc,rate)}
    if (location=="FALSE"){
      den=function(par,x){rate=par[3]; dexp(x,rate)}
      cum=function(par,x){rate=par[3]; pexp(x,rate)}
    }}
  if(g=="rayleigh"){
    den=function(par,x){a=par[3]; loc=par[4]; 2*(x-loc)/a*exp(-((x-loc)/a)^2)}
    cum=function(par,x){a=par[3]; loc=par[4]; 1-exp(-((x-loc)/a)^2)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; 2*x/a*exp(-(x/a)^2)}
      cum=function(par,x){a=par[3]; 1-exp(-(x/a)^2)}
    }}
  if(g=="burrxii"){
    den=function(par,x){a=par[3]; d=par[4]; loc=par[5]; d*a*(1+(x-loc)^d)^(-a-1)*((x-loc)^(d-1))}
    cum=function(par,x){a=par[3]; d=par[4]; loc=par[5]; 1-(1+(x-loc)^d)^(-a)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*((x)^(d-1))}
      cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
    }}
  if(g=="frechet"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (a*exp(-((x-loc)/b)^(-a))*((x-loc)/b)^(-a-1))/(b)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; exp(-((x-loc)/b)^(-a))}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (a*exp(-((x)/b)^(-a))*((x)/b)^(-a-1))/(b)}
      cum=function(par,x){a=par[3]; b=par[4]; exp(-((x)/b)^(-a))}
    }}
  if(g=="birnbaum-saunders"){
    den=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (sqrt((x-loc)/b)+sqrt(b/(x-loc)))/(2*a*(x-loc))*dnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    cum=function(par,x){a=par[3]; b=par[4]; loc=par[5]; pnorm((sqrt((x-loc)/b)-sqrt(b/(x-loc)))/a)}
    if (location=="FALSE"){
      den=function(par,x){a=par[3]; b=par[4]; (sqrt((x)/b)+sqrt(b/(x)))/(2*a*(x))*dnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
      cum=function(par,x){a=par[3]; b=par[4]; pnorm((sqrt((x)/b)-sqrt(b/(x)))/a)}
    }}
  if(g=="weibull"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dweibull(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pweibull(x-loc,shape,scale)}
    if (location=="FALSE"){
      den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    den=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; dgamma(x-loc,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; pgamma(x-loc,shape,scale)}
    if (location=="FALSE"){
      den=function(par,x){shape=par[3]; scale=par[4]; dgamma(x,shape,scale)}
      cum=function(par,x){shape=par[3]; scale=par[4]; pgamma(x,shape,scale)}
    }}
  pdf0<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    a=par[1]
    b=par[2]
    a*b*f0*c0**(a-1)/((1-c0)**(a+1))*exp(-b*(c0/(1-c0))**a)
  }
  pdf<-pdf0(param,mydata)
  if(log==TRUE){pdf<-log(pdf)}
  return(pdf)
}

qmokumg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[4]; loc=par[5]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[4]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[4]; df2=par[5]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[4]; loc=par[5]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[4]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[4]; loc=par[5]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[4]; d=par[5]; loc=par[6]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; d=par[5]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    d=par[3]
    quan0(par,(1-((1-x)/(d+(1-d)*(1-x)))^(1/b))^(1/a))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rmokumg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[4]; loc=par[5]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[4]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[4]; df2=par[5]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[4]; loc=par[5]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[4]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[4]; loc=par[5]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[4]; d=par[5]; loc=par[6]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; d=par[5]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    d=par[3]
    quan0(par,(1-((1-x)/(d+(1-d)*(1-x)))^(1/b))^(1/a))
  }
  return(quan(param,runif(n)))
}
qmbetag<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[4]; loc=par[5]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[4]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[4]; df2=par[5]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[4]; loc=par[5]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[4]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[4]; loc=par[5]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[4]; d=par[5]; loc=par[6]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; d=par[5]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    d=par[3]
    quan0(par,qbeta(x,shape1=a,shape2=b)/(d+(1-d)*qbeta(x,shape1=a,shape2=b)))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rmbetag<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[4]; loc=par[5]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[4]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[4]; df2=par[5]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[4]; loc=par[5]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[4]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[4]; loc=par[5]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[4]; d=par[5]; loc=par[6]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; d=par[5]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    d=par[3]
    quan0(par,qbeta(x,shape1=a,shape2=b)/(d+(1-d)*qbeta(x,shape1=a,shape2=b)))
  }
  return(quan(param,runif(n)))
}
qgbetag<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[4]; loc=par[5]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[4]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[4]; df2=par[5]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[4]; loc=par[5]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[4]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[4]; loc=par[5]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[4]; d=par[5]; loc=par[6]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; d=par[5]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    d=par[3]
    quan0(par,(qbeta(x,shape1=a,shape2=b))^(1/d))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rgbetag<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[4]; loc=par[5]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[4]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[4]; df2=par[5]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[4]; loc=par[5]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[4]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[4]; loc=par[5]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[4]; d=par[5]; loc=par[6]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; d=par[5]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    d=par[3]
    quan0(par,(qbeta(x,shape1=a,shape2=b))^(1/d))
  }
  return(quan(param,runif(n)))
}
qexpkumg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[4]; loc=par[5]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[4]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[4]; df2=par[5]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[4]; loc=par[5]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[4]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[4]; loc=par[5]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[4]; d=par[5]; loc=par[6]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; d=par[5]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    d=par[3]
    quan0(par,(1-(1-x^(1/d))^(1/b))^(1/a))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rexpkumg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[4]; loc=par[5]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[4]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[4]; df2=par[5]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[4]; loc=par[5]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[4]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[4]; loc=par[5]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[4]; d=par[5]; loc=par[6]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; d=par[5]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    d=par[3]
    quan0(par,(1-(1-x^(1/d))^(1/b))^(1/a))
  }
  return(quan(param,runif(n)))
}
qbetaexpg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[4]; loc=par[5]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[4]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[4]; df2=par[5]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[4]; loc=par[5]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[4]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[4]; loc=par[5]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[4]; d=par[5]; loc=par[6]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; d=par[5]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    d=par[3]
    quan0(par,1-(qbeta(x,shape1=a,shape2=b))^(1/d))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rbetaexpg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[4]; loc=par[5]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[4]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[4]; df2=par[5]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[4]; loc=par[5]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[4]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[4]; loc=par[5]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[4]; d=par[5]; loc=par[6]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; d=par[5]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    d=par[3]
    quan0(par,1-(qbeta(x,shape1=a,shape2=b))^(1/d))
  }
  return(quan(param,runif(n)))
}
qologlogg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[4]; loc=par[5]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[4]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[4]; df2=par[5]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[4]; loc=par[5]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[4]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[4]; loc=par[5]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[4]; d=par[5]; loc=par[6]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; d=par[5]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    d=par[3]
    quan0(par,(((1-(1-(1-x)^(1/b))^(1/a))/(1-(1-x)^(1/b))^(1/a))^(1/d))/(1+((1-(1-(1-x)^(1/b))^(1/a))/(1-(1-x)^(1/b))^(1/a))^(1/d)))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rologlogg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; loc=par[6]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[4]; sdlog=par[5]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[4]; loc=par[5]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[4]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[4]; df2=par[5]; loc=par[6]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[4]; df2=par[5]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[4]; loc=par[5]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[4]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[4]; loc=par[5]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[4]; d=par[5]; loc=par[6]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; d=par[5]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[4]; b=par[5]; loc=par[6]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[4]; b=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[4]; scale=par[5]; loc=par[6]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[4]; scale=par[5]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    d=par[3]
    quan0(par,(((1-(1-(1-x)^(1/b))^(1/a))/(1-(1-x)^(1/b))^(1/a))^(1/d))/(1+((1-(1-(1-x)^(1/b))^(1/a))/(1-(1-x)^(1/b))^(1/a))^(1/d)))
  }
  return(quan(param,runif(n)))
}
qweibullextg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,(-log(1-x)/b)^(1/a)/(1+(-log(1-x)/b)^(1/a)))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rweibullextg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,(-log(1-x)/b)^(1/a)/(1+(-log(1-x)/b)^(1/a)))
  }
  return(quan(param,runif(n)))
}
qloggammag2<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,exp(-qgamma(1-x,shape=a,scale=1)/b))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rloggammag2<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,exp(-qgamma(1-x,shape=a,scale=1)/b))
  }
  return(quan(param,runif(n)))
}
qexpexppg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,(-log(1-x*(1-exp(-b)))/b)^(1/a))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rexpexppg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,(-log(1-x*(1-exp(-b)))/b)^(1/a))
  }
  return(quan(param,runif(n)))
}
qexpgg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,1-(1-x^(1/b))^(1/a))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rexpgg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,1-(1-x^(1/b))^(1/a))
  }
  return(quan(param,runif(n)))
}
qweibullg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,1-exp(-b*(-log(1-x))^(1/a)))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rweibullg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,1-exp(-b*(-log(1-x))^(1/a)))
  }
  return(quan(param,runif(n)))
}
qgmbetaexpg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,(-log(1-x^(1/a))/b)/(1+(-log(1-x^(1/a))/b)))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rgmbetaexpg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,(-log(1-x^(1/a))/b)/(1+(-log(1-x^(1/a))/b)))
  }
  return(quan(param,runif(n)))
}
qgtransg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,(1+b-sqrt((1+b)^2-4*b*x^(1/a)))/(2*b))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rgtransg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,(1+b-sqrt((1+b)^2-4*b*x^(1/a)))/(2*b))
  }
  return(quan(param,runif(n)))
}
qbetag<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,qbeta(x,shape1=a,shape2=b))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rbetag<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,qbeta(x,shape1=a,shape2=b))
  }
  return(quan(param,runif(n)))
}
qgexppg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    theta=par[1]
    eta=par[2]
    quan0(par,1+log(((1-exp(-theta))*x-x*eta+exp(-theta))/(1-x*eta))/theta)
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rgexppg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    theta=par[1]
    eta=par[2]
    quan0(par,1+log(((1-exp(-theta))*x-x*eta+exp(-theta))/(1-x*eta))/theta)
  }
  return(quan(param,runif(n)))
}
qgxlogisticg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par,1-exp(-((1-x)/x)^(-1/a)))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rgxlogisticg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par,1-exp(-((1-x)/x)^(-1/a)))
  }
  return(quan(param,runif(n)))
}
qmog<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par,1-(1-x)/(a+(1-x)*(1-a)))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rmog<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par,1-(1-x)/(a+(1-x)*(1-a)))
  }
  return(quan(param,runif(n)))
}
qgammag<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par,qgamma(x,shape=a)/(1+qgamma(x,shape=a)))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rgammag<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par,qgamma(x,shape=a)/(1+qgamma(x,shape=a)))
  }
  return(quan(param,runif(n)))
}
qgammag1<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par, 1-exp(-qgamma(x,shape=a)))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rgammag1<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par, 1-exp(-qgamma(x,shape=a)))
  }
  return(quan(param,runif(n)))
}
qgammag2<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par, exp(-qgamma(x,shape=a)))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rgammag2<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par, exp(-qgamma(x,shape=a)))
  }
  return(quan(param,runif(n)))
}
qexpg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par,x^(1/a))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rexpg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par,x^(1/a))
  }
  return(quan(param,runif(n)))
}
qtexpsg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par,-log(1-x*(1-exp(-a)))/a)
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rtexpsg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; loc=par[4]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[2]; sdlog=par[3]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[2]; loc=par[3]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[2]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[2]; df2=par[3]; loc=par[4]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[2]; df2=par[3]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[2]; loc=par[3]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[2]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[2]; loc=par[3]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[2]; d=par[3]; loc=par[4]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; d=par[3]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[2]; b=par[3]; loc=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[2]; b=par[3]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[2]; scale=par[3]; loc=par[4]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[2]; scale=par[3]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    quan0(par,-log(1-x*(1-exp(-a)))/a)
  }
  return(quan(param,runif(n)))
}
qkumg<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,(1-(1-x)^(1/b))^(1/a))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rkumg<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,(1-(1-x)^(1/b))^(1/a))
  }
  return(quan(param,runif(n)))
}
qloggammag1<-function(p, g, param, location=TRUE, log.p = FALSE, lower.tail = TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,1-exp(-qgamma(x,shape=a,scale=1)/b))
  }
  out<-quan(param,p)
  if(log.p==TRUE & lower.tail == FALSE) out<-quan(param,exp(p-1))
  if(log.p==TRUE & lower.tail == TRUE) out<-quan(param,exp(-p))
  if(log.p==FALSE & lower.tail == FALSE) out<-quan(param,1-p)
  return(out)
}

rloggammag1<-function(n, g, param, location=TRUE){
  if(g!="birnbaum-saunders" & g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma"
     & g!="log-normal" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet"
     & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
  { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }
  if(n %% 1 !=0) { stop ("parameter n must be integer.") }
  if(g=="log-normal"){
    quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; loc=par[5]; qlnorm(x,meanlog,sdlog)+loc}
    if (location==FALSE){
      quan0=function(par,x){meanlog=par[3]; sdlog=par[4]; qlnorm(x,meanlog,sdlog)}
    }}
  if(g=="chisq"){
    quan0=function(par,x){df=par[3]; loc=par[4]; qchisq(x,df)+loc}
    if (location==FALSE){
      quan0=function(par,x){df=par[3]; qchisq(x,df)}
    }}
  if(g=="f"){
    quan0=function(par,x){df1=par[3]; df2=par[4]; loc=par[5]; qf(x,df1,df2)+loc}
    if (location==FALSE){
      quan0=function(par,x){df1=par[3]; df2=par[4]; qf(x,df1,df2)}
    }}
  if(g=="chen"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (log(1-log(1-x)/b))^(1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (log(1-log(1-x)/b))^(1/a)}
    }}
  if(g=="lfr"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; (-a+sqrt(a^2-2*b*log(1-x)))/b+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; (-a+sqrt(a^2-2*b*log(1-x)))/b}
    }}
  if(g=="log-logistic"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(1/x-1)^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(1/x-1)^(-1/a)}
    }}
  if(g=="lomax"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; ((1-x)^(-1/b)-1)/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; ((1-x)^(-1/b)-1)/a}
    }}
  if(g=="gompertz"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; log(1-a/b*log(1-x))/a+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; log(1-a/b*log(1-x))/a}
    }}
  if(g=="exp"){
    quan0=function(par,x){rate=par[3]; loc=par[4]; qexp(x,rate)+loc}
    if (location==FALSE){
      quan0=function(par,x){rate=par[3]; qexp(x,rate)}
    }}
  if(g=="rayleigh"){
    quan0=function(par,x){a=par[3]; loc=par[4]; a*(-log(1-x))^(1/2)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; a*(-log(1-x))^(1/2)}
    }}
  if(g=="burrxii"){
    quan0=function(par,x){a=par[3]; d=par[4]; loc=par[5]; (1/(1-x)-1)^(1/d)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; d=par[4]; (1/(1-x)-1)^(1/d)}
    }}
  if(g=="frechet"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b*(-log(x))^(-1/a)+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b*(-log(x))^(-1/a)}
    }}
  if(g=="birnbaum-saunders"){
    quan0=function(par,x){a=par[3]; b=par[4]; loc=par[5]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2+loc}
    if (location==FALSE){
      quan0=function(par,x){a=par[3]; b=par[4]; b/4*(a*qnorm(x)+sqrt((a*qnorm(x))^2+4))^2}
    }}
  if(g=="weibull"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qweibull(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qweibull(x,shape,scale)}
    }}
  if(g=="gamma"){
    quan0=function(par,x){shape=par[3]; scale=par[4]; loc=par[5]; qgamma(x,shape,scale)+loc}
    if (location==FALSE){
      quan0=function(par,x){shape=par[3]; scale=par[4]; qgamma(x,shape,scale)}
    }}
  quan<-function(par,x){
    a=par[1]
    b=par[2]
    quan0(par,1-exp(-qgamma(x,shape=a,scale=1)/b))
  }
  return(quan(param,runif(n)))
}


qqbetaexpg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsbetaexpg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qbetaexpg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqbetag<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsbetag(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qbetag(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqexpexppg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsexpexppg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qexpexppg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqexpg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsexpg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qexpg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqexpgg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsexpgg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qexpgg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqexpkumg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsexpkumg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qexpkumg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqgammag<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsgammag(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qgammag(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqgammag1<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsgammag1(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qgammag1(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqgammag2<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsgammag2(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qgammag2(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqgbetag<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsgbetag(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qgbetag(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqgexppg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsgexppg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qgexppg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqgmbetaexpg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsgmbetaexpg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qgmbetaexpg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqgtransg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsgtransg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qgtransg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqgxlogisticg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsgxlogisticg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qgxlogisticg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqkumg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpskumg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qkumg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqloggammag1<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsloggammag1(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qloggammag1(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqloggammag2<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsloggammag2(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qloggammag2(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqmbetag<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsmbetag(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qmbetag(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqmog<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsmog(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qmog(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqmokumg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsmokumg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qmokumg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqologlogg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsologlogg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qologlogg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqtexpsg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpstexpsg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qtexpsg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqweibullextg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsweibullextg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qweibullextg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}

qqweibullg<-function(mydata, g, location = TRUE, method){
  if (g != "birnbaum-saunders" & g != "exp" & g != "rayleigh" &
      g != "weibull" & g != "gompertz" & g != "gamma" & g !=
      "log-normal" & g != "chisq" & g != "f" & g != "burrxii" &
      g != "frechet" & g != "lomax" & g != "log-logistic" &
      g != "lfr" & g != "chen") {
    stop("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.")
  }
  n<-length(mydata)
  sx<-sort(mydata)
  k<-1.25
  rankx<-(seq(1,n)-.5)/n
  hat<-mpsweibullg(mydata, g=g, location=location, method=method, sig.level=.05)$MPS
  quan<-qweibullg(rankx, g=g, hat, location=location, log.p=FALSE, lower.tail=TRUE)
  out<-plot(sort(quan), sx, main="Q-Q plot", xlab="Theoretical Quantiles", ylab=
              "Sample Quantiles", cex=k, cex.lab=k, cex.axis=k, col='black', lwd=2)
  lines(mydata, mydata, col='steelblue', cex=0.5, lwd=2)
}
