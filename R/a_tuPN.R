#' Time-uncertain polynomial regression
#'
#' \code{tupolyn} performs estimation of parameters of N-th order polynomial regression of time-uncertain time series.
#'
#' @param y A vector of observations.
#' @param ti.mu A vector of estimates of timing of observations.
#' @param ti.sd A vector of standard deviations of timing.
#' @param polyorder Order of the polynomial regression.
#' @param n.sim A number of simulations.
#' @param CV TRUE/FALSE cross-validation indicator.
#' @param ... optional arguments: \cr
#' - n.chains: number of MCMC chains, the default number of chains is set to 2.\cr
#' - Thin: thinning factor, the default values is set to 4.\cr
#' - n.cores: number of cores used in cross-validation. No value or 'MAX' applies all the available cores in computation.\cr
#'
#' @examples
#' # Note: Most of models included in tuts package are computationally intensive. In the example
#' # below I set parameters to meet CRAN's testing requirement of maximum 5 sec per example.
#' # A more practical example would contain N=50 in the first line of the code and n.sim=10000.
#'
#' #1. Import or simulate the data (a simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=10,Harmonics=c(4,0,0), sin.ampl=c(10,0, 0), cos.ampl=c(0,0,0),
#'             trend=0,y.sd=2, ti.sd=0.2)
#' y=DATA$observed$y.obs
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(0.2, length(ti.mu))
#'
#' #2. Fit the model:
#' polyorder=2
#' n.sim=1000
#' PN=tupolyn(y=y,ti.mu=ti.mu,ti.sd=ti.sd,polyorder=polyorder,n.sim=n.sim,CV=TRUE,n.cores=2)
#'
#' #3. Generate summary results (optional parameters are listed in brackets):
#' summary(PN)                               # Summary results (burn, CI).
#'
#' #4. Generate plots and diagnostics of the model (optional parameters are listed in brackets):
#' plot(PN,type='predTUTS',CI=0.95)          # One step out of salmple predictions (CI, burn).
#' plot(PN,type='cv',burn=0.3)               # 5 fold cross-validation (CI, burn).
#' plot(PN,type='GR',CI=0.95)                # Gelman-Rubin diagnostic (CI).
#' plot(PN,type='mcmc')                      # MCMC diagnostics.
#' plot(PN,type='volatility')                # Volatility realizaitons.
#' @export
tupolyn=function(y,ti.mu,ti.sd, n.sim, polyorder,CV=FALSE, ...){
# Data checking and basic operations ---------------------------------------------------------------------
  if (length(y)*4!=length(ti.mu)*2+length(ti.sd)*2){stop("Vectors y, ti.mu and ti.sd should be of equal lengths.")}
if(is.numeric(y)==FALSE ){stop("y must be a vector of rational numbers.")}
if(is.numeric(ti.mu)==FALSE | sum((ti.mu)<0)>0 ){
  stop("ti.mu must be a vector of positive rational numbers.")}
if(is.numeric(ti.sd)==FALSE | sum((ti.sd)<0)>0 ){
  stop("ti.sd must be a vector of positive rational numbers.")}
if (sum(is.na(c(y,ti.mu,ti.sd)))>0){stop("Remove NAs.")}
if (n.sim!=abs(round(n.sim))){stop("n.sim must be a positive integer.")}
  if (polyorder!=abs(round(polyorder))){stop("polyorder must be a positive integer.")}
if (!is.logical(CV)){stop("CV must be a logical value.")}
# Optional parameters
  dots = list(...)
  if(missing(...)){Thin=4; n.chains=2; polyorder=3 ; n.cores='MAX'}
  if(!is.numeric(dots$Thin)){
    Thin=4
  } else{
    Thin=round(abs(dots$Thin))
  }

  if(!is.numeric(dots$n.cores)){
    n.cores='MAX'
  } else{
    n.cores=dots$n.cores
  }

  if(!is.numeric(dots$n.chains)){
    n.chains=2
  } else{
    n.chains=round(abs(dots$n.chains))
  }


y=y[order(ti.mu,decreasing = FALSE)]; ti.sd=ti.sd[order(ti.mu,decreasing = FALSE)]
ti.mu=ti.mu[order(ti.mu,decreasing = FALSE)]

# JAGS model --------------------------------------------------------------------------------------------
modelstring= "model {
for(i in 1:n) {
  y[i]~dnorm(mu[i], precision )
  mu[i]<- const+inprod(ti.sim.MX[i,],alpha)
}
for (i in 1:n) {
  ti.sim.tmp[i]~dnorm(ti.mu[i],1/ti.sd[i])
}
ti.sim<-sort(ti.sim.tmp)

for (i in 1:n) {
  for(j in 1:(polyorder)){
    ti.sim.MX[i,j]<- ti.sim[i]^j
  }
}
for(i in 1:polyorder){
  alpha[i]~dnorm(0,0.01)
}
const~dnorm(0,0.01)
precision ~dgamma(0.01,0.01)
}"

# JAGS data ---------------------------------------------------------------------------------------------
data=list(y=y, ti.mu=ti.mu,ti.sd=ti.sd, n=length(ti.mu),polyorder=polyorder)

inits=parallel.seeds("base::BaseRNG", n.chains)
for(k in (1:n.chains)){
  inits[[k]]$ti.sim.tmp=ti.mu
}

model=jags.model(textConnection(modelstring), data=data, inits=inits,n.chains=n.chains)
update(model,n.iter=n.sim,thin=Thin)

output=coda.samples(model=model, variable.names=c("const","alpha","precision","ti.sim")
                      ,n.iter=n.sim, thin=Thin,n.chains=n.chains)
DIC = dic.samples(model=model,n.iter=n.sim,thin=Thin)

Sim.Objects=JAGS.objects(output)
Sim.Objects$JAGS=output
Sim.Objects$DIC=DIC
# Cross Validation -------------------------------------------------------------------------------------------------
if(CV==TRUE){
  print(noquote('Cross-validation of the model....'))
  folds = 5
  fold= sample(rep(1:folds,length=length(y)))
  for (i in 2:length(fold)){
    if (fold[i-1]==fold[i]){
      Sample=c(1:5)
      Sample=Sample[Sample !=fold[i]]
      fold[i]=sample(Sample,size=1)
    }
  }
  TI.SIM=apply(Sim.Objects$ti.sim,2,'quantile',0.5)
  Cores=min(c(parallel::detectCores()-1,folds))
  if(n.cores=="MAX"){
    Cores=Cores
  } else{
    if(Cores>n.cores) {
      Cores=n.cores} else {
        Cores=Cores
      }
  }
  cl = parallel::makeCluster(Cores)
  doParallel::registerDoParallel(cl)

  BSFCV=function(i,y,ti.mu,ti.sd,modelstring,n.sim,fold){
    Y=y; Y[fold==i]=NA;
    data=list(y=Y, ti.mu=TI.SIM,ti.sd=ti.sd, n=length(y),polyorder=polyorder)
    model=jags.model(textConnection(modelstring), data=data,n.chains=1)
    update(model,n.iter=n.sim,thin=Thin)
    output=coda.samples(model=model,variable.names=c("y"), n.iter=n.sim, thin=Thin)
    return(output)
  }

  CVRES=foreach(i=1:folds,.export=c('jags.model','coda.samples')) %dopar%
    BSFCV(i,fold=fold,y=y,ti.mu=TI.SIM,ti.sd=ti.sd, modelstring=modelstring, n.sim=n.sim)

  stopCluster(cl)

  if(is.null(dim(JAGS.objects(CVRES[[1]])$y[,fold==(1)]))){
    DIM1=length(JAGS.objects(CVRES[[1]])$y[,fold==(1)])
  }else{
    DIM1=  dim(JAGS.objects(CVRES[[1]])$y[,fold==(1)])[1]
  }

  pred_y = array(NA, dim=c(DIM1,length(ti.mu)))
  colnames(pred_y)= paste("y[",1:length(y),"]",sep="")



  for (i in 1:folds){
    pred_y[,fold==(i)] = JAGS.objects(CVRES[[i]])$y[,fold==(i)]
  }
  Sim.Objects$CVpred=pred_y
  }
Sim.Objects$y=y
Sim.Objects$ti.mu=ti.mu
Sim.Objects$polyorder=polyorder
class(Sim.Objects)='tuts_polyn'
return(Sim.Objects)
}

#' Prints summary tables of tuts_polyn objects
#'
#' \code{summary.tuts_polyn} prints summary tables of tuts_polyn objects.
#'
#' @param object A tuts_polyn object.
#'
#' @param ... A list of optional parameters: \cr
#'  - burn: burn-in parameter ranging from 0 to 0.7, the default value is 0.\cr
#'  - CI: confidence interval, the default value is set to 0.99.
#' @examples
#' # Note: Most of models included in tuts package are computationally intensive. In the example
#' # below I set parameters to meet CRAN’s testing requirement of maximum 5 sec per example.
#' # A more practical example would contain N=50 in the first line of the code and n.sim=10000.
#'
#' #1. Import or simulate the data (simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=10,Harmonics=c(4,0,0), sin.ampl=c(10,0, 0), cos.ampl=c(0,0,0),
#'             trend=0,y.sd=2, ti.sd=0.2)
#' y=DATA$observed$y.obs
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(0.2, length(ti.mu))
#'
#' #2. Fit the model:
#' polyorder=2
#' n.sim=1000
#' PN=tupolyn(y=y,ti.mu=ti.mu,ti.sd=ti.sd,polyorder=polyorder,n.sim=n.sim, CV=FALSE)
#'
#' #3. Generate summary results (optional parameters are listed in brackets):
#' summary(PN)                               # Summary results (burn, CI).
#' @export
summary.tuts_polyn = function(object, ...) {
  dots = list(...)
  if(missing(...)){burn=0; CI=0.99}

  if(!is.numeric(dots$CI)){
    CI=0.99
  } else{
    if(dots$CI<=0.5 | dots$CI> 1){stop('Credible interval is bounded between 0.5 and 1')}
    CI=dots$CI
  }

  if(!is.numeric(dots$burn)){
    burn=0
  } else{
    burn=dots$burn
    if(burn<0 | burn>0.7){stop('burn is bounded between 0 and 0.7')
    }
  }
  n.sim=dim(object$const)[1]
  if (burn==0){BURN=1}else{BURN=floor(burn*n.sim)}

  # ----------------------------------------------------------------------------
  cat('\n')
  cat('Regression Parameters and estimates of timing:\n')
  cat('----------------------------------------------\n')

  const=object$const[BURN:length(object$const)]
  const.lwr=quantile(const,(1-CI)/2)
  const.med=quantile(const,0.5)
  const.upr=quantile(const,1-(1-CI)/2)
  constName="const"

  lwr=med=upr=NAMES=NA
  polyorder=object$polyorder
  alpha=object$alpha[BURN:dim(object$alpha)[1],]
  alpha.lwr=apply(alpha,2,'quantile',(1-CI)/2)
  alpha.med=apply(alpha,2,'quantile',0.5)
  alpha.upr=apply(alpha,2,'quantile', 1-(1-CI)/2)
  alphaNames=colnames(alpha)

  precision=object$precision[BURN:length(object$precision)]
  precision.lwr=quantile(precision,(1-CI)/2)
  precision.med=quantile(precision,0.5)
  precision.upr=quantile(precision,1-(1-CI)/2)
  precisionName="precision"

  ti=object$ti.sim[BURN:dim(object$ti.sim)[1],]
  ti.lwr=apply(ti,2,'quantile',(1-CI)/2)
  ti.med=apply(ti,2,'quantile',0.5)
  ti.upr=apply(ti,2,'quantile',1-(1-CI)/2)
  tiNames=names(ti.med)

  LWR=c(const.lwr,alpha.lwr,precision.lwr,ti.lwr)
  MED=c(const.med,alpha.med,precision.med,ti.med)
  UPR=c(const.upr,alpha.upr,precision.upr,ti.upr)
  TABLE2=data.frame(LWR,MED,UPR)
  row.names(TABLE2)=c(constName,alphaNames,precisionName,tiNames)

  colnames(TABLE2)=c(paste(round((1-CI)/2,3)*100,"%",sep=""),'50%',paste(round(1-(1-CI)/2,3)*100,"%",sep=""))
  print(TABLE2)
  # ----------------------------------------------------------------------------
  cat('\n')
  cat('Deviance information criterion:\n')
  cat('-------------------------------\n')
  print(object$DIC)
  cat('-------------------------------\n')
}
####################################################################################################
#' Plots and visual diagnostics of tuts_polyn objects
#'
#' \code{plot.tuts_polyn} generates plots and visual diagnostics of tuts_polyn objects.
#'
#' @param x A tuts_polyn object.
#' @param type plot type with the following options:\cr
#'  - 'predTUTS' plots one step predictions of the model. \cr
#'  - 'GR' plots Gelman-Rubin diagnostics. \cr
#'  - 'cv' plots 5-fold cross validation. \cr
#'  - 'mcmc' plots diagnostics of MCMC/JAGS objects. \cr
#'  - 'volatility' plots volatility realizations. \cr
#' @param ... list of optional parameters:\cr
#'  - burn: burn-in parameter ranging from 0 to 0.7 with default value set to 0. \cr
#'  - CI: credible interval ranging from 0.3 to 1 with default value set to 0.95.
#' @examples
#' # Note: Most of models included in tuts package are computationally intensive. In the example
#' # below I set parameters to meet CRAN’s testing requirement of maximum 5 sec per example.
#' # A more practical example would contain N=50 in the first line of the code and n.sim=10000.
#'
#' #1. Import or simulate the data (simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=10,Harmonics=c(4,0,0), sin.ampl=c(10,0, 0), cos.ampl=c(0,0,0),
#'             trend=0,y.sd=2, ti.sd=0.2)
#' y=DATA$observed$y.obs
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(0.2, length(ti.mu))
#'
#' #2. Fit the model:
#' polyorder=2
#' n.sim=1000
#' PN=tupolyn(y=y,ti.mu=ti.mu,ti.sd=ti.sd,polyorder=polyorder,n.sim=n.sim,CV=TRUE,n.cores=2)
#'
#' #3. Generate plots and diagnostics of the model (optional parameters are listed in brackets):
#' plot(PN,type='predTUTS',CI=0.95)          # One step out of salmple predictions (CI, burn).
#' plot(PN,type='cv',burn=0.3)               # 5 fold cross-validation (CI, burn).
#' plot(PN,type='GR',CI=0.95)                # Gelman-Rubin diagnostic (CI).
#' plot(PN,type='mcmc')                      # MCMC diagnostics.
#' plot(PN,type='volatility')                # Volatility realizaitons.
#' @export
plot.tuts_polyn = function(x, type, ...) {
  if (sum(type==c('predTUTS','GR','cv','mcmc','volatility'))==0){
    stop('type should be set as either predTUTS, GR, cv, mcmc or volatility')
  }

  dots = list(...)
  if(missing(...)){burn=0; CI=0.99}

  if(!is.numeric(dots$CI)){
    CI=0.99
  } else{
    if(dots$CI<=0.5 | dots$CI> 1){stop('Credible interval is bounded between 0.5 and 1')}
    CI=dots$CI
  }

  if(!is.numeric(dots$burn)){
    burn=0
  } else{
    burn=dots$burn
    if(burn<0 | burn>0.7){stop('burn is bounded between 0 and 0.7')
    }
  }
  n.sim=dim(x$const)[1]
  if (burn==0){BURN=1}else{BURN=floor(burn*n.sim)}

  graphics::par(mfrow=c(1,1))
  # ----------------------------------------------------------------------------
  if(type=='cv') {
    if (sum(names(x)=="CVpred")<1){stop("Object does not contain cross validation")}
    PRED=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',0.5)
    MAIN="Cross-Validation: One step out of sample predictions"

    graphics::plot(x=x$y,y=PRED,xlab="Actual",ylab="Predicted", main=MAIN, pch=18)
    graphics::abline(0,1,col='blue')

    RSQ=cor(x$y,PRED)^2* 100

    LAB = bquote(italic(R)^2 == .(paste(format(RSQ, digits = 0),"%",sep="")))
    graphics::text(x=(min(x$y)+0.9*(max(x$y)-min(x$y))),y=(min(PRED)+0.1*(max(PRED)-min(PRED))),LAB)
  }
  # ----------------------------------------------------------------------------
  if(type=='predTUTS') {
    if (sum(names(x)=="CVpred")<1){stop("Object does not contain cross validation")}
    PRED.LWR=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',(1-CI)/2)
    PRED.MED=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',0.5)
    PRED.UPR=apply(x$CVpred[BURN:dim(x$CVpred)[1],],2,'quantile',1-(1-CI)/2)

    ti.sim=apply(x$ti.sim[BURN:dim(x$ti.sim)[1],],2,'quantile',0.5)

    MAIN=paste("One step out of sample predictions at CI= ", CI*100,"%",sep='')
    graphics::plot(y=x$y,x=x$ti.mu,type='l',main=MAIN,ylab="Observations",xlab='time',lwd=2,
         ylim=c(min(x$CVpred),1.2*max(x$CVpred)), xlim=c(min(x$ti.mu,ti.sim),max(x$ti.mu,ti.sim)))
    graphics::lines(y=PRED.LWR,x=ti.sim,type='l',col='blue',lwd=1,lty=2)
    graphics::lines(y=PRED.MED,x=ti.sim,type='l',col='blue',lwd=1,lty=1)
    graphics::lines(y=PRED.UPR,x=ti.sim,type='l',col='blue',lwd=1,lty=2)

    graphics::legend("topright",legend = c("Observed","Upper CI","Median","Lower CI"),
           col=c("black","blue","blue","blue"),lwd=c(2,1,1,1),lty=c(1,2,1,2))
  }
  # ----------------------------------------------------------------------------
  if(type=='GR') {
    options(warn=-1)
    if(burn>0){ABURNIN=TRUE} else{ABURNIN=FALSE}

    GELMAN=gelman.diag(x$JAGS,confidence =CI,autoburnin=ABURNIN)
    TIOBJ=rownames(GELMAN$psrf)%in% c(paste("ti.sim[",1:length(x$ti.mu),"]",sep=""))

    REG=GELMAN$psrf[!TIOBJ,]
    labels =rownames(REG)
    TI=GELMAN$psrf[TIOBJ,]

    graphics::par(mfrow=c(2,1))
    graphics::plot(REG[,1],ylim=c(0,max(REG)+1.5),xaxt='n',ylab="Factor",
         xlab='Parameters of polynomial regression and precision of the model',
         main=paste("Gelman-Rubin diagnostics: Potential scale reduction factors \n with the upper confidence bounds at ",
                    CI*100,"%",sep="")

    )
    graphics::text(seq(1,length(labels),by=1), graphics::par("usr")[3]-max(REG)*0.2, srt = 00, adj= 0.5, xpd = TRUE, labels =labels, cex=0.65)
    graphics::arrows(x0=1:dim(REG)[1],y0=REG[,1],x1=1:dim(REG)[1],y1=REG[,2],code=3,length=0.04,angle=90,col='darkgray')
    graphics::abline(h=1);
    graphics::legend("topright", c("Estimate"),lty=c(NA),pch=c(1),lwd=c(1), col=c("black"),border="white")

    graphics::plot(TI[,1],ylim=c(0,(max(TI)*1.8)),ylab="Factor", xlab='Simulated timings of observations')
    graphics::arrows(x0=1:dim(TI)[1],y0=TI[,1],x1=1:dim(TI)[1],y1=TI[,2],code=3,length=0.04,angle=90,col='darkgray')
    graphics::abline(h=1);
    graphics::legend("topright", c("Estimate"),lty=c(NA),pch=c(1),lwd=c(1), col=c("black"),border="white")

    graphics::par(mfrow=c(1,1))
    options(warn=0)
  }
  if(type=='volatility') {

    graphics::plot(sqrt(sqrt(1/sqrt(x$precision[BURN:dim(x$precision)[1],]))),type='l', xlab='Sim ID',ylab='Std Deviation',main='Standard Deviaiton')
  }
  if(type=='mcmc') {
    options(warn=-1)
    mcmcplot(x$JAGS)
    options(warn=0)
  }
}

