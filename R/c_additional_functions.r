#' JAGS output an internal function used within the tuts package
#'
#' \code{JAGS.objects} an internal function used within the tuts package.
#' @param JAGS.output An JAGS MCMC Object.
JAGS.objects=function(JAGS.output){
  chain=1
  Output.JAGS=list()
  JAGS.names=varnames(JAGS.output[[chain]])
  JAGS.objects=unique(gsub("\\[.*","",JAGS.names))
  for (i in 1:length(JAGS.objects)){
    ncol=sum((gsub("\\[.*","",JAGS.names))==JAGS.objects[i])
    objects=paste("JAGS.",JAGS.objects[i],"=matrix(data=NA,nrow=dim(JAGS.output[[1]])[1],ncol=ncol)",sep="")
    eval(parse(text=objects))
    object.names=paste("JAGS.",JAGS.objects[i],".names=rep(NA,ncol)",sep="")
    eval(parse(text=object.names))
  }
  for (Obj in JAGS.objects){
    j=1
    for (i in 1:length(JAGS.names)){
      if (gsub("\\[.*","",JAGS.names[i])==Obj){
        eval(parse(text=paste("JAGS.",Obj,".names[j]=JAGS.names[i]",sep="")))
        eval(parse(text=paste("JAGS.",Obj,"[,j]=JAGS.output[[1]][,JAGS.names[i]]",sep="")))
        j=j+1
      }
    }
    eval(parse(text=paste("colnames(JAGS.",Obj,")=JAGS.",Obj,".names",sep="")))
  }
  for (Obj in JAGS.objects){
    Output.JAGS[[(length(Output.JAGS)+1)]]=eval(parse(text=paste("JAGS.",Obj,sep="")))
  }
  # add optput
  names(Output.JAGS)=JAGS.objects
  return(Output.JAGS)
}

##############################################################################################
#' Wrapper of the models contained in the tuts package
#'
#' \code{tuwrap} tuwrap compares results obtained from
#'  fitting multiple models of time-uncertain time series.
#'
#' @param y A vector of observations.
#' @param ti.mu A vector of estimates/observed timings of observations.
#' @param ti.sd A vector of standard deviations of timings.
#' @param n.sim A number of simulations.
#' @param ... optional arguments: \cr
#' - CV: TRUE/FALSE cross-validation indicator, the default value is set to FALSE.
#' - n.chains: number of MCMC chains, the default number of chains is set to 2.\cr
#' - Thin: thinning factor, the default values is set to 4.\cr
#' - m: maximum number of significant frequencies in the data, the default value is set to 5. \cr
#' - polyorder: the polynomial regression component, the default odrer is set to 3. \cr
#' - freqs: set to a positive integer k returns a vector of k equally spaced frequencies in the Nyquist
#'   range. freqs can be provided as a vector of custom frequencies of interest. Set to 'internal'
#'   (the default value) generates a vector of equally spaced frequencies in the Nyquist range.
#' - n.cores: number of cores used in cross-validation. No value or 'MAX' applies all the available cores in computation.\cr
#'
#' @examples
#' # Note: Most of models included in tuts package are computationally intensive. In the example
#' # below I set parameters to meet CRAN’s testing requirement of maximum 5 sec per example.
#' # A more practical example would contain N=50 in the first line of the code and n.sim=10000.
#'
#' #1. Import or simulate the data (a simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=5,Harmonics=c(4,0,0), sin.ampl=c(10,0, 0), cos.ampl=c(0,0,0),
#'             trend=0,y.sd=2, ti.sd=0.2)
#' y=DATA$observed$y.obs
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(0.2, length(ti.mu))
#'
#' #2. Fit the models:
#' n.sim=100
#' WRAP=tuwrap(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim,CV=FALSE,n.cores=2)
#'
#' #3. Generate summary results:
#' summary(WRAP)
#'
#' # Note: Accessing individual summaries, diagnostics and plots is presented in manuals
#' #       of models contained in the wrapper.SOme examples:
#'
#' plot(WRAP$BFS,type='periodogram')
#' summary(WRAP$BFS,CI=0.99)
#'
#' @export

tuwrap=function(y,ti.mu,ti.sd,n.sim, ...){
  dots = list(...)

  if (is.null(dots$n.chains)) {n.chains=2}else{n.chains=dots$n.chains}
  if (is.null(dots$CV)) {CV=FALSE}else{CV=dots$CV}
  if (is.null(dots$freqs)) {freqs="internal"}else{freqs=dots$freqs}
  if (is.null(dots$polyorder)) {polyorder=3}else{polyorder=dots$polyorder}
  if (is.null(dots$n.cores)) {n.cores='MAX'}else{n.cores=dots$n.cores}

  if (is.null(dots$type)) {
   if(!sum(y==round(y) & y>=0)==length(y)){
     type="continuous"
    } else {
      RESP=readline(prompt="Are the observations discrete? (y/n):")
      if (RESP=='y' | RESP=='Y'){
        type="counting"
      } else{
        type="continuous"
      }
    }
  }

  output=list()
  if (type=="continuous"){
    BFS=tubfs(y=y,ti.mu=ti.mu,ti.sd=ti.sd,freqs=freqs,n.sim=n.sim,n.chains=n.chains, CV=CV,n.cores=n.cores)
    output$BFS=BFS

    TUAR1=tuar1(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim,n.chains=n.chains, CV=CV,n.cores=n.cores)
    output$TUAR1=TUAR1

    AR1REDF=tuar1redf(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim,n.chains=n.chains, CV=CV,n.cores=n.cores)
    output$AR1REDF=AR1REDF

    for (i in 1: polyorder){
      eval(parse(text = paste(
        "PN",i,"=tupolyn(y=y,ti.mu=ti.mu,ti.sd=ti.sd,polyorder=i,
            n.sim=n.sim,n.chains=n.chains, CV=CV,n.cores=n.cores)"
        ,sep='')))
      eval(parse(text = paste("output$PN",i,"=PN",i,sep='')))
    }


  } else if(type=="counting"){
    output$TUPOISBSF=tupoisbsf(y=y,ti.mu=ti.mu,ti.sd=ti.sd,freqs='internal',n.sim=n.sim,n.chains=n.chains, CV=CV,n.cores=n.cores)
    output$TUPOISPN3=tupoispn(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim,n.chains=n.chains, CV=CV,n.cores=n.cores)
  } else{
    stop("Type must be either a 'continuous' or 'counting' process.")
  }
  class(output)="tuts_wrap"
  return(output)
}

# =========
#' Model comaprison of multiple time-uncertain models based on DIC criterion
#'
#' \code{summary.tuts_wrap} function returns DIC criteria of multiple models contained in tuts_wrap objects.
#'
#' @param object A tuwrap object.
#' @param ... optional arguments, not in use in the current version on tuts.
#'
#' @examples
#' # Note: Most of models included in tuts package are computationally intensive. In the example
#' # below I set parameters to meet CRAN’s testing requirement of maximum 5 sec per example.
#' # A more practical example would contain N=50 in the first line of the code and n.sim=10000.
#'
#' #1. Import or simulate the data (a simulation is chosen for illustrative purposes):
#' DATA=simtuts(N=5,Harmonics=c(4,0,0), sin.ampl=c(10,0, 0), cos.ampl=c(0,0,0),
#'             trend=0,y.sd=2, ti.sd=0.2)
#' y=DATA$observed$y.obs
#' ti.mu=DATA$observed$ti.obs.tnorm
#' ti.sd= rep(0.2, length(ti.mu))
#'
#' #2. Fit the models:
#' n.sim=100
#' WRAP=tuwrap(y=y,ti.mu=ti.mu,ti.sd=ti.sd,n.sim=n.sim)
#'
#' #3. Generate summary results:
#' summary(WRAP)
#'
#' @export
#'
summary.tuts_wrap = function(object, ...) {
  dots = list(...)
  NAMES=names(object)
  DIC=array(NA,dim=c(length(NAMES),3))
  rownames(DIC)=NAMES; colnames(DIC)=c("Deviance","Penalty", "Penalized Dev")

  for(i in 1:length(NAMES)){
    Deviance=eval(parse(text = paste( "round(sum(object$",NAMES[i],"$DIC$deviance),2)",sep='')))
    Penalty=eval(parse(text = paste( "round(sum(object$",NAMES[i],"$DIC$penalty),2)",sep='')))
    DIC[i,1]=Deviance
    DIC[i,2]=Penalty
    DIC[i,3]=Deviance+Penalty
  }
  print(DIC)
}
