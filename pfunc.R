### R function to fit SOR with log-bilinear form of odds ratio functions. 
fitsor<-function(y, x, weight=matrix(1, length(y), length(unique(y))), strata=rep(1, length(y)), Z=NULL)
{

# Arguments:
#      y --- the vector of outcome
#      x --- predictor matrix
#      weight --- optional. the sampling weights, supplied if known. Otherwise, left as is. 
#      strata --- optional. Strata employed in sampling.

#  Start with some preliminary calculations and error checks

	call <- match.call()
        if (any(is.na(y)) | any(is.na(x)) ) stop("Data should contain no missing data")


      ##  standardize the variables for stable computation.
      x=as.matrix(x)
      p=ncol(x)
      Xstd= cbind(y,x)
      nstrata=length(unique(strata))
      ## write out the dataset to hard disk and run Fortran program.
     write.table(t(c(nrow(Xstd),p, nstrata)), file="dataforor.dat", col.names=F, row.names=F)
     if (nstrata==1) write.table(Xstd, file="dataforor.dat", row.names=FALSE, col.names=FALSE, append=TRUE) else
       write.table(cbind(Xstd, strata), file="dataforor.dat", row.names=FALSE, col.names=FALSE, append=TRUE)
     if (nstrata==1) write.table(weight, file="dataforor.dat", row.names=FALSE, col.names=FALSE, append=TRUE)
     ##system("Outsample.exe", invisible=FALSE)
     system("cmd.exe", input="Outsample")

     ## then read the outpout
     gamma=read.table("gamma_file.out")
     yunique=0 
     gammavar=as.matrix(read.table("variance.out"))
     negllk=as.matrix(read.table("llk_file.out"))
     
	# save, and return
	sorobj <- list(call = call, lambdagamma=gamma, yunique=yunique, gammase=sqrt(diag(gammavar)), gammavar=gammavar, negllk=negllk
                         )
	class(sorobj) <- "sor"
	return(sorobj)
} # END OF FITSOR FUNCTION



####################################################################################
###### parametric OR model ########################################################
#######################################################################################

fitpor<-function(y, x, strata=rep(1, length(y)), dist="gaussian", par0)
{
 
 res <- optim(par0, llk.por, y=y,x=x, method="L-BFGS-B", hessian=T,strata)
 res

}


## Assuming parametric baseline function
llk.por<-function(para, y,x, dist="gaussian", strata=NULL){
  x=as.matrix(x)
  if (dist=="gaussian") {
   if (is.null(strata)){
    sigma=exp(para[length(para)])
    para[2:ncol(x)]=para[2:ncol(x)]*sigma^2
    llk= dnorm(y, mean= x%*%para[-length(para)], sd=sigma, log=T)} 
   else {
     nstrata=length(unique(strata))
     beta0=para[1:nstrata]
    gamma=para[(nstrata+1):(nstrata+ncol(x))]
    sigma=exp(para[(nstrata+ncol(x)+1):(nstrata+ncol(x)+nstrata)])
    mu=beta0[strata]+c(x%*%gamma)*sigma[strata]^2
    llk= dnorm(y, mean=mu, sd=sigma[strata], log=T)
        }
                     } else 
  if (dist=="poisson") {
     if (is.null(strata)) {
       llk= dpois(y, lambda=exp(x%*%para), log=T) 
       } else 
      {
      nstrata=length(unique(strata))
     beta0=para[1:nstrata]
    gamma=para[(nstrata+1):(nstrata+ncol(x))]
    mu=exp(beta0[strata]+c(x%*%gamma))
    llk= dpois(y, lambda=mu, log=T)
       }
                       } else
  stop("Incorrect Distribution")
  sum(-llk)
}


###################################################################################
##### Hausman & wise method #######################################################
###################################################################################
fun.hw<-function(y, x, stratapoints=NULL, weights=NULL, dist="gaussian", par0)
{
 
 res <- optim(par0, llk.hw, y=y,x=x, method="L-BFGS-B", hessian=T,stratapoints=stratapoints, weights=weights,dist=dist)
 res

}

llk.hw<-function(para, y,x, dist="gaussian", stratapoints=NULL, weights=NULL){
  x=as.matrix(x)

  if (dist=="gaussian") {
    beta = para[1:ncol(x)]
    mu = as.matrix(x)%*%beta
    sigma=exp(para[ncol(x)+1])
    if (is.null(weights)) weights=c(1, exp(para[(ncol(x)+2):length(para)]))
    nc = numeric(length(y))
    wt= weights
    for (i in 1:(length(wt)-1)) wt[i]=wt[i]-wt[i+1]
    for (i in 1:length(nc))  {
       nc[i]=-log(sum(wt[1:length(stratapoints)]*pnorm(stratapoints, mean=mu[i], sd=sigma)) + weights[length(weights)])
       if (y[i] > max(stratapoints))  nc[i] = nc[i]+ log(weights[length(weights)]) else
        {
          for (j in 1:length(stratapoints)) {
           if (y[i] < stratapoints[j]) {nc[i]=nc[i]+log(weights[j]);break} 
                                         }              
        }
                            }
       llk= dnorm(y, mean=mu, sd=sigma, log=T) + nc
                     } else
   if (dist=="poisson")  {
      beta = para[1:ncol(x)]
    mu = exp(as.matrix(x)%*%beta)
    if (is.null(weights)) weights=c(1, exp(para[(ncol(x)+1):length(para)]))
    nc = numeric(length(y))
    wt= weights
    for (i in 1:(length(wt)-1)) wt[i]=wt[i]-wt[i+1]
    for (i in 1:length(nc))  {
       ##print(wt); print(ppois(stratapoints, lambda=mu[i]))
       nc[i]=-log(sum(wt[1:length(stratapoints)]*ppois(stratapoints, lambda=mu[i])) + weights[length(weights)])
       if (y[i] > max(stratapoints))  nc[i] = nc[i]+ log(weights[length(weights)]) else
        {
          for (j in 1:length(stratapoints)) {
           if (y[i] <= stratapoints[j]) {nc[i]=nc[i]+log(weights[j]);break} 
                                         }              
        }
                            }
       llk= dpois(y, lambda=mu, log=T) + nc
        } else
    stop("Incorrect Distribution")
                     
  sum(-llk)
}


###################################################################################
##### Cosslet Method #######################################################
###################################################################################
fun.coss<-function(y, x, stratapoints=NULL, weights=NULL, dist="gaussian", par0)
{
 
 res <- optim(par0, llk.coss, y=y,x=x, method="L-BFGS-B", hessian=T,stratapoints=stratapoints, dist=dist, weights=weights)
 res

}

llk.coss<-function(para, y,x, dist="gaussian", stratapoints=NULL, weights=NULL){
  x=as.matrix(x)

  if (dist=="gaussian") {
    beta = para[1:ncol(x)]
    mu = as.matrix(x)%*%beta
    sigma=exp(para[ncol(x)+1])
    if (is.null(weights)) stop("Must Specify Weights") ##weights=c(1, exp(para[(ncol(x)+2):length(para)]))
    nc = numeric(length(y))
    for (i in 1:length(nc))  {
         wt= weights[i,]
        for (j in 1:(length(wt)-1)) wt[j]=wt[j]-wt[j+1]
       nc[i]=-log(sum(wt[1:length(stratapoints)]*pnorm(stratapoints, mean=mu[i], sd=sigma)) + weights[i,length(wt)])
       if (y[i] > max(stratapoints))  nc[i] = nc[i]+ log(weights[i,length(wt)]) else
        {
          for (j in 1:length(stratapoints)) {
           if (y[i] < stratapoints[j]) {nc[i]=nc[i]+log(weights[i,j]);break} 
                                         }              
        }
                            }
       llk= dnorm(y, mean=mu, sd=sigma, log=T) + nc
                     } else
   if (dist=="poisson") {
     beta = para[1:ncol(x)]
    mu = exp(as.matrix(x)%*%beta)
    if (is.null(weights)) stop("Must Specify Weights") ##weights=c(1, exp(para[(ncol(x)+2):length(para)]))
    nc = numeric(length(y))
    for (i in 1:length(nc))  {
         wt= weights[i,]
        for (j in 1:(length(wt)-1)) wt[j]=wt[j]-wt[j+1]
       nc[i]=-log(sum(wt[1:length(stratapoints)]*ppois(stratapoints, lambda=mu[i])) + weights[i,length(wt)])
       if (y[i] > max(stratapoints))  nc[i] = nc[i]+ log(weights[i,length(wt)]) else
        {
          for (j in 1:length(stratapoints)) {
           if (y[i] <= stratapoints[j]) {nc[i]=nc[i]+log(weights[i,j]);break} 
                                         }              
        } 
                            }
       llk= dpois(y, lambda=mu, log=T) + nc
        } else
    stop("Incorrect distribution")
                     
  sum(-llk)
}


## Random number generation from density function (bx)/(e^(bx)-1) e^(bxy))

rsor<-function(n, b,x){
  ##u=runif(n)
  ##res=log(u*exp(b*x-1)+1)/(b*x)
  u=seq(0,1, by=0.001)
  res=numeric(n)
  for (i in 1:n) {
      p=exp(b*x[i]*u)/sum(exp(b*x[i]*u))
      res[i]=u[rmultinom(1,1, p)[,1]==1]
             }
  res
}




pred.sor <- function(uy, Xgrid,lambda, gamma) {
  Xgrid<-as.matrix(Xgrid)
  emax=as.matrix(Xgrid)%*%(as.matrix(gamma))%*%t(as.matrix(uy)) 
   emax=emax+ matrix(rep(as.matrix(c(lambda, 0)),each=nrow(Xgrid)),nrow=nrow(Xgrid))
   emax=exp(emax)
   emax=t(apply(emax, 1, function(i) i/sum(i)))
  ## apply(emax*matrix(rep(uy, each=length(Xgrid)), nrow =length(Xgrid)), 1, sum)
}


fun.bootgamma <- function (popdata, xc, weight, dist="gaussian", inter=FALSE) {
    Y10= quantile(popdata$Y, 0.1); Y90=quantile(popdata$Y, 0.9)
    N=nrow(popdata)
    if (inter) {
      gammaboot=matrix(0,20,ncol(popdata)) 
      cossboot=matrix(0,20,ncol(popdata)) 
     } else 
    {
      gammaboot=matrix(0,20,ncol(popdata)-1) 
      cossboot=matrix(0,20,ncol(popdata)-1) 
    }
     for (j in 1:20) {
     print(j)
      simdata = rbind(popdata[sample((1:N)[popdata$Y<=Y10 & popdata$X1<xc], size=n*0.2),], 
                          popdata[sample((1:N)[popdata$Y> Y10 & popdata$Y<= Y90 & popdata$X1<xc], size=n*0.1),], 
            popdata[sample((1:N)[ popdata$Y>Y90 & popdata$X1<xc ], size=n*0.2),],
         popdata[sample((1:N)[popdata$Y<=Y10 & popdata$X1>=xc], size=n*0.2),], 
       popdata[sample((1:N)[popdata$Y> Y10 & popdata$Y<= Y90 & popdata$X1>=xc], size=n*0.1),], 
           popdata[sample((1:N)[ popdata$Y>Y90 & popdata$X1>=xc], size=n*0.2),])
     uy= sort(unique(round(simdata$Y,1)))
     wt = matrix(1, nrow=nrow(simdata), ncol=length(uy))
     wt[simdata$X1< xc, uy <= Y10 ] = weight[1]
     wt[simdata$X1 < xc, uy> Y10 & uy<=Y90 ] = weight[2]
     wt[simdata$X1 < xc, uy > Y90  ] = weight[3]
     wt[simdata$X1 >= xc, uy <= Y10  ] = weight[4]
     wt[simdata$X1 >= xc, uy> Y10 & uy<=Y90  ] = weight[5]
      wt[simdata$X1 >= xc, uy > Y90 ] = weight[6]
   cosswt=matrix(1, nrow=nrow(simdata), ncol=3)
   for (k in 1:nrow(simdata)) {
   if (simdata$X1[k] <xc) cosswt[k, ]=weight[1:3] else 
    cosswt[k,]=weight[4:6]
          }
   ## with known weights
    ## temp=fitsor(simdata, digits=c(2,2), predictorMatrix=1-diag(2), visitSequence=2, weight=wt)
   if (inter) {
       temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2, simdata$X1*simdata$X2), weight=wt)
    tempcoss=fun.coss(y=simdata$Y, x=cbind(1,simdata$X1, simdata$X2, simdata$X1*simdata$X2), stratapoints=c(Y10, Y90), weights=cosswt, dist=dist, par0=c(0,0.1,0.1,0.1, 0))
    nlg=length(temp$lambdagamma)
     gammaboot[j,]=as.matrix(temp$lambdagamma)[1,((nlg-2):nlg)]
     cossboot[j,]=tempcoss$par[2:4]
         } else 
     {
       temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2), weight=wt)
    tempcoss=fun.coss(y=simdata$Y, x=cbind(1,simdata$X1, simdata$X2), stratapoints=c(Y10, Y90), weights=cosswt, dist=dist, par0=c(0,0.1,0.1,0))
    nlg=length(temp$lambdagamma)
     gammaboot[j,]=as.matrix(temp$lambdagamma)[1,((nlg-1):nlg)]
     cossboot[j,]=tempcoss$par[2:3]
     }
   }
  return(list(gammaboot=gammaboot, cossboot=cossboot))
}


### Truncated negative binomial (random sampling and onsite sampling) regression using R
fun.dnb<-function (y, lalpha, lmu) {
 ## parametrization in Englin paper. 
 ## lalpha=lalpha-lmu 
 theta=1/exp(lalpha); if (theta==0) theta=exp(-700)
 ##if (is.na(gamma(theta))) stop("gamma(theta) is NA")
 return(lgamma(y+theta) - lgamma(y+1)-lgamma(theta)+(lmu+lalpha)*y+ log(1+exp(lmu+lalpha))*(-y-theta))
}

fun.dTnb<-function (y, lalpha, lmu)  return(fun.dnb(y,lalpha,lmu)-log(1-(1+exp(lalpha+lmu))^(-1/exp(lalpha))))
fun.dOSnb<-function (y, lalpha, lmu)  return(fun.dnb(y,lalpha,lmu)+log(y)-lmu)


fun.nllknb<-function(par, y,x, dist="nb"){
  beta=par[-length(par)]; lalpha=par[length(par)]
  lmu=c(x%*%beta)
  res<-numeric(length(y)); 
  for (i in 1:length(res)) {
    if (dist=="nb") { res[i]= fun.dnb(y[i],lalpha,lmu[i])} else
    if (dist=="Tnb") {  res[i]= fun.dTnb(y[i],lalpha,lmu[i])} else
    if (dist=="OSnb") {  res[i]= fun.dOSnb(y[i],lalpha,lmu[i])}
    }
  return(-sum(res)) 
}


fun.fitnb<- function(ymodel, data, par0, dist="nb"){
    ymodel.data= model.frame(ymodel,data=data,na.action=NULL)
  x<- as.matrix(model.matrix(ymodel, data=ymodel.data))
  y<- as.vector(model.extract(ymodel.data, "response"))
  if (dist=="nb") res <- optim(par0,fun.nllknb, method="BFGS",y=y,x=x, hessian=T, control=list(trace=T)) else 
  if (dist=="Tnb") res <- optim(par0,fun.nllknb, method="BFGS",y=y,x=x, dist="Tnb", hessian=T, control=list(trace=T,factr=1e-20)) else
  if (dist=="OSnb") res <- optim(par0,fun.nllknb, method="BFGS",y=y,x=x, dist="OSnb", hessian=T, control=list(trace=T,factr=1e-20))
  return(res)
}


## SOR with Log-linear-softplus form OR regression in R 

Rfitsor<-function(ymodel, data, par0,  hessian=T) {
   ymodel.data= model.frame(ymodel,data=data,na.action=NULL)
  x<- as.matrix(model.matrix(ymodel, data=ymodel.data))
  x<-x[,-1] ## remove intercept
  y<- as.vector(model.extract(ymodel.data, "response"))
  uy=sort(unique(y))
  res <- optim(par0,fun.nllsor, method="BFGS",data=list(y=y,x=x,uy=uy), hessian=hessian, control=list(trace=T, maxit = 30000, factr=1e-20))
  res
}
fun.nllsor <- function(par, data){
  y <- data$y
  x <- as.matrix(data$x)
  uy <- data$uy
  n=length(y)
  lambda <- par[1:(length(par)-ncol(x)-1)]
  gamma <- par[(length(par)-ncol(x)):(length(par)-1)]
  alpha <-exp(par[length(par)])
  res <- llk.sor(y=y, uy=uy,x=x,  gamma=gamma, lambda=lambda, alpha=alpha) 
  -res
}

llk.sor <- function(y, uy, x, gamma, lambda, alpha) {
  emax=as.matrix(x)%*%(as.matrix(gamma))
  emax=emax-log(1+alpha*exp(emax))+log(1+alpha); 
  emax=emax%*%t(as.matrix(uy)) ##temp$lambdagamma[1,length(temp$lambdagamma)]
   emax=emax+ matrix(rep(as.matrix(c(lambda,0)),each=length(y)),nrow=length(y))  ## as.vector(temp$lambdagamma[1,-length(temp$lambdagamma)])
   emax=exp(emax)
   emax=t(apply(emax, 1, function(i) i/sum(i)))
  res=0
  for (i in 1:length(y)) res=res+log(emax[i,y[i]==uy])
  res
}


pred.sorNB <- function(uy, Xgrid,lambda, gamma,alpha) {
  Xgrid<-as.matrix(Xgrid)
  emax=as.matrix(Xgrid)%*%(as.matrix(gamma))
  emax=emax-log(1+alpha*exp(emax))+log(1+alpha)
  emax=emax%*%t(as.matrix(uy)) 
   emax=emax+ matrix(rep(as.matrix(c(lambda, 0)),each=nrow(Xgrid)),nrow=nrow(Xgrid))
   emax=exp(emax)
   emax=t(apply(emax, 1, function(i) i/sum(i)))
  ## apply(emax*matrix(rep(uy, each=length(Xgrid)), nrow =length(Xgrid)), 1, sum)
}


