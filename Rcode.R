setwd("Change to Directory where all files are located")
source("pfunc.R")

#####################################################################################################
#### Part 1: Code to compare HW, COSSLET and SOR for handling endogenous selective sampling
#####################################################################################################

library(truncnorm)
library(truncreg)
library(extraDistr)
library(VGAM)

#########################
## outcome is normal 
#############################
N=1000000
X2 <- rnorm(N) 
X1 <- rbinom(N, size=1, p=0.15) 
Y <- rnorm(N, mean=0.1*X1+0.1*X2) 
popdata <-data.frame(X1=X1,X2=X2, Y=Y)
N <- nrow(popdata)
lm(Y~X1+X2, data=popdata)

## Endogenous outcome-dependent selective sampling
n=500
nsim=500
srs.lm = srs.por = srs.hw=srs.hwfi=matrix(0, nsim,4)
srs.sor = matrix(0, nsim,4)
ods.lm= ods.por= ods.hw=ods.hwfi=matrix(0, nsim,4)
ods.sor=ods.sorfi= matrix(0, nsim, 4)
odsyx.lm= odsyx.por=odsyx.coss= matrix(0, nsim,4)
odsyx.sor= odsyx.sorfi=matrix(0, nsim, 4)
d.srslm=d.srspor=d.srshw=d.srshwfi=d.srssor=d.odslm=d.odspor=d.odshw=d.odshwfi=d.odssor=d.odssorfi=numeric(nsim)
d.odsyxlm=d.odsyxpor=d.odsyxsor=d.odsyxsorfi=d.odsyxcoss=numeric(nsim)
Y10= quantile(Y, 0.1); Y90=quantile(Y, 0.9)

## Simple random sampling
for (i in 1:nsim) {
  print(i)
   ## simple random sample with linear regression 
   simdata.srs= popdata[sample(1:N, size=n),]
   sim.reg<-lm(Y~X1+X2, data=simdata.srs)
   srs.lm[i,1:2] <- sim.reg$coef[2:3] 
   srs.lm[i,3:4] <- summary(sim.reg)$coef[2:3,2]
   d.srslm[i]= sqrt(det(vcov(sim.reg)[2:3,2:3]))



   ## simple random ssample with POR   
   fit.por <- optim(c(0,1,1, 0), llk.por, y=simdata.srs$Y, x=cbind(1, simdata.srs$X1, simdata.srs$X2), method="L-BFGS-B", hessian=T)
   srs.por[i,1:2] <- fit.por$par[2:3]
   srs.por[i,3:4] <- sqrt(diag(solve(fit.por$hessian)))[2:3]
   d.srspor[i]= sqrt(det(solve(fit.por$hessian)[2:3,2:3]))


   ## Simple random sampling with HW
   fithw=fun.hw(y=simdata.srs$Y, x=cbind(1,simdata.srs$X1, simdata.srs$X2), stratapoints=c(Y10, Y90), weights=c(1,1,1), dist="gaussian", par0=c(0,0.1,0.1,0))
   srs.hwfi[i,]=c(fithw$par[2:3], sqrt(diag(solve(fithw$hessian)))[2:3])
   d.srshwfi[i]= sqrt(det(solve(fithw$hessian)[2:3,2:3]))

   ## SRS with HW without weights
   fithw=fun.hw(y=simdata.srs$Y, x=cbind(1,simdata.srs$X1, simdata.srs$X2), stratapoints=c(Y10, Y90),  dist="gaussian", par0=c(0,0.1,0.1,0, 0,0))
   srs.hw[i,]= c(fithw$par[2:3], sqrt(diag(solve(fithw$hessian)))[2:3])
   d.srshw[i]= sqrt(det(solve(fithw$hessian)[2:3,2:3]))

   ## simple random sample with sor 
   temp.srs=fitsor(y=round(simdata.srs$Y,1), x=cbind(simdata.srs$X1, simdata.srs$X2))
   srs.sor[i,1:2]= as.matrix(temp.srs$lambdagamma)[1,(length(temp.srs$lambdagamma)-1):length(temp.srs$lambdagamma)]
   srs.sor[i,3:4]= c(temp.srs$gammase) 
   d.srssor[i]= sqrt(det(temp.srs$gammavar))
}


   #####################  ODS on y only ########################
for (i in 1:nsim) {
  print(i)
   ## ODS with OLS
   simdata = rbind(popdata[sample((1:N)[Y<=Y10], size=n*0.4),], popdata[sample((1:N)[Y> Y10 & Y<= Y90 ], size=n*0.2),], popdata[sample((1:N)[ Y>Y90 ], size=n*0.4),])
   sim.reg<-lm(Y~X1+X2, data=simdata)
   ods.lm[i,1:2] <- sim.reg$coef[2:3] 
   ods.lm[i,3:4] <- summary(sim.reg)$coef[2:3,2]
   d.odslm[i]= sqrt(det(vcov(sim.reg)[2:3,2:3])) 
  

   ## ods with POR   
   fit.por <- optim(c(0,1,1, 0), llk.por, y=simdata$Y, x=cbind(1, simdata$X1, simdata$X2), method="L-BFGS-B", hessian=T)
   ods.por[i,1:2] <- fit.por$par[2:3]
   ods.por[i,3:4] <- sqrt(diag(solve(fit.por$hessian)))[2:3]
    d.odspor[i]= sqrt(det(solve(fit.por$hessian)[2:3,2:3]))

      ## ODS  with HW+ weights
    fithw=fun.hw(y=simdata$Y, x=cbind(1,simdata$X1, simdata$X2), stratapoints=c(Y10, Y90), weights=c(4,0.25,4), dist="gaussian", par0=c(0,0.1,0.1,0))
    ods.hwfi[i,]=c(fithw$par[2:3], sqrt(diag(solve(fithw$hessian)))[2:3])
    d.odshwfi[i]= sqrt(det(solve(fithw$hessian)[2:3,2:3]))

   ## ODS with HW without weights
   fithw=fun.hw(y=simdata$Y, x=cbind(1,simdata$X1, simdata$X2), stratapoints=c(Y10, Y90),  dist="gaussian", par0=c(0,0.1,0.1,0, 0,0))
   ods.hw[i,]= c(fithw$par[2:3], sqrt(diag(solve(fithw$hessian)))[2:3])
   d.odshw[i]= sqrt(det(solve(fithw$hessian)[2:3,2:3]))

    ## ODS with sor 
   temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2))
   ods.sor[i,1:2]= as.matrix(temp$lambdagamma)[1,(length(temp$lambdagamma)-1):length(temp$lambdagamma)]
   ods.sor[i,3:4]= c(temp$gammase) 
   d.odssor[i]= sqrt(det(temp$gammavar))

   ## ODS with SOR and weights
    uy= sort(unique(round(simdata$Y,1)))
   wt <- matrix(1, nrow(simdata), length(uy))
   wt[, uy <= Y10] = 0.4/sum(Y<=Y10)*N
   wt[, uy> Y10 & uy<=Y90] = 0.2/sum(Y>Y10 & Y<=Y90)*N
   wt[, uy > Y90] = 0.4/sum(Y > Y90)*N
   temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2), weight=wt)
   ods.sorfi[i,1:2]= as.matrix(temp$lambdagamma)[1,(length(temp$lambdagamma)-1):length(temp$lambdagamma)]
   ods.sorfi[i,3:4]= c(temp$gammase)
   d.odssorfi[i]= sqrt(det(temp$gammavar))
}


   #####################  ODS on both y and x ########################
for (i in 1:nsim) {
   print(i)

     xc=0.5 
   if (i==1) {
   weight = c(0.2/sum(Y<=Y10 & X1<xc), 0.1/sum(Y> Y10 & Y<= Y90 & X1<xc), 0.2/sum( Y>Y90 & X1<xc ),
          0.2/sum(Y<=Y10 & X1>=xc), 0.1/sum(Y> Y10 & Y<= Y90 & X1>=xc), 0.2/sum(Y>Y90 & X1>=xc))*N
   weight= weight/sum(weight)
             }
   simdata = rbind(popdata[sample((1:N)[Y<=Y10 & X1<xc], size=n*0.2),], popdata[sample((1:N)[Y> Y10 & Y<= Y90 & X1<xc], size=n*0.1),], 
                popdata[sample((1:N)[ Y>Y90 & X1<xc ], size=n*0.2),],popdata[sample((1:N)[Y<=Y10 & X1>=xc], size=n*0.2),], 
           popdata[sample((1:N)[Y> Y10 & Y<= Y90 & X1>=xc], size=n*0.1),], popdata[sample((1:N)[ Y>Y90 & X1>=xc], size=n*0.2),])
   
   ## ODS with regular analysis
   sim.reg<-lm(Y~X1+X2, data=simdata)
   odsyx.lm[i,1:2] <- sim.reg$coef[2:3] ##/(summary(sim.reg)$sigma^2)
   odsyx.lm[i,3:4] <- summary(sim.reg)$coef[2:3,2]
   d.odsyxlm[i]= sqrt(det(vcov(sim.reg)[2:3,2:3])) ##sqrt(det(sim.reg$vcov[2:3,2:3])) ## 

   ## ODS  with COSS + weights   
  cosswt=matrix(1, nrow=nrow(simdata), ncol=3)
  for (j in 1:nrow(simdata)) {
   if (simdata$X1[j] <xc) cosswt[j, ]=weight[1:3] else 
    cosswt[j,]=weight[4:6]
          }
    fitcoss=fun.coss(y=simdata$Y, x=cbind(1,simdata$X1, simdata$X2), stratapoints=c(Y10, Y90), weights=cosswt, dist="gaussian", par0=c(0,0.1,0.1,0))
   odsyx.coss[i,]=c(fitcoss$par[2:3], sqrt(diag(solve(fitcoss$hessian)))[2:3])
    d.odsyxcoss[i]= sqrt(det(solve(fitcoss$hessian)[2:3,2:3]))


  ## ods with SOR +weights
   uy= sort(unique(round(simdata$Y,1)))
   wt = matrix(1, nrow=nrow(simdata), ncol=length(uy))
   wt[simdata$X1 < xc, uy <= Y10 ] = weight[1]
   wt[simdata$X1 < xc, uy> Y10 & uy<=Y90 ] = weight[2]
   wt[simdata$X1 < xc, uy > Y90  ] = weight[3]
   wt[simdata$X1 >= xc, uy <= Y10  ] = weight[4]
   wt[simdata$X1 >= xc, uy> Y10 & uy<=Y90  ] = weight[5]
   wt[simdata$X1 >= xc, uy > Y90 ] = weight[6]
   ## with known weights
   temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2), weight=wt)
   odsyx.sorfi[i,1:2]= as.matrix(temp$lambdagamma)[1,(length(temp$lambdagamma)-1):length(temp$lambdagamma)]
   odsyx.sorfi[i,3:4]= c(temp$gammase)
     d.odsyxsorfi[i]= sqrt(det(temp$gammavar))

  ## obtain bootstrap for standard error
   Xboot=rbind(cbind(simdata$X1, simdata$X2)[sample((1:nrow(simdata))[simdata$X1 <=xc], sum(popdata$X1<=xc), replace=T),],
               cbind(simdata$X1, simdata$X2)[sample((1:nrow(simdata))[simdata$X1 > xc], sum(popdata$X1>xc), replace=T),]) 
   nlg=length(temp$lambdagamma)
   Pboot=pred.sor(uy, Xboot,lambda=as.matrix(temp$lambdagamma)[1,-((nlg-1):nlg)], gamma=as.matrix(temp$lambdagamma)[1,((nlg-1):nlg)] ) 
   Yboot=apply(Pboot, 1,  function(i) uy[as.vector(rmultinom(1,size=1, i))==1])
   popdataboot=data.frame(X=Xboot, Y=Yboot); names(popdataboot)=c("X1","X2","Y")
   res=fun.bootgamma(popdataboot, xc, weight)
   odsyx.sorfi[i,3:4]= sqrt(diag(var(res$gammaboot)))
   d.odsyxsorfi[i]= sqrt(det(var(res$gammaboot)))
   odsyx.coss[i,]=c(fitcoss$par[2:3], sqrt(diag(var(res$cossboot))))
   d.odsyxcoss[i]= sqrt(det(var(res$cossboot)))

    ## ODS-YX with POR
      x2c=0
   simdata = rbind(popdata[sample((1:N)[Y<=Y10 & X2<x2c], size=n*0.2),], popdata[sample((1:N)[Y> Y10 & Y<= Y90 & X2<x2c], size=n*0.1),], 
                popdata[sample((1:N)[ Y>Y90 & X2<x2c ], size=n*0.2),],popdata[sample((1:N)[Y<=Y10 & X2>=x2c], size=n*0.2),], 
           popdata[sample((1:N)[Y> Y10 & Y<= Y90 & X2>=x2c], size=n*0.1),], popdata[sample((1:N)[ Y>Y90 & X2>=x2c], size=n*0.2),])  
   strata <- ifelse(simdata$X2<x2c,1,2) 
   fit.por <- optim(c(0,0,0.1, 0.1, 0,0), llk.por, y=simdata$Y, x=cbind(simdata$X1, simdata$X2), strata=strata, method="L-BFGS-B", hessian=T)
   odsyx.por[i,1:2] <- fit.por$par[3:4]
   odsyx.por[i,3:4] <- sqrt(diag(solve(fit.por$hessian)))[3:4]
   d.odsyxpor[i]= sqrt(det(solve(fit.por$hessian)[3:4,3:4])) 

   ## ODS-YX with SOR without weights
   temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2),strata=strata)
   odsyx.sor[i,1:2]= as.matrix(temp$lambdagamma)[1,(length(temp$lambdagamma)-1):length(temp$lambdagamma)]
   odsyx.sor[i,3:4]= c(temp$gammase) ## c(temp$lambdagamma[1,length(temp$lambdagamma)], temp$gammase)
   d.odsyxsor[i]= sqrt(det(temp$gammavar))

}



apply(srs.lm,2,mean)
apply(srs.hw,2,mean)
apply(srs.por,2,mean)
apply(srs.sor,2,mean)
apply(ods.lm,2,mean)
apply(ods.hw,2,mean)
apply(ods.hwfi,2,mean)
apply(ods.por,2,mean)
apply(ods.sor,2,mean)
apply(ods.sorfi,2,mean)
apply(odsyx.lm,2,mean)
apply(odsyx.coss,2,mean)
apply(odsyx.por,2,mean)
apply(odsyx.sor,2,mean)
apply(odsyx.sorfi,2,mean)

apply(srs.lm,2,sd)
apply(srs.por,2,sd)
apply(srs.hw,2,sd)
apply(srs.sor,2,sd)
apply(ods.lm,2,sd)
apply(ods.hw,2,sd)
apply(ods.hwfi,2,sd)
apply(ods.por,2,sd)
apply(ods.sor,2,sd)
apply(ods.sorfi,2,sd)
apply(odsyx.lm,2,sd)
apply(odsyx.por,2,sd)
apply(odsyx.coss,2,sd)
apply(odsyx.sor,2,sd)
apply(odsyx.sorfi,2,sd)


sum(srs.lm[,1]-1.96*srs.lm[,3] <= 0.1 &  srs.lm[,1]+ 1.96*srs.lm[,3] >= 0.1 )/nsim
sum(srs.hw[,1]-1.96*srs.hw[,3] <= 0.1 &  srs.hw[,1]+ 1.96*srs.hw[,3] >= 0.1 )/nsim
sum(srs.por[,1]-1.96*srs.por[,3] <= 0.1 &  srs.por[,1]+ 1.96*srs.por[,3] >= 0.1 )/nsim
sum(srs.sor[,1]-1.96*srs.sor[,3] <= 0.1 &  srs.sor[,1]+ 1.96*srs.sor[,3] >= 0.1 )/nsim
sum(ods.lm[,1]-1.96*ods.lm[,3] <= 0.1 &  ods.lm[,1]+ 1.96*ods.lm[,3] >= 0.1 )/nsim
sum(ods.hw[,1]-1.96*ods.hw[,3] <= 0.1 &  ods.hw[,1]+ 1.96*ods.hw[,3] >= 0.1 )/nsim
sum(ods.hwfi[,1]-1.96*ods.hwfi[,3] <= 0.1 &  ods.hwfi[,1]+ 1.96*ods.hwfi[,3] >= 0.1 )/nsim
sum(ods.por[,1]-1.96*ods.por[,3] <= 0.1 &  ods.por[,1]+ 1.96*ods.por[,3] >= 0.1 )/nsim
sum(ods.sor[,1]-1.96*ods.sor[,3] <= 0.1 &  ods.sor[,1]+ 1.96*ods.sor[,3] >= 0.1 )/nsim
sum(ods.sorfi[,1]-1.96*ods.sorfi[,3] <= 0.1 &  ods.sorfi[,1]+ 1.96*ods.sorfi[,3] >= 0.1 )/nsim
sum(odsyx.lm[,1]-1.96*odsyx.lm[,3] <= 0.1 &  odsyx.lm[,1]+ 1.96*odsyx.lm[,3] >= 0.1 )/nsim
sum(odsyx.coss[,1]-1.96*odsyx.coss[,3] <= 0.1 &  odsyx.coss[,1]+ 1.96*odsyx.coss[,3] >= 0.1 )/nsim
sum(odsyx.por[,1]-1.96*odsyx.por[,3] <= 0.1 &  odsyx.por[,1]+ 1.96*odsyx.por[,3] >= 0.1 )/nsim
sum(odsyx.sor[,1]-1.96*odsyx.sor[,3] <= 0.1 &  odsyx.sor[,1]+ 1.96*odsyx.sor[,3] >= 0.1 )/nsim
sum(odsyx.sorfi[,1]-1.96*odsyx.sorfi[,3] <= 0.1 &  odsyx.sorfi[,1]+ 1.96*odsyx.sorfi[,3] >= 0.1 )/nsim


sum(srs.lm[,2]-1.96*srs.lm[,4] <= 0.1 &  srs.lm[,2]+ 1.96*srs.lm[,4] >= 0.1 )/nsim
sum(srs.hw[,2]-1.96*srs.hw[,4] <= 0.1 &  srs.hw[,2]+ 1.96*srs.hw[,4] >= 0.1 )/nsim
sum(srs.por[,2]-1.96*srs.por[,4] <= 0.1 &  srs.por[,2]+ 1.96*srs.por[,4] >= 0.1 )/nsim
sum(srs.sor[,2]-1.96*srs.sor[,4] <= 0.1 &  srs.sor[,2]+ 1.96*srs.sor[,4] >= 0.1 )/nsim
sum(ods.lm[,2]-1.96*ods.lm[,4] <= 0.1 &  ods.lm[,2]+ 1.96*ods.lm[,4] >= 0.1 )/nsim
sum(ods.hw[,2]-1.96*ods.hw[,4] <= 0.1 &  ods.hw[,2]+ 1.96*ods.hw[,4] >= 0.1 )/nsim
sum(ods.hwfi[,2]-1.96*ods.hwfi[,4] <= 0.1 &  ods.hwfi[,2]+ 1.96*ods.hwfi[,4] >= 0.1 )/nsim
sum(ods.por[,2]-1.96*ods.por[,4] <= 0.1 &  ods.por[,2]+ 1.96*ods.por[,4] >= 0.1 )/nsim
sum(ods.sor[,2]-1.96*ods.sor[,4] <= 0.1 &  ods.sor[,2]+ 1.96*ods.sor[,4] >= 0.1 )/nsim
sum(ods.sorfi[,2]-1.96*ods.sorfi[,4] <= 0.1 &  ods.sorfi[,2]+ 1.96*ods.sorfi[,4] >= 0.1 )/nsim
sum(odsyx.lm[,2]-1.96*odsyx.lm[,4] <= 0.1 &  odsyx.lm[,2]+ 1.96*odsyx.lm[,4] >= 0.1 )/nsim
sum(odsyx.coss[,2]-1.96*odsyx.coss[,4] <= 0.1 &  odsyx.coss[,2]+ 1.96*odsyx.coss[,4] >= 0.1 )/nsim
sum(odsyx.por[,2]-1.96*odsyx.por[,4] <= 0.1 &  odsyx.por[,2]+ 1.96*odsyx.por[,4] >= 0.1 )/nsim
sum(odsyx.sor[,2]-1.96*odsyx.sor[,4] <= 0.1 &  odsyx.sor[,2]+ 1.96*odsyx.sor[,4] >= 0.1 )/nsim
sum(odsyx.sorfi[,2]-1.96*odsyx.sorfi[,4] <= 0.1 &  odsyx.sorfi[,2]+ 1.96*odsyx.sorfi[,4] >= 0.1 )/nsim

1-sum(srs.lm[,1]-1.96*srs.lm[,3] <= 0 &  srs.lm[,1]+ 1.96*srs.lm[,3] >= 0 )/nsim
1-sum(srs.hw[,1]-1.96*srs.hw[,3] <= 0 &  srs.hw[,1]+ 1.96*srs.hw[,3] >= 0 )/nsim
1-sum(srs.por[,1]-1.96*srs.por[,3] <= 0 &  srs.por[,1]+ 1.96*srs.por[,3] >= 0 )/nsim
1-sum(srs.sor[,1]-1.96*srs.sor[,3] <= 0 &  srs.sor[,1]+ 1.96*srs.sor[,3] >= 0 )/nsim
1-sum(ods.lm[,1]-1.96*ods.lm[,3] <= 0 &  ods.lm[,1]+ 1.96*ods.lm[,3] >= 0 )/nsim
1-sum(ods.hw[,1]-1.96*ods.hw[,3] <= 0 &  ods.hw[,1]+ 1.96*ods.hw[,3] >= 0 )/nsim
1-sum(ods.hwfi[,1]-1.96*ods.hwfi[,3] <= 0 &  ods.hwfi[,1]+ 1.96*ods.hwfi[,3] >= 0 )/nsim
1-sum(ods.por[,1]-1.96*ods.por[,3] <= 0 &  ods.por[,1]+ 1.96*ods.por[,3] >= 0 )/nsim
1-sum(ods.sor[,1]-1.96*ods.sor[,3] <= 0 &  ods.sor[,1]+ 1.96*ods.sor[,3] >= 0 )/nsim
1-sum(ods.sorfi[,1]-1.96*ods.sorfi[,3] <= 0 &  ods.sorfi[,1]+ 1.96*ods.sorfi[,3] >= 0 )/nsim
1-sum(odsyx.lm[,1]-1.96*odsyx.lm[,3] <= 0 &  odsyx.lm[,1]+ 1.96*odsyx.lm[,3] >= 0 )/nsim
1-sum(odsyx.coss[,1]-1.96*odsyx.coss[,3] <= 0 &  odsyx.coss[,1]+ 1.96*odsyx.coss[,3] >= 0 )/nsim
1-sum(odsyx.por[,1]-1.96*odsyx.por[,3] <= 0 &  odsyx.por[,1]+ 1.96*odsyx.por[,3] >= 0 )/nsim
1-sum(odsyx.sor[,1]-1.96*odsyx.sor[,3] <= 0 &  odsyx.sor[,1]+ 1.96*odsyx.sor[,3] >= 0 )/nsim
1-sum(odsyx.sorfi[,1]-1.96*odsyx.sorfi[,3] <= 0 &  odsyx.sorfi[,1]+ 1.96*odsyx.sorfi[,3] >= 0 )/nsim

1-sum(srs.lm[,2]-1.96*srs.lm[,4] <= 0 &  srs.lm[,2]+ 1.96*srs.lm[,4] >= 0 )/nsim
1-sum(srs.hw[,2]-1.96*srs.hw[,4] <= 0 &  srs.hw[,2]+ 1.96*srs.hw[,4] >= 0 )/nsim
1-sum(srs.por[,2]-1.96*srs.por[,4] <= 0 &  srs.por[,2]+ 1.96*srs.por[,4] >= 0 )/nsim
1-sum(srs.sor[,2]-1.96*srs.sor[,4] <= 0 &  srs.sor[,2]+ 1.96*srs.sor[,4] >= 0 )/nsim
1-sum(ods.lm[,2]-1.96*ods.lm[,4] <= 0 &  ods.lm[,2]+ 1.96*ods.lm[,4] >= 0 )/nsim
1-sum(ods.hw[,2]-1.96*ods.hw[,4] <= 0 &  ods.hw[,2]+ 1.96*ods.hw[,4] >= 0 )/nsim
1-sum(ods.hwfi[,2]-1.96*ods.hwfi[,4] <= 0 &  ods.hwfi[,2]+ 1.96*ods.hwfi[,4] >= 0 )/nsim
1-sum(ods.por[,2]-1.96*ods.por[,4] <= 0 &  ods.por[,2]+ 1.96*ods.por[,4] >= 0 )/nsim
1-sum(ods.sor[,2]-1.96*ods.sor[,4] <= 0 &  ods.sor[,2]+ 1.96*ods.sor[,4] >= 0 )/nsim
1-sum(ods.sorfi[,2]-1.96*ods.sorfi[,4] <= 0 &  ods.sorfi[,2]+ 1.96*ods.sorfi[,4] >= 0 )/nsim
1-sum(odsyx.lm[,2]-1.96*odsyx.lm[,4] <= 0 &  odsyx.lm[,2]+ 1.96*odsyx.lm[,4] >= 0 )/nsim
1-sum(odsyx.coss[,2]-1.96*odsyx.coss[,4] <= 0 &  odsyx.coss[,2]+ 1.96*odsyx.coss[,4] >= 0 )/nsim
1-sum(odsyx.por[,2]-1.96*odsyx.por[,4] <= 0 &  odsyx.por[,2]+ 1.96*odsyx.por[,4] >= 0 )/nsim
1-sum(odsyx.sor[,2]-1.96*odsyx.sor[,4] <= 0 &  odsyx.sor[,2]+ 1.96*odsyx.sor[,4] >= 0 )/nsim
1-sum(odsyx.sorfi[,2]-1.96*odsyx.sorfi[,4] <= 0 &  odsyx.sorfi[,2]+ 1.96*odsyx.sorfi[,4] >= 0 )/nsim


1-sum((srs.lm[,1]-1.96*srs.lm[,3] <= 0 &  srs.lm[,1]+ 1.96*srs.lm[,3] >= 0) |( srs.lm[,2]-1.96*srs.lm[,4] <= 0 &  srs.lm[,2]+ 1.96*srs.lm[,4] >= 0) )/nsim
1-sum((srs.hw[,1]-1.96*srs.hw[,3] <= 0 &  srs.hw[,1]+ 1.96*srs.hw[,3] >= 0) | (srs.hw[,2]-1.96*srs.hw[,4] <= 0 &  srs.hw[,2]+ 1.96*srs.hw[,4] >= 0)  )/nsim
1-sum((srs.por[,1]-1.96*srs.por[,3] <= 0 &  srs.por[,1]+ 1.96*srs.por[,3] >= 0) | (srs.por[,2]-1.96*srs.por[,4] <= 0 &  srs.por[,2]+ 1.96*srs.por[,4] >= 0) )/nsim
1-sum((ods.hw[,1]-1.96*ods.hw[,3] <= 0 &  ods.hw[,1]+ 1.96*ods.hw[,3] >= 0) |(ods.hw[,2]-1.96*ods.hw[,4] <= 0 &  ods.hw[,2]+ 1.96*ods.hw[,4] >= 0) )/nsim
1-sum((ods.hwfi[,1]-1.96*ods.hwfi[,3] <= 0 &  ods.hwfi[,1]+ 1.96*ods.hwfi[,3] >= 0) |(ods.hwfi[,2]-1.96*ods.hwfi[,4] <= 0 &  ods.hwfi[,2]+ 1.96*ods.hwfi[,4] >= 0) )/nsim
1-sum((ods.sor[,1]-1.96*ods.sor[,3] <= 0 &  ods.sor[,1]+ 1.96*ods.sor[,3] >= 0) |(ods.sor[,2]-1.96*ods.sor[,4] <= 0 &  ods.sor[,2]+ 1.96*ods.sor[,4] >= 0) )/nsim
1-sum((ods.sorfi[,1]-1.96*ods.sorfi[,3] <= 0 &  ods.sorfi[,1]+ 1.96*ods.sorfi[,3] >= 0) |(ods.sorfi[,2]-1.96*ods.sorfi[,4] <= 0 &  ods.sorfi[,2]+ 1.96*ods.sorfi[,4] >= 0) )/nsim
1-sum((odsyx.coss[,1]-1.96*odsyx.coss[,3] <= 0 &  odsyx.coss[,1]+ 1.96*odsyx.coss[,3] >= 0) |( odsyx.coss[,2]-1.96*odsyx.coss[,4] <= 0 &  odsyx.coss[,2]+ 1.96*odsyx.coss[,4] >= 0 ))/nsim
1-sum((odsyx.sor[,1]-1.96*odsyx.sor[,3] <= 0 &  odsyx.sor[,1]+ 1.96*odsyx.sor[,3] >= 0) |( odsyx.sor[,2]-1.96*odsyx.sor[,4] <= 0 &  odsyx.sor[,2]+ 1.96*odsyx.sor[,4] >= 0 ))/nsim
1-sum((odsyx.sorfi[,1]-1.96*odsyx.sorfi[,3] <= 0 &  odsyx.sorfi[,1]+ 1.96*odsyx.sorfi[,3] >= 0) |( odsyx.sorfi[,2]-1.96*odsyx.sorfi[,4] <= 0 &  odsyx.sorfi[,2]+ 1.96*odsyx.sorfi[,4] >= 0))/nsim

mean(d.srslm)
mean(d.srshw)
mean(d.srspor)
mean(d.srssor)
mean(d.odslm)
mean(d.odshw)
mean(d.odshwfi)
mean(d.odspor)
mean(d.odssor)
mean(d.odssorfi)
mean(d.odsyxpor)
mean(d.odsyxsor)
mean(d.odsyxsorfi)

#####################################
## outcome is truncated normal 
#####################################
N=1000000
X2 <- rnorm(N) 
X1 <- rbinom(N, size=1, p=0.15) 
Y<-rtruncnorm(N, a=0, mean=0.1*X1+0.1*X2) ## use Truncated normal distribution
popdata <-data.frame(X1=X1,X2=X2, Y=Y)
N <- nrow(popdata)
truncreg(Y~X1+X2, data=popdata, point=0, direction="left")

## Endogenous outcome-dependent selective sampling
n=500
nsim=500
srs.lm = srs.por = srs.hw=srs.hwfi=matrix(0, nsim,4)
srs.sor = matrix(0, nsim,4)
ods.lm= ods.por= ods.hw=ods.hwfi=matrix(0, nsim,4)
ods.sor=ods.sorfi= matrix(0, nsim, 4)
odsyx.lm= odsyx.por=odsyx.coss= matrix(0, nsim,4)
odsyx.sor= odsyx.sorfi=matrix(0, nsim, 4)
d.srslm=d.srspor=d.srshw=d.srshwfi=d.srssor=d.odslm=d.odspor=d.odshw=d.odshwfi=d.odssor=d.odssorfi=numeric(nsim)
d.odsyxlm=d.odsyxpor=d.odsyxsor=d.odsyxsorfi=d.odsyxcoss=numeric(nsim)
Y10= quantile(Y, 0.1); Y90=quantile(Y, 0.9)

## Simple random sampling
for (i in 1:nsim) {
  print(i)
   ## simple random sample with linear regression 
   simdata.srs= popdata[sample(1:N, size=n),]
   sim.reg<-truncreg(Y~X1+X2, data=simdata.srs, point=0, direction="left")
   srs.lm[i,1:2] <- sim.reg$coef[2:3] ##/sim.reg$coef["sigma"]
   srs.lm[i,3:4] <- summary(sim.reg)$coef[2:3,2]
   d.srslm[i]= sqrt(det(vcov(sim.reg)[2:3,2:3])) 



   ## simple random ssample with POR   
   fit.por <- optim(c(0,1,1, 0), llk.por, y=simdata.srs$Y, x=cbind(1, simdata.srs$X1, simdata.srs$X2), method="L-BFGS-B", hessian=T)
   srs.por[i,1:2] <- fit.por$par[2:3]
   srs.por[i,3:4] <- sqrt(diag(solve(fit.por$hessian)))[2:3]
   d.srspor[i]= sqrt(det(solve(fit.por$hessian)[2:3,2:3]))


   ## Simple random sampling with HW
   fithw=fun.hw(y=simdata.srs$Y, x=cbind(1,simdata.srs$X1, simdata.srs$X2), stratapoints=c(Y10, Y90), weights=c(1,1,1), dist="gaussian", par0=c(0,0.1,0.1,0))
   srs.hwfi[i,]=c(fithw$par[2:3], sqrt(diag(solve(fithw$hessian)))[2:3])
   d.srshwfi[i]= sqrt(det(solve(fithw$hessian)[2:3,2:3]))

   ## SRS with HW without weights
   fithw=fun.hw(y=simdata.srs$Y, x=cbind(1,simdata.srs$X1, simdata.srs$X2), stratapoints=c(Y10, Y90),  dist="gaussian", par0=c(0,0.1,0.1,0, 0,0))
   srs.hw[i,]= c(fithw$par[2:3], sqrt(diag(solve(fithw$hessian)))[2:3])
   d.srshw[i]= sqrt(det(solve(fithw$hessian)[2:3,2:3]))

   ## simple random sample with sor 
   temp.srs=fitsor(y=round(simdata.srs$Y,1), x=cbind(simdata.srs$X1, simdata.srs$X2))
   srs.sor[i,1:2]= as.matrix(temp.srs$lambdagamma)[1,(length(temp.srs$lambdagamma)-1):length(temp.srs$lambdagamma)]
   srs.sor[i,3:4]= c(temp.srs$gammase) 
   d.srssor[i]= sqrt(det(temp.srs$gammavar))
}


   #####################  ODS on y only ########################
for (i in 1:nsim) {
  print(i)
   ## ODS with OLS
   simdata = rbind(popdata[sample((1:N)[Y<=Y10], size=n*0.4),], popdata[sample((1:N)[Y> Y10 & Y<= Y90 ], size=n*0.2),], popdata[sample((1:N)[ Y>Y90 ], size=n*0.4),])
   sim.reg<-truncreg(Y~X1+X2, data=simdata, point=0, direction="left")
   ods.lm[i,1:2] <- sim.reg$coef[2:3] ##/sim.reg$coef["sigma"] 
   ods.lm[i,3:4] <- summary(sim.reg)$coef[2:3,2]
   d.odslm[i]= sqrt(det(vcov(sim.reg)[2:3,2:3])) 
  

   ## ods with POR   
   fit.por <- optim(c(0,1,1, 0), llk.por, y=simdata$Y, x=cbind(1, simdata$X1, simdata$X2), method="L-BFGS-B", hessian=T)
   ods.por[i,1:2] <- fit.por$par[2:3]
   ods.por[i,3:4] <- sqrt(diag(solve(fit.por$hessian)))[2:3]
    d.odspor[i]= sqrt(det(solve(fit.por$hessian)[2:3,2:3]))

      ## ODS  with HW+ weights
    fithw=fun.hw(y=simdata$Y, x=cbind(1,simdata$X1, simdata$X2), stratapoints=c(Y10, Y90), weights=c(4,0.25,4), dist="gaussian", par0=c(0,0.1,0.1,0))
    ods.hwfi[i,]=c(fithw$par[2:3], sqrt(diag(solve(fithw$hessian)))[2:3])
    d.odshwfi[i]= sqrt(det(solve(fithw$hessian)[2:3,2:3]))

   ## ODS with HW without weights
   fithw=fun.hw(y=simdata$Y, x=cbind(1,simdata$X1, simdata$X2), stratapoints=c(Y10, Y90),  dist="gaussian", par0=c(0,0.1,0.1,0, 0,0))
   ods.hw[i,]= c(fithw$par[2:3], sqrt(diag(solve(fithw$hessian)))[2:3])
   d.odshw[i]= sqrt(det(solve(fithw$hessian)[2:3,2:3]))

    ## ODS with sor 
   temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2))
   ods.sor[i,1:2]= as.matrix(temp$lambdagamma)[1,(length(temp$lambdagamma)-1):length(temp$lambdagamma)]
   ods.sor[i,3:4]= c(temp$gammase) ## c(temp$lambdagamma[1,length(temp$lambdagamma)], temp$gammase)
   d.odssor[i]= sqrt(det(temp$gammavar))

   ## ODS with SOR and weights
    uy= sort(unique(round(simdata$Y,1)))
   wt <- matrix(1, nrow(simdata), length(uy))
   wt[, uy <= Y10] = 0.4/sum(Y<=Y10)*N
   wt[, uy> Y10 & uy<=Y90] = 0.2/sum(Y>Y10 & Y<=Y90)*N
   wt[, uy > Y90] = 0.4/sum(Y > Y90)*N
   temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2), weight=wt)
   ods.sorfi[i,1:2]= as.matrix(temp$lambdagamma)[1,(length(temp$lambdagamma)-1):length(temp$lambdagamma)]
   ods.sorfi[i,3:4]= c(temp$gammase)
   d.odssorfi[i]= sqrt(det(temp$gammavar))
}


   #####################  ODS on both y and x ########################
for (i in 1:nsim) {
   print(i)
     xc=0.5 
   if (i==1) {
   weight = c(0.2/sum(Y<=Y10 & X1<xc), 0.1/sum(Y> Y10 & Y<= Y90 & X1<xc), 0.2/sum( Y>Y90 & X1<xc ),
          0.2/sum(Y<=Y10 & X1>=xc), 0.1/sum(Y> Y10 & Y<= Y90 & X1>=xc), 0.2/sum(Y>Y90 & X1>=xc))*N
   weight= weight/sum(weight)
             }
   simdata = rbind(popdata[sample((1:N)[Y<=Y10 & X1<xc], size=n*0.2),], popdata[sample((1:N)[Y> Y10 & Y<= Y90 & X1<xc], size=n*0.1),], 
                popdata[sample((1:N)[ Y>Y90 & X1<xc ], size=n*0.2),],popdata[sample((1:N)[Y<=Y10 & X1>=xc], size=n*0.2),], 
           popdata[sample((1:N)[Y> Y10 & Y<= Y90 & X1>=xc], size=n*0.1),], popdata[sample((1:N)[ Y>Y90 & X1>=xc], size=n*0.2),])
   
   ## ODS with regular analysis
   sim.reg<-truncreg(Y~X1+X2, data=simdata, point=0, direction="left")
   odsyx.lm[i,1:2] <- sim.reg$coef[2:3] 
   odsyx.lm[i,3:4] <- summary(sim.reg)$coef[2:3,2]
   d.odsyxlm[i]= sqrt(det(vcov(sim.reg)[2:3,2:3])) 

   ## ODS  with COSS + weights   
  cosswt=matrix(1, nrow=nrow(simdata), ncol=3)
  for (j in 1:nrow(simdata)) {
   if (simdata$X1[j] <xc) cosswt[j, ]=weight[1:3] else 
    cosswt[j,]=weight[4:6]
          }
    fitcoss=fun.coss(y=simdata$Y, x=cbind(1,simdata$X1, simdata$X2), stratapoints=c(Y10, Y90), weights=cosswt, dist="gaussian", par0=c(0,0.1,0.1,0))
   odsyx.coss[i,]=c(fitcoss$par[2:3], sqrt(diag(solve(fitcoss$hessian)))[2:3])
    d.odsyxcoss[i]= sqrt(det(solve(fitcoss$hessian)[2:3,2:3]))


  ## ods with SOR +weights
   uy= sort(unique(round(simdata$Y,1)))
   wt = matrix(1, nrow=nrow(simdata), ncol=length(uy))
   wt[simdata$X1 < xc, uy <= Y10 ] = weight[1]
   wt[simdata$X1 < xc, uy> Y10 & uy<=Y90 ] = weight[2]
   wt[simdata$X1 < xc, uy > Y90  ] = weight[3]
   wt[simdata$X1 >= xc, uy <= Y10  ] = weight[4]
   wt[simdata$X1 >= xc, uy> Y10 & uy<=Y90  ] = weight[5]
   wt[simdata$X1 >= xc, uy > Y90 ] = weight[6]
   ## with known weights
   temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2), weight=wt)
   odsyx.sorfi[i,1:2]= as.matrix(temp$lambdagamma)[1,(length(temp$lambdagamma)-1):length(temp$lambdagamma)]
   odsyx.sorfi[i,3:4]= c(temp$gammase)
     d.odsyxsorfi[i]= sqrt(det(temp$gammavar))

  ## obtain bootstrap for standard error
   Xboot=rbind(cbind(simdata$X1, simdata$X2)[sample((1:nrow(simdata))[simdata$X1 <=xc], sum(popdata$X1<=xc), replace=T),],
               cbind(simdata$X1, simdata$X2)[sample((1:nrow(simdata))[simdata$X1 > xc], sum(popdata$X1>xc), replace=T),]) 
   nlg=length(temp$lambdagamma)
   Pboot=pred.sor(uy, Xboot,lambda=as.matrix(temp$lambdagamma)[1,-((nlg-1):nlg)], gamma=as.matrix(temp$lambdagamma)[1,((nlg-1):nlg)] ) 
   Yboot=apply(Pboot, 1,  function(i) uy[as.vector(rmultinom(1,size=1, i))==1])
   popdataboot=data.frame(X=Xboot, Y=Yboot); names(popdataboot)=c("X1","X2","Y")
   res=fun.bootgamma(popdataboot, xc, weight)
   odsyx.sorfi[i,3:4]= sqrt(diag(var(res$gammaboot)))
   d.odsyxsorfi[i]= sqrt(det(var(res$gammaboot)))
   odsyx.coss[i,]=c(fitcoss$par[2:3], sqrt(diag(var(res$cossboot))))
   d.odsyxcoss[i]= sqrt(det(var(res$cossboot)))

    ## ODS-YX with POR
      x2c=0
   simdata = rbind(popdata[sample((1:N)[Y<=Y10 & X2<x2c], size=n*0.2),], popdata[sample((1:N)[Y> Y10 & Y<= Y90 & X2<x2c], size=n*0.1),], 
                popdata[sample((1:N)[ Y>Y90 & X2<x2c ], size=n*0.2),],popdata[sample((1:N)[Y<=Y10 & X2>=x2c], size=n*0.2),], 
           popdata[sample((1:N)[Y> Y10 & Y<= Y90 & X2>=x2c], size=n*0.1),], popdata[sample((1:N)[ Y>Y90 & X2>=x2c], size=n*0.2),])  
   strata <- ifelse(simdata$X2<x2c,1,2) 
   fit.por <- optim(c(0,0,0.1, 0.1, 0,0), llk.por, y=simdata$Y, x=cbind(simdata$X1, simdata$X2), strata=strata, method="L-BFGS-B", hessian=T)
   odsyx.por[i,1:2] <- fit.por$par[3:4]
   odsyx.por[i,3:4] <- sqrt(diag(solve(fit.por$hessian)))[3:4]
   d.odsyxpor[i]= sqrt(det(solve(fit.por$hessian)[3:4,3:4])) 

   ## ODS-YX with SOR without weights
   temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2),strata=strata)
   odsyx.sor[i,1:2]= as.matrix(temp$lambdagamma)[1,(length(temp$lambdagamma)-1):length(temp$lambdagamma)]
   odsyx.sor[i,3:4]= c(temp$gammase) 
   d.odsyxsor[i]= sqrt(det(temp$gammavar))

}



apply(srs.lm,2,mean)
apply(srs.hw,2,mean)
apply(srs.por,2,mean)
apply(srs.sor,2,mean)
apply(ods.lm,2,mean)
apply(ods.hw,2,mean)
apply(ods.hwfi,2,mean)
apply(ods.por,2,mean)
apply(ods.sor,2,mean)
apply(ods.sorfi,2,mean)
apply(odsyx.lm,2,mean)
apply(odsyx.coss,2,mean)
apply(odsyx.por,2,mean)
apply(odsyx.sor,2,mean)
apply(odsyx.sorfi,2,mean)

apply(srs.lm,2,sd)
apply(srs.por,2,sd)
apply(srs.hw,2,sd)
apply(srs.sor,2,sd)
apply(ods.lm,2,sd)
apply(ods.hw,2,sd)
apply(ods.hwfi,2,sd)
apply(ods.por,2,sd)
apply(ods.sor,2,sd)
apply(ods.sorfi,2,sd)
apply(odsyx.lm,2,sd)
apply(odsyx.por,2,sd)
apply(odsyx.coss,2,sd)
apply(odsyx.sor,2,sd)
apply(odsyx.sorfi,2,sd)


sum(srs.lm[,1]-1.96*srs.lm[,3] <= 0.1 &  srs.lm[,1]+ 1.96*srs.lm[,3] >= 0.1 )/nsim
sum(srs.hw[,1]-1.96*srs.hw[,3] <= 0.1 &  srs.hw[,1]+ 1.96*srs.hw[,3] >= 0.1 )/nsim
sum(srs.por[,1]-1.96*srs.por[,3] <= 0.1 &  srs.por[,1]+ 1.96*srs.por[,3] >= 0.1 )/nsim
sum(srs.sor[,1]-1.96*srs.sor[,3] <= 0.1 &  srs.sor[,1]+ 1.96*srs.sor[,3] >= 0.1 )/nsim
sum(ods.lm[,1]-1.96*ods.lm[,3] <= 0.1 &  ods.lm[,1]+ 1.96*ods.lm[,3] >= 0.1 )/nsim
sum(ods.hw[,1]-1.96*ods.hw[,3] <= 0.1 &  ods.hw[,1]+ 1.96*ods.hw[,3] >= 0.1 )/nsim
sum(ods.hwfi[,1]-1.96*ods.hwfi[,3] <= 0.1 &  ods.hwfi[,1]+ 1.96*ods.hwfi[,3] >= 0.1 )/nsim
sum(ods.por[,1]-1.96*ods.por[,3] <= 0.1 &  ods.por[,1]+ 1.96*ods.por[,3] >= 0.1 )/nsim
sum(ods.sor[,1]-1.96*ods.sor[,3] <= 0.1 &  ods.sor[,1]+ 1.96*ods.sor[,3] >= 0.1 )/nsim
sum(ods.sorfi[,1]-1.96*ods.sorfi[,3] <= 0.1 &  ods.sorfi[,1]+ 1.96*ods.sorfi[,3] >= 0.1 )/nsim
sum(odsyx.lm[,1]-1.96*odsyx.lm[,3] <= 0.1 &  odsyx.lm[,1]+ 1.96*odsyx.lm[,3] >= 0.1 )/nsim
sum(odsyx.coss[,1]-1.96*odsyx.coss[,3] <= 0.1 &  odsyx.coss[,1]+ 1.96*odsyx.coss[,3] >= 0.1 )/nsim
sum(odsyx.por[,1]-1.96*odsyx.por[,3] <= 0.1 &  odsyx.por[,1]+ 1.96*odsyx.por[,3] >= 0.1 )/nsim
sum(odsyx.sor[,1]-1.96*odsyx.sor[,3] <= 0.1 &  odsyx.sor[,1]+ 1.96*odsyx.sor[,3] >= 0.1 )/nsim
sum(odsyx.sorfi[,1]-1.96*odsyx.sorfi[,3] <= 0.1 &  odsyx.sorfi[,1]+ 1.96*odsyx.sorfi[,3] >= 0.1 )/nsim


sum(srs.lm[,2]-1.96*srs.lm[,4] <= 0.1 &  srs.lm[,2]+ 1.96*srs.lm[,4] >= 0.1 )/nsim
sum(srs.hw[,2]-1.96*srs.hw[,4] <= 0.1 &  srs.hw[,2]+ 1.96*srs.hw[,4] >= 0.1 )/nsim
sum(srs.por[,2]-1.96*srs.por[,4] <= 0.1 &  srs.por[,2]+ 1.96*srs.por[,4] >= 0.1 )/nsim
sum(srs.sor[,2]-1.96*srs.sor[,4] <= 0.1 &  srs.sor[,2]+ 1.96*srs.sor[,4] >= 0.1 )/nsim
sum(ods.lm[,2]-1.96*ods.lm[,4] <= 0.1 &  ods.lm[,2]+ 1.96*ods.lm[,4] >= 0.1 )/nsim
sum(ods.hw[,2]-1.96*ods.hw[,4] <= 0.1 &  ods.hw[,2]+ 1.96*ods.hw[,4] >= 0.1 )/nsim
sum(ods.hwfi[,2]-1.96*ods.hwfi[,4] <= 0.1 &  ods.hwfi[,2]+ 1.96*ods.hwfi[,4] >= 0.1 )/nsim
sum(ods.por[,2]-1.96*ods.por[,4] <= 0.1 &  ods.por[,2]+ 1.96*ods.por[,4] >= 0.1 )/nsim
sum(ods.sor[,2]-1.96*ods.sor[,4] <= 0.1 &  ods.sor[,2]+ 1.96*ods.sor[,4] >= 0.1 )/nsim
sum(ods.sorfi[,2]-1.96*ods.sorfi[,4] <= 0.1 &  ods.sorfi[,2]+ 1.96*ods.sorfi[,4] >= 0.1 )/nsim
sum(odsyx.lm[,2]-1.96*odsyx.lm[,4] <= 0.1 &  odsyx.lm[,2]+ 1.96*odsyx.lm[,4] >= 0.1 )/nsim
sum(odsyx.coss[,2]-1.96*odsyx.coss[,4] <= 0.1 &  odsyx.coss[,2]+ 1.96*odsyx.coss[,4] >= 0.1 )/nsim
sum(odsyx.por[,2]-1.96*odsyx.por[,4] <= 0.1 &  odsyx.por[,2]+ 1.96*odsyx.por[,4] >= 0.1 )/nsim
sum(odsyx.sor[,2]-1.96*odsyx.sor[,4] <= 0.1 &  odsyx.sor[,2]+ 1.96*odsyx.sor[,4] >= 0.1 )/nsim
sum(odsyx.sorfi[,2]-1.96*odsyx.sorfi[,4] <= 0.1 &  odsyx.sorfi[,2]+ 1.96*odsyx.sorfi[,4] >= 0.1 )/nsim

1-sum(srs.lm[,1]-1.96*srs.lm[,3] <= 0 &  srs.lm[,1]+ 1.96*srs.lm[,3] >= 0 )/nsim
1-sum(srs.hw[,1]-1.96*srs.hw[,3] <= 0 &  srs.hw[,1]+ 1.96*srs.hw[,3] >= 0 )/nsim
1-sum(srs.por[,1]-1.96*srs.por[,3] <= 0 &  srs.por[,1]+ 1.96*srs.por[,3] >= 0 )/nsim
1-sum(srs.sor[,1]-1.96*srs.sor[,3] <= 0 &  srs.sor[,1]+ 1.96*srs.sor[,3] >= 0 )/nsim
1-sum(ods.lm[,1]-1.96*ods.lm[,3] <= 0 &  ods.lm[,1]+ 1.96*ods.lm[,3] >= 0 )/nsim
1-sum(ods.hw[,1]-1.96*ods.hw[,3] <= 0 &  ods.hw[,1]+ 1.96*ods.hw[,3] >= 0 )/nsim
1-sum(ods.hwfi[,1]-1.96*ods.hwfi[,3] <= 0 &  ods.hwfi[,1]+ 1.96*ods.hwfi[,3] >= 0 )/nsim
1-sum(ods.por[,1]-1.96*ods.por[,3] <= 0 &  ods.por[,1]+ 1.96*ods.por[,3] >= 0 )/nsim
1-sum(ods.sor[,1]-1.96*ods.sor[,3] <= 0 &  ods.sor[,1]+ 1.96*ods.sor[,3] >= 0 )/nsim
1-sum(ods.sorfi[,1]-1.96*ods.sorfi[,3] <= 0 &  ods.sorfi[,1]+ 1.96*ods.sorfi[,3] >= 0 )/nsim
1-sum(odsyx.lm[,1]-1.96*odsyx.lm[,3] <= 0 &  odsyx.lm[,1]+ 1.96*odsyx.lm[,3] >= 0 )/nsim
1-sum(odsyx.coss[,1]-1.96*odsyx.coss[,3] <= 0 &  odsyx.coss[,1]+ 1.96*odsyx.coss[,3] >= 0 )/nsim
1-sum(odsyx.por[,1]-1.96*odsyx.por[,3] <= 0 &  odsyx.por[,1]+ 1.96*odsyx.por[,3] >= 0 )/nsim
1-sum(odsyx.sor[,1]-1.96*odsyx.sor[,3] <= 0 &  odsyx.sor[,1]+ 1.96*odsyx.sor[,3] >= 0 )/nsim
1-sum(odsyx.sorfi[,1]-1.96*odsyx.sorfi[,3] <= 0 &  odsyx.sorfi[,1]+ 1.96*odsyx.sorfi[,3] >= 0 )/nsim

1-sum(srs.lm[,2]-1.96*srs.lm[,4] <= 0 &  srs.lm[,2]+ 1.96*srs.lm[,4] >= 0 )/nsim
1-sum(srs.hw[,2]-1.96*srs.hw[,4] <= 0 &  srs.hw[,2]+ 1.96*srs.hw[,4] >= 0 )/nsim
1-sum(srs.por[,2]-1.96*srs.por[,4] <= 0 &  srs.por[,2]+ 1.96*srs.por[,4] >= 0 )/nsim
1-sum(srs.sor[,2]-1.96*srs.sor[,4] <= 0 &  srs.sor[,2]+ 1.96*srs.sor[,4] >= 0 )/nsim
1-sum(ods.lm[,2]-1.96*ods.lm[,4] <= 0 &  ods.lm[,2]+ 1.96*ods.lm[,4] >= 0 )/nsim
1-sum(ods.hw[,2]-1.96*ods.hw[,4] <= 0 &  ods.hw[,2]+ 1.96*ods.hw[,4] >= 0 )/nsim
1-sum(ods.hwfi[,2]-1.96*ods.hwfi[,4] <= 0 &  ods.hwfi[,2]+ 1.96*ods.hwfi[,4] >= 0 )/nsim
1-sum(ods.por[,2]-1.96*ods.por[,4] <= 0 &  ods.por[,2]+ 1.96*ods.por[,4] >= 0 )/nsim
1-sum(ods.sor[,2]-1.96*ods.sor[,4] <= 0 &  ods.sor[,2]+ 1.96*ods.sor[,4] >= 0 )/nsim
1-sum(ods.sorfi[,2]-1.96*ods.sorfi[,4] <= 0 &  ods.sorfi[,2]+ 1.96*ods.sorfi[,4] >= 0 )/nsim
1-sum(odsyx.lm[,2]-1.96*odsyx.lm[,4] <= 0 &  odsyx.lm[,2]+ 1.96*odsyx.lm[,4] >= 0 )/nsim
1-sum(odsyx.coss[,2]-1.96*odsyx.coss[,4] <= 0 &  odsyx.coss[,2]+ 1.96*odsyx.coss[,4] >= 0 )/nsim
1-sum(odsyx.por[,2]-1.96*odsyx.por[,4] <= 0 &  odsyx.por[,2]+ 1.96*odsyx.por[,4] >= 0 )/nsim
1-sum(odsyx.sor[,2]-1.96*odsyx.sor[,4] <= 0 &  odsyx.sor[,2]+ 1.96*odsyx.sor[,4] >= 0 )/nsim
1-sum(odsyx.sorfi[,2]-1.96*odsyx.sorfi[,4] <= 0 &  odsyx.sorfi[,2]+ 1.96*odsyx.sorfi[,4] >= 0 )/nsim


1-sum((srs.lm[,1]-1.96*srs.lm[,3] <= 0 &  srs.lm[,1]+ 1.96*srs.lm[,3] >= 0) |( srs.lm[,2]-1.96*srs.lm[,4] <= 0 &  srs.lm[,2]+ 1.96*srs.lm[,4] >= 0) )/nsim
1-sum((srs.hw[,1]-1.96*srs.hw[,3] <= 0 &  srs.hw[,1]+ 1.96*srs.hw[,3] >= 0) | (srs.hw[,2]-1.96*srs.hw[,4] <= 0 &  srs.hw[,2]+ 1.96*srs.hw[,4] >= 0)  )/nsim
1-sum((srs.por[,1]-1.96*srs.por[,3] <= 0 &  srs.por[,1]+ 1.96*srs.por[,3] >= 0) | (srs.por[,2]-1.96*srs.por[,4] <= 0 &  srs.por[,2]+ 1.96*srs.por[,4] >= 0) )/nsim
1-sum((ods.hw[,1]-1.96*ods.hw[,3] <= 0 &  ods.hw[,1]+ 1.96*ods.hw[,3] >= 0) |(ods.hw[,2]-1.96*ods.hw[,4] <= 0 &  ods.hw[,2]+ 1.96*ods.hw[,4] >= 0) )/nsim
1-sum((ods.hwfi[,1]-1.96*ods.hwfi[,3] <= 0 &  ods.hwfi[,1]+ 1.96*ods.hwfi[,3] >= 0) |(ods.hwfi[,2]-1.96*ods.hwfi[,4] <= 0 &  ods.hwfi[,2]+ 1.96*ods.hwfi[,4] >= 0) )/nsim
1-sum((ods.sor[,1]-1.96*ods.sor[,3] <= 0 &  ods.sor[,1]+ 1.96*ods.sor[,3] >= 0) |(ods.sor[,2]-1.96*ods.sor[,4] <= 0 &  ods.sor[,2]+ 1.96*ods.sor[,4] >= 0) )/nsim
1-sum((ods.sorfi[,1]-1.96*ods.sorfi[,3] <= 0 &  ods.sorfi[,1]+ 1.96*ods.sorfi[,3] >= 0) |(ods.sorfi[,2]-1.96*ods.sorfi[,4] <= 0 &  ods.sorfi[,2]+ 1.96*ods.sorfi[,4] >= 0) )/nsim
1-sum((odsyx.coss[,1]-1.96*odsyx.coss[,3] <= 0 &  odsyx.coss[,1]+ 1.96*odsyx.coss[,3] >= 0) |( odsyx.coss[,2]-1.96*odsyx.coss[,4] <= 0 &  odsyx.coss[,2]+ 1.96*odsyx.coss[,4] >= 0 ))/nsim
1-sum((odsyx.sor[,1]-1.96*odsyx.sor[,3] <= 0 &  odsyx.sor[,1]+ 1.96*odsyx.sor[,3] >= 0) |( odsyx.sor[,2]-1.96*odsyx.sor[,4] <= 0 &  odsyx.sor[,2]+ 1.96*odsyx.sor[,4] >= 0 ))/nsim
1-sum((odsyx.sorfi[,1]-1.96*odsyx.sorfi[,3] <= 0 &  odsyx.sorfi[,1]+ 1.96*odsyx.sorfi[,3] >= 0) |( odsyx.sorfi[,2]-1.96*odsyx.sorfi[,4] <= 0 &  odsyx.sorfi[,2]+ 1.96*odsyx.sorfi[,4] >= 0))/nsim

mean(d.srslm)
mean(d.srshw)
mean(d.srspor)
mean(d.srssor)
mean(d.odslm)
mean(d.odshw)
mean(d.odshwfi)
mean(d.odspor)
mean(d.odssor)
mean(d.odssorfi)
mean(d.odsyxpor)
mean(d.odsyxsor)
mean(d.odsyxsorfi)

### Make a forest plot to summarize results
library(forestplot)

mean1 <- c(NA, NA, 0.104, 0.104, 0.103, 0.103, 0.267, 0.101, 0.102, 0.102, 0.102, 0.102,0.033,0.099,0.106,0.104,0.101)
sd1 <- c(NA, NA, 0.127, 0.129, 0.129, 0.129, 0.204, 0.079, 0.079, 0.079, 0.079,0.079, 0.047,0.018,0.079,0.079,0.018)
mean2 <- c(NA, NA, 0.089, 0.035,0.112,0.088,1.136, 0.028,0.031,0.108,0.096,0.096, 0.149,0.033,0.124,0.108,0.097)
sd2 <- c(NA, NA, 0.203,0.079,0.215,0.207,1.775, 0.039,0.046,0.135,0.131,0.131,0.509,0.011,0.129,0.128,0.031)
sim.ess <- 
  structure(list(
    mean  =mean1 , 
    lower =mean1-sd1,
    upper =mean1+sd1),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -11L), 
    class = "data.frame")


tabletext<-cbind(
  c("Sampling", "Scheme", "RS", "",  "", "", "ES-Y", "", "","","","",
    "ES-YX", "", "", "", ""),
  c("Estimation", "Method", "OLS/TNREG", "HW", "POR", "SOR", "OLS/TNREG", "HW", "HW-FI","POR","SOR","SOR-FI",
    "OLS/TNREG", "COSS-FI", "POR","SOR","SOR-FI"))

forestplot(tabletext, legend=c("Normal","Truncated Normal"), 
           fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI), boxsize=.2, line.margin=0.3, 
      mean=cbind(mean1,mean2), lower=cbind(mean1-sd1,mean2-sd2), upper=cbind(mean1+sd1,mean2+sd2),new_page = TRUE,
        is.summary=c(TRUE,TRUE,rep(FALSE,length(mean1)-2)),
           clip=c(-0.2,0.3), xticks=seq(-0.2,0.3,by=0.1),
           xlog=FALSE, zero=0.1, lwd.zero=3, txt_gp=fpTxtGp(cex=1.5),
           col=fpColors(box=c("blue","darkred"), line=c("blue","darkred")), title="Continuous Outcome")



#########################################################################################
############# Factors affecting efficiency and bias ##################################
##########################################################################################


N=1000000
X2 <- rnorm(N) 
X1 <- rbinom(N, size=1, p=0.15) 
Y <- rnorm(N, mean=0.1*X1+0.1*X2)
popdata <-data.frame(X1=X1,X2=X2, Y=Y)


## Varying alpha (cutoff values), fixing psrs (proportion of random sample)
n=500
nsim=500
srs.lm = srs.por = srs.hw=srs.hwfi=matrix(0, nsim,4)
srs.sor = matrix(0, nsim,4)
ods.lm= ods.por= ods.hw=ods.hwfi=matrix(0, nsim,4)
ods.sor=ods.sorfi= matrix(0, nsim, 4)
odsyx.lm= odsyx.por=odsyx.coss= matrix(0, nsim,4)
odsyx.sor= odsyx.sorfi=matrix(0, nsim, 4)
d.srslm=d.srspor=d.srshw=d.srshwfi=d.srssor=d.odslm=d.odspor=d.odshw=d.odshwfi=d.odssor=d.odssorfi=numeric(nsim)
d.odsyxlm=d.odsyxpor=d.odsyxsor=d.odsyxsorfi=d.odsyxcoss=numeric(nsim)
alpha=0.1  ## 0.1, 0.2, 0.3, 0.4, 0.5  
Y10= quantile(Y, alpha); Y90=quantile(Y, 1-alpha)    ## vayring cutoff values for endogenous samples.
psrs=0.25 ## fix the proportion of random samples at 0.25 

for (i in 1:nsim) {
  print(i)
   ## simple random sample with linear regression 
   simdata.srs= popdata[sample(1:N, size=n),]
   sim.reg<-lm(Y~X1+X2, data=simdata.srs)
   srs.lm[i,1:2] <- sim.reg$coef[2:3] ##/sim.reg$coef["sigma"] 
   srs.lm[i,3:4] <- summary(sim.reg)$coef[2:3,2]
   d.srslm[i]= sqrt(det(vcov(sim.reg)[2:3,2:3]))  

   #####################  ODS on y only ########################
   ## ODS with OLS
   simdata = rbind(popdata[sample(1:N, size=n*psrs),], popdata[sample((1:N)[Y<=Y10], size=n*(1-psrs)/2),], popdata[sample((1:N)[ Y>Y90 ], size=n*(1-psrs)/2),])
   sim.reg<-lm(Y~X1+X2, data=simdata)
   ods.lm[i,1:2] <- sim.reg$coef[2:3]  
   ods.lm[i,3:4] <- summary(sim.reg)$coef[2:3,2]


   ## ODS with sor 
   temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2))
   ods.sor[i,1:2]= as.matrix(temp$lambdagamma)[1,(length(temp$lambdagamma)-1):length(temp$lambdagamma)]
   ods.sor[i,3:4]= c(temp$gammase)
   d.odssor[i]= sqrt(det(temp$gammavar))


   #####################  ODS on both y and x ########################
   ## ODS with regular analysis

     xc=0.5 
    simdata = rbind(popdata[sample((1:N)[Y<=Y10 & X1<xc], size=n*(1-psrs)*0.25),], popdata[sample((1:N)[X1<xc], size=n*psrs*0.5),], 
                popdata[sample((1:N)[ Y>Y90 & X1<xc ], size=n*(1-psrs)*0.25),],popdata[sample((1:N)[Y<=Y10 & X1>=xc], size=n*(1-psrs)*0.25),], 
           popdata[sample((1:N)[X1>=xc], size=n*psrs*0.5),], popdata[sample((1:N)[ Y>Y90 & X1>=xc], size=n*(1-psrs)*0.25),])
  if (i==1) {
   weight = c(sum(simdata$Y<=Y10 & simdata$X1<xc)/sum(Y<=Y10 & X1<xc), sum(simdata$Y>Y10 & simdata$Y<=Y90 & simdata$X1<xc)/sum(Y> Y10 & Y<= Y90 & X1<xc), 
             sum(simdata$Y > Y90 & simdata$X1<xc)/sum( Y>Y90 & X1<xc ),
           sum(simdata$Y<=Y10 & simdata$X1>xc)/sum(Y<=Y10 & X1>xc), sum(simdata$Y>Y10 & simdata$Y<=Y90 & simdata$X1>xc)/sum(Y> Y10 & Y<= Y90 & X1>xc), 
             sum(simdata$Y > Y90 & simdata$X1>xc)/sum( Y>Y90 & X1>xc ))*N
   weight= weight/sum(weight, na.rm=T)
             }
   sim.reg<-lm(Y~X1+X2, data=simdata)
   odsyx.lm[i,1:2] <- sim.reg$coef[2:3] 
   odsyx.lm[i,3:4] <- summary(sim.reg)$coef[2:3,2]



  ## ods-yx with SOR +weights

   uy= sort(unique(round(simdata$Y,1)))
   if (psrs==0) {
        wt = matrix(0, nrow=nrow(simdata), ncol=length(uy))
     wt[simdata$X1 < xc, uy <= round(Y10,1) ] = weight[1]
    ## wt[simdata$X1 < xc, uy> Y10 & uy<=Y90 ] = weight[2]
     wt[simdata$X1 < xc, uy >= round(Y90,1)  ] = weight[3]
     wt[simdata$X1 >= xc, uy <= round(Y10,1)  ] = weight[4]
     ##wt[simdata$X1 >= xc, uy> Y10 & uy<=Y90  ] = weight[5]
     wt[simdata$X1 >= xc, uy >=round(Y90,1) ] = weight[6]
   } else 
     {
     wt = matrix(0, nrow=nrow(simdata), ncol=length(uy))
     wt[simdata$X1 < xc, uy <= Y10 ] = weight[1]
     wt[simdata$X1 < xc, uy> Y10 & uy<=Y90 ] = weight[2]
     wt[simdata$X1 < xc, uy > Y90  ] = weight[3]
     wt[simdata$X1 >= xc, uy <= Y10  ] = weight[4]
     wt[simdata$X1 >= xc, uy> Y10 & uy<=Y90  ] = weight[5]
     wt[simdata$X1 >= xc, uy > Y90 ] = weight[6]
     }
   if (any(wt==0)) stop("Weight cannot be zero.")
   ## with known weights
   temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2), weight=wt)
   odsyx.sorfi[i,1:2]= as.matrix(temp$lambdagamma)[1,(length(temp$lambdagamma)-1):length(temp$lambdagamma)]
   odsyx.sorfi[i,3:4]= c(temp$gammase)
   d.odsyxsorfi[i]= sqrt(det(temp$gammavar))
}

niter=500
det(cov(srs.lm[1:niter,1:2]))^0.5/det(cov(ods.sor[1:niter,1:2]))^0.5
det(cov(srs.lm[1:niter,1:2]))^0.5/det(cov(odsyx.sorfi[1:niter,1:2]))^0.5

mean(abs(apply(ods.lm[1:niter,1:2],2,mean)-c(0.1,0.1))/c(0.1,0.1))
mean(abs(apply(odsyx.lm[1:niter,1:2],2,mean)-c(0.1,0.1))/c(0.1,0.1))
mean(abs(apply(ods.sor[1:niter,1:2],2,mean)-c(0.1,0.1))/c(0.1,0.1))
mean(abs(apply(odsyx.sorfi[1:niter,1:2],2,mean)-c(0.1,0.1))/c(0.1,0.1))

## make a plot to summarize results

# results from varying cutoff values
palpha=c(0.5,0.6,0.7,0.8,0.9)
eff.odsy= c(1.00,  1.03, 1.46, 1.71,2.2)
eff.odsyx= c(2.44, 3.36, 4.85, 7.41, 12.4)

olsbias.odsy=c(5.5,22.9, 43.5, 95.3, 165.5)
olsbias.odsyx=c(25.4, 36.7, 50.8, 72.7, 115.5)
sorbias.odsy=c(3.8,4.3, 2.5, 5.5, 1.5)
sorbias.odsyx=c(1.5, 2.1, 3.9, 2.7, 1)

par(mfrow=c(2,2))
plot(palpha, eff.odsyx,xlab=expression(paste("Strata Cutoff Point: ", alpha)), ylab="Sample Info.", type="b", ylim=c(1,13), lwd=2,cex=2, cex.lab=1.5, cex.axis=1.5, pch=1)
lines(palpha, eff.odsy, type="b", lty=2, lwd=2,cex=2,pch=2)
abline(1,0,lty=3, lwd=2,cex=2) 
legend(x=0.5, y=12, legend=c("SOR in ES-YX","SOR in ES-Y", "OLS in RS"), lty=1:3,bty="n",cex=1.5, lwd=2,pch=c(1,2,-1))

plot(palpha, olsbias.odsy,xlab=expression(paste("Strata Cutoff Point: ", alpha)), ylab="Percentage Bias", type="b", ylim=c(0,170), lwd=2,cex=2, cex.lab=1.5, cex.axis=1.5, pch=1)
lines(palpha, olsbias.odsyx, type="b", lty=2, lwd=2,cex=2,pch=2)
lines(palpha, sorbias.odsy, type="l", lty=3, lwd=2,cex=2)
lines(palpha, sorbias.odsyx, type="l", lty=4, lwd=2,cex=2)
##abline(0,0,lty=3) 
legend(x=0.5, y=170, legend=c("OLS in ES-Y","OLS in ES-YX", "SOR in ES-Y", "SOR in ES-YX"), lty=1:4,bty="n",cex=1.5,lwd=2, pch=c(1,2,-1,-1))




## Fixing alpha (cutoff values), varying psrs (proportion of random sample)
n=500
nsim=500
srs.lm = srs.por = srs.hw=srs.hwfi=matrix(0, nsim,4)
srs.sor = matrix(0, nsim,4)
ods.lm= ods.por= ods.hw=ods.hwfi=matrix(0, nsim,4)
ods.sor=ods.sorfi= matrix(0, nsim, 4)
odsyx.lm= odsyx.por=odsyx.coss= matrix(0, nsim,4)
odsyx.sor= odsyx.sorfi=matrix(0, nsim, 4)
d.srslm=d.srspor=d.srshw=d.srshwfi=d.srssor=d.odslm=d.odspor=d.odshw=d.odshwfi=d.odssor=d.odssorfi=numeric(nsim)
d.odsyxlm=d.odsyxpor=d.odsyxsor=d.odsyxsorfi=d.odsyxcoss=numeric(nsim)
alpha=0.1  
Y10= quantile(Y, alpha); Y90=quantile(Y, 1-alpha)    ## vayring cutoff values for endogenous samples.
psrs= 0 ##  0.25, 0.5, 0.75,0.99  ## varying proportions of random samples

for (i in 1:nsim) {
  print(i)
   ## simple random sample with linear regression 
   simdata.srs= popdata[sample(1:N, size=n),]
   sim.reg<-lm(Y~X1+X2, data=simdata.srs)
   srs.lm[i,1:2] <- sim.reg$coef[2:3] ##/sim.reg$coef["sigma"] 
   srs.lm[i,3:4] <- summary(sim.reg)$coef[2:3,2]
   d.srslm[i]= sqrt(det(vcov(sim.reg)[2:3,2:3]))  

   #####################  ODS on y only ########################
   ## ODS with OLS
   simdata = rbind(popdata[sample(1:N, size=n*psrs),], popdata[sample((1:N)[Y<=Y10], size=n*(1-psrs)/2),], popdata[sample((1:N)[ Y>Y90 ], size=n*(1-psrs)/2),])
   sim.reg<-lm(Y~X1+X2, data=simdata)
   ods.lm[i,1:2] <- sim.reg$coef[2:3]  
   ods.lm[i,3:4] <- summary(sim.reg)$coef[2:3,2]


   ## ODS with sor 
   temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2))
   ods.sor[i,1:2]= as.matrix(temp$lambdagamma)[1,(length(temp$lambdagamma)-1):length(temp$lambdagamma)]
   ods.sor[i,3:4]= c(temp$gammase)
   d.odssor[i]= sqrt(det(temp$gammavar))


   #####################  ODS on both y and x ########################
   ## ODS with regular analysis

     xc=0.5 
    simdata = rbind(popdata[sample((1:N)[Y<=Y10 & X1<xc], size=n*(1-psrs)*0.25),], popdata[sample((1:N)[X1<xc], size=n*psrs*0.5),], 
                popdata[sample((1:N)[ Y>Y90 & X1<xc ], size=n*(1-psrs)*0.25),],popdata[sample((1:N)[Y<=Y10 & X1>=xc], size=n*(1-psrs)*0.25),], 
           popdata[sample((1:N)[X1>=xc], size=n*psrs*0.5),], popdata[sample((1:N)[ Y>Y90 & X1>=xc], size=n*(1-psrs)*0.25),])
  if (i==1) {
   weight = c(sum(simdata$Y<=Y10 & simdata$X1<xc)/sum(Y<=Y10 & X1<xc), sum(simdata$Y>Y10 & simdata$Y<=Y90 & simdata$X1<xc)/sum(Y> Y10 & Y<= Y90 & X1<xc), 
             sum(simdata$Y > Y90 & simdata$X1<xc)/sum( Y>Y90 & X1<xc ),
           sum(simdata$Y<=Y10 & simdata$X1>xc)/sum(Y<=Y10 & X1>xc), sum(simdata$Y>Y10 & simdata$Y<=Y90 & simdata$X1>xc)/sum(Y> Y10 & Y<= Y90 & X1>xc), 
             sum(simdata$Y > Y90 & simdata$X1>xc)/sum( Y>Y90 & X1>xc ))*N
   weight= weight/sum(weight, na.rm=T)
             }
   sim.reg<-lm(Y~X1+X2, data=simdata)
   odsyx.lm[i,1:2] <- sim.reg$coef[2:3] 
   odsyx.lm[i,3:4] <- summary(sim.reg)$coef[2:3,2]



  ## ods-yx with SOR +weights

   uy= sort(unique(round(simdata$Y,1)))
   if (psrs==0) {
        wt = matrix(0, nrow=nrow(simdata), ncol=length(uy))
     wt[simdata$X1 < xc, uy <= round(Y10,1) ] = weight[1]
    ## wt[simdata$X1 < xc, uy> Y10 & uy<=Y90 ] = weight[2]
     wt[simdata$X1 < xc, uy >= round(Y90,1)  ] = weight[3]
     wt[simdata$X1 >= xc, uy <= round(Y10,1)  ] = weight[4]
     ##wt[simdata$X1 >= xc, uy> Y10 & uy<=Y90  ] = weight[5]
     wt[simdata$X1 >= xc, uy >=round(Y90,1) ] = weight[6]
   } else 
     {
     wt = matrix(0, nrow=nrow(simdata), ncol=length(uy))
     wt[simdata$X1 < xc, uy <= Y10 ] = weight[1]
     wt[simdata$X1 < xc, uy> Y10 & uy<=Y90 ] = weight[2]
     wt[simdata$X1 < xc, uy > Y90  ] = weight[3]
     wt[simdata$X1 >= xc, uy <= Y10  ] = weight[4]
     wt[simdata$X1 >= xc, uy> Y10 & uy<=Y90  ] = weight[5]
     wt[simdata$X1 >= xc, uy > Y90 ] = weight[6]
     }
   if (any(wt==0)) stop("Weight cannot be zero.")
   ## with known weights
   temp=fitsor(y=round(simdata$Y,1), x=cbind(simdata$X1, simdata$X2), weight=wt)
   odsyx.sorfi[i,1:2]= as.matrix(temp$lambdagamma)[1,(length(temp$lambdagamma)-1):length(temp$lambdagamma)]
   odsyx.sorfi[i,3:4]= c(temp$gammase)
   d.odsyxsorfi[i]= sqrt(det(temp$gammavar))
}

niter=500
det(cov(srs.lm[1:niter,1:2]))^0.5/det(cov(ods.sor[1:niter,1:2]))^0.5
det(cov(srs.lm[1:niter,1:2]))^0.5/det(cov(odsyx.sorfi[1:niter,1:2]))^0.5

mean(abs(apply(ods.lm[1:niter,1:2],2,mean)-c(0.1,0.1))/c(0.1,0.1))
mean(abs(apply(odsyx.lm[1:niter,1:2],2,mean)-c(0.1,0.1))/c(0.1,0.1))
mean(abs(apply(ods.sor[1:niter,1:2],2,mean)-c(0.1,0.1))/c(0.1,0.1))
mean(abs(apply(odsyx.sorfi[1:niter,1:2],2,mean)-c(0.1,0.1))/c(0.1,0.1))

## Results varying the proportion of random samples

psrs=c(0, 0.25,0.5,0.75,1)
eff.odsy= c(2.94, 2.2, 1.96,1.64,0.98)
eff.odsyx= c(14.37, 12.4, 7.51,4.73,2.10)

olsbias.odsy=c(210, 165.5, 110,52.0,5.4)
olsbias.odsyx=c(153, 115.5,78.3,36.8,5.2)
sorbias.odsy=c(2.2, 1.5,3.3,0.8,5.2 )
sorbias.odsyx=c(3.0, 1,2.5,1.9,1.8)


plot(1-psrs, eff.odsyx,xlab="Prop. of Nonrandom Sample: 1-p", ylab="Sample Info.", type="b", ylim=c(1,15), lwd=2,cex=2, cex.lab=1.5, cex.axis=1.5, pch=1)
lines(1-psrs, eff.odsy, type="b", lty=2, lwd=2,cex=2,pch=2)
abline(1,0,lty=3, lwd=2,cex=2) 
legend(x=0, y=13, legend=c("SOR in ES-YX","SOR in ES-Y", "OLS in RS"), lty=1:3,bty="n", lwd=2,cex=1.5, pch=c(1,2,-1))

plot(1-psrs, olsbias.odsy,xlab="Prop. of Nonrandom Sample: 1-p", ylab="Percentage Bias", type="b", ylim=c(0,230), lwd=2,cex=2, cex.lab=1.5, cex.axis=1.5, pch=1)
lines(1-psrs, olsbias.odsyx, type="b", lty=2, lwd=2,cex=2,pch=2)
lines(1-psrs, sorbias.odsy, type="l", lty=3, lwd=2,cex=2)
lines(1-psrs, sorbias.odsyx, type="l", lty=4, lwd=2,cex=2)
##abline(0,0,lty=3) 
legend(x=0, y=230, legend=c("OLS in ES-Y","OLS in ES-YX", "SOR in ES-Y", "SOR in ES-YX"), lwd=2,lty=1:4,bty="n",cex=1.5, pch=c(1,2,-1,-1))






#####################################################################################################
#### Part 2: Code to compare OLS, Heckit and SOR for selective sampling that results in missing outcomes
#####################################################################################################


library(sampleSelection)
library(mvtnorm)

nsim=500
ols.res<-matrix(0,nsim,2)
heck.res<-matrix(0,nsim,2)
sor.res <-matrix(0, nsim,2)
n=10000

## 1. linear probit selection probability
for (i in 1:nsim) {
  print(i)

  x2<- rnorm(n)
  x1 <- rbinom(n, size=1, p=0.15) 
  y<-0.1*x1+ 0.1*x2+ rnorm(n) 

  g=rbinom(n, 1, pnorm(y))
  yo <- ifelse(g==1, y, NA)
  ols.res[i,]<- coef(lm(yo~x1+x2))[-1]
  heck.res[i,] <-summary(selection(g~x1+x2,yo~x1+x2, method="mle"))$estimate[5:6,1]
  temp=fitsor(y=round(yo[!is.na(yo)],1), x=cbind(x1[!is.na(yo)],x2[!is.na(yo)] ))
  sor.res[i,]= as.matrix(temp$lambdagamma)[1, (length(temp$lambdagamma)-1):length(temp$lambdagamma)]
}

apply(ols.res,2,mean)
apply(heck.res,2,mean)
apply(sor.res,2,mean)
apply(ols.res,2,sd)
apply(heck.res,2,sd)
apply(sor.res,2,sd)
 
## 2. quadratic probit selection probability
for (i in 1:nsim) {
  print(i)

  x2<- rnorm(n)
  x1 <- rbinom(n, size=1, p=0.15) 
  y<-0.1*x1+ 0.1*x2+ rnorm(n) 

  g=rbinom(n, 1, pnorm(-1+10*0.1*y+10*0.1*y^2))
  yo <- ifelse(g==1, y, NA)
  ols.res[i,]<- coef(lm(yo~x1+x2))[-1]
  heck.res[i,] <-summary(selection(g~x1+x2,yo~x1+x2, method="mle"))$estimate[5:6,1]
  temp=fitsor(y=round(yo[!is.na(yo)],1), x=cbind(x1[!is.na(yo)],x2[!is.na(yo)] ))
  sor.res[i,]= as.matrix(temp$lambdagamma)[1, (length(temp$lambdagamma)-1):length(temp$lambdagamma)]
}

apply(ols.res,2,mean)
apply(heck.res,2,mean)
apply(sor.res,2,mean)
apply(ols.res,2,sd)
apply(heck.res,2,sd)
apply(sor.res,2,sd)

## 3. Extreme value selected
for (i in 1:nsim) {
  print(i)

  x2<- rnorm(n)
  x1 <- rbinom(n, size=1, p=0.15) 
  y<-0.1*x1+ 0.1*x2+ rnorm(n)

  ## Sort by y
  x1<-x1[order(y)]; x2<-x2[order(y)]
  y<-y[order(y)]
  yo <- y

  L=1000; R=9000
  yo[L:R] <-NA
  g = ifelse(is.na(yo),0,1)
  ols.res[i,]<- coef(lm(yo~x1+x2))[-1]
  heck.res[i,] <-summary(selection(g~x1+x2,yo~x1+x2, method="2step"))$estimate[5:6,1] ## mle does not work for some data, use 2step procedure instead. 
  temp=fitsor(y=round(yo[!is.na(yo)],1), x=cbind(x1[!is.na(yo)],x2[!is.na(yo)] ))
  sor.res[i,]= as.matrix(temp$lambdagamma)[1, (length(temp$lambdagamma)-1):length(temp$lambdagamma)]
}

apply(ols.res,2,mean)
apply(heck.res,2,mean)
apply(sor.res,2,mean)
apply(ols.res,2,sd)
apply(heck.res,2,sd)
apply(sor.res,2,sd)


## 4. Middle value selected
for (i in 1:nsim) {
  print(i)

  x2<- rnorm(n)
  x1 <- rbinom(n, size=1, p=0.15) 
  y<-0.1*x1+ 0.1*x2+ rnorm(n)

  ## Sort by y
  x1<-x1[order(y)]; x2<-x2[order(y)]
  y<-y[order(y)]
  yo <- y

  L=2500; R=7500
  yo[1:L] <- NA
  yo[R:n] <- NA
  g = ifelse(is.na(yo),0,1)
  ols.res[i,]<- coef(lm(yo~x1+x2))[-1]
  heck.res[i,] <-summary(selection(g~x1+x2,yo~x1+x2, method="mle"))$estimate[5:6,1]
  temp=fitsor(y=round(yo[!is.na(yo)],1), x=cbind(x1[!is.na(yo)],x2[!is.na(yo)] ))
  sor.res[i,]= as.matrix(temp$lambdagamma)[1, (length(temp$lambdagamma)-1):length(temp$lambdagamma)]
}

apply(ols.res,2,mean)
apply(heck.res,2,mean)
apply(sor.res,2,mean)
apply(ols.res,2,sd)
apply(heck.res,2,sd)
apply(sor.res,2,sd)



## 5. Linear probit selection stratified on X2
for (i in 1:nsim) {
  print(i)
  x2<- rnorm(n)
  x1 <- rbinom(n, size=1, p=0.15) 
  y<-0.1*x1+ 0.1*x2+ rnorm(n) 

  g=rbinom(n, 1, pnorm(y-2*y*as.numeric(x2>0)))
  yo <- ifelse(g==1, y, NA)
  ols.res[i,]<- coef(lm(yo~x1+x2))[-1]
  heck.res[i,] <-summary(selection(g~x1+x2,yo~x1+x2, method="2step"))$estimate[5:6,1]
  temp=fitsor(y=round(yo[!is.na(yo)],1), x=cbind(x1[!is.na(yo)],x2[!is.na(yo)] ), strata=ifelse(x2[!is.na(yo)]<0,1,2) )
  sor.res[i,]= as.matrix(temp$lambdagamma)[1, (length(temp$lambdagamma)-1):length(temp$lambdagamma)]
}

apply(ols.res,2,mean)
apply(heck.res,2,mean)
apply(sor.res,2,mean)
apply(ols.res,2,sd)
apply(heck.res,2,sd)
apply(sor.res,2,sd)


## 6. Quadratic probit selection stratified on X2
for (i in 1:nsim) {
  print(i)
  x2<- rnorm(n)
  x1 <- rbinom(n, size=1, p=0.15) 
  y<-0.1*x1+ 0.1*x2+ rnorm(n) 

  g=rbinom(n, 1, pnorm(-1+y+y^2+2*as.numeric(x2>0)-2*y*as.numeric(x2>0)-2*y^2*as.numeric(x2>0)))
  yo <- ifelse(g==1, y, NA)
  ols.res[i,]<- coef(lm(yo~x1+x2))[-1]
  heck.res[i,] <-summary(selection(g~x1+x2,yo~x1+x2, method="mle"))$estimate[5:6,1]
  temp=fitsor(y=round(yo[!is.na(yo)],1), x=cbind(x1[!is.na(yo)],x2[!is.na(yo)] ), strata=ifelse(x2[!is.na(yo)]<0,1,2) )
  sor.res[i,]= as.matrix(temp$lambdagamma)[1, (length(temp$lambdagamma)-1):length(temp$lambdagamma)]
}

apply(ols.res,2,mean)
apply(heck.res,2,mean)
apply(sor.res,2,mean)
apply(ols.res,2,sd)
apply(heck.res,2,sd)
apply(sor.res,2,sd)



## 7. Extreme value selected stratified on x2
for (i in 1:nsim) {
  print(i)

  x2<- rnorm(n)
  x1 <- rbinom(n, size=1, p=0.15) 
  y<-0.1*x1+ 0.1*x2+ rnorm(n)

  ## Sort by y
  x1<-x1[order(y)]; x2<-x2[order(y)]
  y<-y[order(y)]
  yo <- y

  yo[x2<0][500:4500] <- NA
  yo[x2>0][1000:4000]<- NA
  g = ifelse(is.na(yo),0,1)
  ols.res[i,]<- coef(lm(yo~x1+x2))[-1]
  heck.res[i,] <-summary(selection(g~x1+x2,yo~x1+x2, method="mle"))$estimate[5:6,1]
  temp=fitsor(y=round(yo[!is.na(yo)],1), x=cbind(x1[!is.na(yo)],x2[!is.na(yo)] ), strata=ifelse(x2[!is.na(yo)]<0,1,2) )
  sor.res[i,]= as.matrix(temp$lambdagamma)[1, (length(temp$lambdagamma)-1):length(temp$lambdagamma)]
}

apply(ols.res,2,mean)
apply(heck.res,2,mean)
apply(sor.res,2,mean)
apply(ols.res,2,sd)
apply(heck.res,2,sd)
apply(sor.res,2,sd)


## 8. Middle value selected stratified on X2
for (i in 1:nsim) {
  print(i)

  x2<- rnorm(n)
  x1 <- rbinom(n, size=1, p=0.15) 
  y<-0.1*x1+ 0.1*x2+ rnorm(n)

  ## Sort by y
  x1<-x1[order(y)]; x2<-x2[order(y)]
  y<-y[order(y)]
  yo <- y

  yo[x2<0][1:1500]<-NA; yo[x2<0][3500:sum(x2<0)] <- NA
  yo[x2>0][1:1000]<-NA; yo[x2>0][4000:sum(x2>0)]<- NA

  g = ifelse(is.na(yo),0,1)
  ols.res[i,]<- coef(lm(yo~x1+x2))[-1]
  heck.res[i,] <-summary(selection(g~x1+x2,yo~x1+x2, method="mle"))$estimate[5:6,1]
  temp=fitsor(y=round(yo[!is.na(yo)],1), x=cbind(x1[!is.na(yo)],x2[!is.na(yo)] ), strata=ifelse(x2[!is.na(yo)]<0,1,2) )
  sor.res[i,]= as.matrix(temp$lambdagamma)[1, (length(temp$lambdagamma)-1):length(temp$lambdagamma)]
}

apply(ols.res,2,mean)
apply(heck.res,2,mean)
apply(sor.res,2,mean)
apply(ols.res,2,sd)
apply(heck.res,2,sd)
apply(sor.res,2,sd)



#####################################################################################################
#### Part 3: Code to compare methods on Truncated and on-site samples. 
#####################################################################################################
library(truncnorm)
library(truncreg)
library(extraDistr)
library(VGAM)
library(sandwich)
library(foreign)
library(MASS)
library(haven)
library(countreg)


######################## 3.1. Fitting Store Shopper Data with different models ################################

#################  START of Stata Code #####################################################
## Parametric count models (Poisson regression, truncated poisson regression, negative binomial regression,
##   and truncated negative binomial regression) are fitted using Stata, in order to obtain robust standard errors. 
## Before running the Stata code below, first import "visits.csv" into Stata. 
########################################################################


poisson totalvisit incomec distc distc2 agec  marriedc kidsc
poisson totalvisit incomec distc distc2 agec  marriedc kidsc,  vce(robust)


tpoisson totalvisit incomec distc distc2 agec  marriedc kidsc, ll(0) 
tpoisson totalvisit incomec distc distc2 agec  marriedc kidsc, ll(0) vce(robust)


nbreg totalvisit incomec distc distc2 agec  marriedc kidsc 
nbreg totalvisit incomec distc distc2 agec  marriedc kidsc,  vce(robust)


tnbreg totalvisit incomec distc distc2 agec  marriedc kidsc, ll(0)
tnbreg totalvisit incomec distc distc2 agec  marriedc kidsc, ll(0) vce(robust)

#################### END of Stata code  #############################################################################

retail<- read.csv("visits.csv", header = TRUE)

### analysis of store shopper data using SOR model with log-bilinear form of OR function. 
retail$totalvisitc= (retail$totalvisit-mean(retail$totalvisit))/sd(retail$totalvisit) ## For stable computation of SORE model. 
retail.sor=fitsor(y=retail$totalvisitc, x=retail[,c("incomec", "distc","distc2", "agec","MARRIEDc","kidsc")])
sor.res<-cbind(as.matrix(retail.sor$lambdagamma)[1, (length(retail.sor$lambdagamma)-5):length(retail.sor$lambdagamma) ],
        retail.sor$gammase,
     as.matrix(retail.sor$lambdagamma)[1, (length(retail.sor$lambdagamma)-5):length(retail.sor$lambdagamma) ]/retail.sor$gammase)
sor.res
## The code below gives Intercept parameter estimate. The standard error estimate can be obtained using delta method or the resampling method as described in resampling experiment below. 
p= c(as.matrix(retail.sor$lambdagamma)[1,1:(length(retail.sor$lambdagamma)-6)],0)
p=exp(p)/sum(exp(p))
log(sum(p*sort(unique(retail$totalvisit)))) 

## SOR regression with log-linear-soft-plus of OR function

nlg=length(retail.sor$lambdagamma)
lambda=as.matrix(retail.sor$lambdagamma)[1,1:(nlg-6)]
gamma = as.matrix(retail.sor$lambdagamma)[1,(nlg-5):nlg]
retail.Rsor= Rfitsor(totalvisitc~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retail, par0=c(lambda,gamma,10), hessian=T)
gammase=sqrt(diag(solve(retail.Rsor$hess)))[(length(retail.sor$lambdagamma)-5):length(retail.sor$lambdagamma) ]
sorLLS.res<-cbind(retail.Rsor$par[(length(retail.sor$lambdagamma)-5):length(retail.sor$lambdagamma) ],
        gammase,retail.Rsor$par[(length(retail.sor$lambdagamma)-5):length(retail.sor$lambdagamma) ]/gammase)
sorLLS.res
p= c(as.matrix(retail.Rsor$par)[1:(length(retail.Rsor$par)-7)],0)
p=exp(p)/sum(exp(p))
log(sum(p*sort(unique(retail$totalvisit))))


###################### Comparing Models Fits ##################################

uy=sort(unique(retail$totalvisit))
## Fitted values using Poisson Regression
retail.pois=glm(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retail, family=poisson)
mean(fitted(retail.pois))
ols.freq<-matrix(0, nrow=nrow(retail), ncol=length(uy))
for (i in 1:nrow(retail)) ols.freq[i,]= dpois(uy, lambda=retail.pois$fit[i])*nrow(retail) 

## Fitted values using Truncated Poisson Regression
retail.Tpois<-vglm(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retail,family=pospoisson())
tp.freq= matrix(0, nrow=nrow(retail), ncol=length(uy))
for (i in 1:nrow(retail)) tp.freq[i,]= dtpois(uy, lambda=fitted.values(retail.Tpois)[i,1],a=0)*nrow(retail) 

## Fitted values using negative binomail regression

retail.nb<-glm.nb(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retail)
nb.freq= matrix(0, nrow=nrow(retail), ncol=length(uy))
for (i in 1:nrow(retail)) nb.freq[i,]= dnbinom(uy, size=retail.nb$theta, mu=retail.nb$fitted[i])*nrow(retail) 

## Fitted values using truncated nb regression.
par0<-c(-21.313, 0.049,-0.609,0.152,-0.002,0.133,-0.095,21.968)
retail.tnb<-fun.fitnb(ymodel=totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retail, par0=par0, dist="Tnb")
lmu<-c(cbind(1, as.matrix(retail[,c("incomec","distc","distc2","agec","MARRIEDc","kidsc")]))%*%retail.tnb$par[-length(retail.tnb$par)])
res<-matrix(0, nrow=nrow(retail), ncol=length(uy))
for (i in 1:nrow(res)) {
  for (j in 1:ncol(res)) res[i,j]<-exp(fun.dTnb(uy[j], retail.tnb$par[length(retail.tnb$par)], lmu[i]))
}
tnb.freq= apply(res,2,mean)*nrow(retail)

## fitted values using SOR model with LLL OR function. 
retailvisit.sor=fitsor(y=retail$totalvisitc, x=retail[,c("incomec", "distc","distc2", "agec","MARRIEDc","kidsc")])

nlg=length(retailvisit.sor$lambdagamma)
sor.freq=pred.sor(uy=sort(unique(retail$totalvisitc)), Xgrid=retail[,c("incomec", "distc","distc2", "agec","MARRIEDc","kidsc")], 
  lambda=as.matrix(retailvisit.sor$lambdagamma)[1,-((nlg-5):nlg)],
  gamma=as.matrix(retailvisit.sor$lambdagamma)[1, (nlg-5):nlg])*nrow(retail)

## fitted values using the SOR with LLS odds ratio.


uy=sort(unique(retail$totalvisit))
p=c(retail.Rsor$par[1:(nlg-6)],0)
  p=exp(p)/sum(exp(p))
  p=p/sum(p)
   int=sum(uy*p) 
sorNB.freq=pred.sorNB(uy=sort(unique(retail$totalvisitc)), Xgrid=retail[,c("incomec", "distc","distc2", "agec","MARRIEDc","kidsc")], 
  lambda=(retail.Rsor$par)[-((nlg-5):(nlg+1))],
  gamma=(retail.Rsor$par)[(nlg-5):(nlg)], alpha=exp(retail.Rsor$par[length(retail.Rsor$par)]))*nrow(retail)


barplot(table(retail$totalvisit), width=1, xlab="Number of Store Visits", ylab="Number of Customers", space=0)
lines(uy-0.5, apply(ols.freq,2,mean), type="b",pch=5, lwd=2, lty=5)
lines(uy-0.5, apply(tp.freq,2,mean), type="b",pch=4, lwd=2,lty=4)
lines(uy-0.5, apply(nb.freq,2,mean), type="b",pch=3, lwd=2,lty=3)
##lines(uy-0.5,prob.tnb*nrow(retail), type="b",pch=2, lwd=2,lty=2)
lines(uy-0.5, tnb.freq, type="b",pch=2, lwd=2,lty=2)
lines(uy-0.5, apply(sor.freq,2,mean), type="b", pch=1, lwd=2)
legend(10,2000, legend=c("SOR Regression", "Zero-truncated NB Regression", "Negative Binomial (NB) Regression","Zero-Truncated Poisson Regression", "Poisson Regression"), pch=1:5, lwd=2, lty=1:5)




######################################################################################
###  3.2  resampling experiments to compare performance of different methods
#####################################################################################


#### evaluation assuming truncated poisson distribution

nsim=500
tpois.srs=tpc.srs=sor.srs=tpois.ods=tpc.ods=sor.ods=matrix(0,nsim,7); tnb.srs=tnbc.srs=tnb.ods=tnbc.ods=matrix(0,nsim,8)
llk.tpoissrs=llk.tpcsrs=llk.tnbsrs=llk.tnbcsrs=llk.sorsrs=llk.tpoisods=llk.tpcods=llk.tnbods=llk.tnbcods=llk.sorods=numeric(nsim)
d.tpoissrs=d.tpcsrs=d.tnbsrs=d.sorsrs=d.tpoisods=d.tpcods=d.sorods=numeric(nsim)
## retail.Tpois<-vglm(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retail,family=pospoisson())
par0= c(0.223787835,  0.043471716, -0.515294397,  0.137664458, -0.006419527,  0.106714292, -0.056811605)
eta=c(as.matrix(cbind(1,retail[,c("incomec", "distc", "distc2", "agec", "MARRIEDc","kidsc")]))%*%as.matrix(par0))

for (i in 1:nsim){
  print(i)

  ### random sampling 
  retailtp=retail
  retailtp$totalvisit=rtpois(nrow(retail), lambda=exp(eta), a=0)   ##exp(predictors(retail.Tpois)),a=0)
  retailrs.Tpois<-vglm(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp,family=pospoisson())
  tpois.srs[i,]=c(coef(retailrs.Tpois)) ##, coef(summary(retailrs.Tpois))[,2])
  llk.tpoissrs[i]=-logLik(retailrs.Tpois)

  retail.tpc<-glm(totalvisit-1~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp,family=poisson)
  tpc.srs[i,]=c(coef(retail.tpc)) ## ,coef(summary(retail.tpc))[,2])
  llk.tpcsrs[i]=-logLik(retail.tpc)

  retail.tnb<-fun.fitnb(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp, dist="Tnb", par0=c(retail.tpc$coef, -10))
  tnb.srs[i,]=c(retail.tnb$par) ##, sqrt(diag(estvcov))) 
    llk.tnbsrs[i]=retail.tnb$value

  retailrs.tnbc<-fun.fitnb(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp, par0=c(retail.tpc$coef, -10), dist="OSnb") 
  tnbc.srs[i,]=c(retailrs.tnbc$par) ##, sqrt(diag(estvcov))) 
  llk.tnbcsrs[i]=retailrs.tnbc$value

  retailrs.sor=fitsor(y=retailtp$totalvisit, x=retailtp[,c("incomec", "distc", "distc2", "agec", "MARRIEDc","kidsc")])
  nlg=length(retailrs.sor$lambdagamma)
   uy= sort(unique(retailtp$totalvisit))
  p=c(as.matrix(retailrs.sor$lambdagamma)[1,1:(nlg-6)],0)
  p=exp(p)/sum(exp(p))
  p0= exp(-exp(coef(retailrs.Tpois)[1]))
   int=sum(uy*p)*(1-p0)
  sor.srs[i,1:7]= c(log(int), as.matrix(retailrs.sor$lambdagamma)[1,(nlg-5):nlg])  
  llk.sorsrs[i]=retailrs.sor$negllk
}

apply(tpois.srs,2,mean)
apply(tpc.srs,2,mean)
apply(tnb.srs,2,mean)
apply(tnbc.srs,2,mean)
apply(sor.srs,2,mean)

apply(tpois.srs,2,sd)
apply(tpc.srs,2,sd)
apply(tnb.srs,2,sd)
apply(tnbc.srs,2,sd)
apply(sor.srs,2,sd)

det(cov(tpois.srs[,2:7]))^(1/6)
det(cov(tnb.srs[,2:7]))^(1/6)
det(cov(sor.srs[,2:7]))^(1/6)



for (i in 1:nsim){
  print(i)

  ## On-site sampling
  retailtp$totalvisit=rpois(nrow(retail), lambda=exp(eta))+1
  retailrs.Tpois<-vglm(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp,family=pospoisson())
  tpois.ods[i,]=c(coef(retailrs.Tpois)) 
  llk.tpoisods[i]=-logLik(retailrs.Tpois)

  retail.tpc<-glm(totalvisit-1~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp,family=poisson)
  tpc.ods[i,]=c(coef(retail.tpc)) 
  llk.tpcods[i]=-logLik(retail.tpc)

   retail.tnb<-fun.fitnb(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp, dist="Tnb", par0=c(retail.tpc$coef, -10))
  ##estvcov= solve(retail.tnb$hessian)
  tnb.ods[i,]=c(retail.tnb$par) 
  llk.tnbods[i]=retail.tnb$value

  retailrs.tnbc<-fun.fitnb(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp, par0=c(retail.tpc$coef, -10), dist="OSnb") 
  tnbc.ods[i,]=c(retailrs.tnbc$par) 
    llk.tnbcods[i]=retailrs.tnbc$value

  retailrs.sor=fitsor(y=retailtp$totalvisit, x=retailtp[,c("incomec", "distc","distc2", "agec","MARRIEDc","kidsc")])
  nlg=length(retailrs.sor$lambdagamma)
  uy= sort(unique(retailtp$totalvisit))
  p=c(as.matrix(retailrs.sor$lambdagamma)[1,1:(nlg-6)],0)
  p=exp(p)/sum(exp(p))
  p=p/uy
  p=p/sum(p)
   p0= exp(-exp(tpc.ods[i]))
   int=sum(uy*p)*(1-p0)
  sor.ods[i,1:7]= c(log(int), as.matrix(retailrs.sor$lambdagamma)[1,(nlg-5):nlg]) ##, retailrs.sor$gammase)
   llk.sorods[i]=retailrs.sor$negllk
}


apply(tpois.ods,2,mean)
apply(tpc.ods,2,mean)
apply(tnb.ods,2,mean)
apply(tnbc.ods,2,mean)
apply(sor.ods,2,mean)

apply(tpois.ods,2,sd)
apply(tpc.ods,2,sd)
apply(tnb.ods,2,sd)
apply(tnbc.ods,2,sd)
apply(sor.ods,2,sd)

det(cov(tpc.ods[,2:7]))^(1/6)
det(cov(tnbc.ods[,2:7]))^(1/6)
det(cov(sor.ods[,2:7]))^(1/6)



#### evaluation assuming truncated NB distribution

nsim=500
tpois.srs=tpc.srs=tpois.ods=tpc.ods=matrix(0,nsim,7); tnb.srs=tnbc.srs=tnb.ods=tnbc.ods=sor.srs=sor.ods=matrix(0,nsim,8)
d.tpoissrs=d.tpcsrs=d.tnbsrs=d.tnbcsrs=d.sorsrs=d.tpoisods=d.tpcods=d.tnbods=d.tnbcods=d.sorods=numeric(nsim)
llk.tpoissrs=llk.tpcsrs=llk.tnbsrs=llk.tnbcsrs=llk.sorsrs=llk.tpoisods=llk.tpcods=llk.tnbods=llk.tnbcods=llk.sorods=nlambda=numeric(nsim)

par0<- c(-21.313, 0.049,-0.609,0.152,-0.002,0.133,-0.095,21.968) 
lmu<-c(cbind(1, as.matrix(retail[,c("incomec","distc","distc2","agec","MARRIEDc","kidsc")]))%*%par0[-length(par0)]) 
lalpha<- par0[length(par0)]   

for (i in 1:nsim){
  print(i)

  ### random sampling 
  retailtp=retail
  retailtp$totalvisit=rztnbinom(nrow(retail), mu=exp(lmu), theta=1/exp(lalpha))
  retailrs.Tpois<-vglm(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp,family=pospoisson(), coefstart=c(0,par0[2:7]))
  tpois.srs[i,]=c(coef(retailrs.Tpois)) ## , coef(summary(retailrs.Tpois))[,2])
  llk.tpoissrs[i]=-logLik(retailrs.Tpois)

  retail.tpc<-glm(totalvisit-1~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp,family=poisson)
  tpc.srs[i,]=c(coef(retail.tpc)) ## ,coef(summary(retail.tpc))[,2])
  llk.tpcsrs[i]=-logLik(retail.tpc)

   
  retailrs.tnb<-fun.fitnb(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp, par0=par0, dist="Tnb") 
  tnb.srs[i,]=c(retailrs.tnb$par) ## , sqrt(diag(estvcov))) 
 llk.tnbsrs[i]=retailrs.tnb$value

  retailrs.tnbc<-fun.fitnb(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp, par0=par0, dist="OSnb") 
  tnbc.srs[i,]=c(retailrs.tnbc$par) ## , sqrt(diag(estvcov))) 
  llk.tnbcsrs[i]=retailrs.tnbc$value

   retailtp$totalvisitc= (retailtp$totalvisit-mean(retailtp$totalvisit))/sd(retailtp$totalvisit)
  retailrs.sor=fitsor(y=retailtp$totalvisitc, x=retailtp[,c("incomec", "distc", "distc2", "agec", "MARRIEDc","kidsc")])
  nlg=length(retailrs.sor$lambdagamma)
   uy= sort(unique(retailtp$totalvisit))
  lambda=as.matrix(retailrs.sor$lambdagamma)[1,1:(nlg-6)]
  gamma = as.matrix(retailrs.sor$lambdagamma)[1,(nlg-5):nlg]
  retailrs.sor= Rfitsor(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp, par0=c(lambda,par0[2:7],0.578), hessian=F)
  p=c(retailrs.sor$par[1:(nlg-6)],0)
  p=exp(p)/sum(exp(p))
  p0= exp(fun.dnb(0, lalpha=retailrs.tnb$par[length(retailrs.tnb$par)], lmu=retailrs.tnb$par[1]))
   int=sum(uy*p)*(1-p0)
  sor.srs[i,1:8]= c(log(int), retailrs.sor$par[(length(retailrs.sor$par)-6):length(retailrs.sor$par)]) 
    llk.sorsrs[i]=retailrs.sor$value
}


apply(tpois.srs,2,mean)
apply(tpc.srs,2,mean)
apply(tnb.srs,2,mean)
apply(tnbc.srs,2,mean)
apply(sor.srs,2,mean)
mean(sor.srs[,8]-sor.srs[,1]) ## for log(alpha)

apply(tpois.srs,2,sd)
apply(tpc.srs,2,sd)
apply(tnb.srs,2,sd)
apply(tnbc.srs,2,sd)
apply(sor.srs,2,sd)
sd(sor.srs[,8]-sor.srs[,1])

det(cov(tnb.srs[,2:7]))^(1/6)
det(cov(sor.srs[,2:7]))^(1/6)


for (i in 1:nsim){
  print(i)

  ## onsite sampling
  uuy=1:100; pmat=matrix(0, nrow=nrow(retail), ncol=length(uuy))
  for (j in 1:nrow(retail)) {
    for (k in 1:length(uuy)) pmat[j,k]=exp(fun.dOSnb(uuy[k], lalpha=lalpha, lmu=lmu[j]))
    }
  pmat=t(apply(pmat, 1, function(i) i/sum(uuy)))
  retailtp=retail
  retailtp$totalvisit=apply(pmat, 1,  function(i) uuy[as.vector(rmultinom(1,size=1, i))==1])
  retailrs.Tpois<-vglm(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp,family=pospoisson())
  tpois.ods[i,]=c(coef(retailrs.Tpois)) ##,coef(summary(retailrs.Tpois))[,2])
  llk.tpoisods[i]=-logLik(retailrs.Tpois)

  retail.tpc<-glm(totalvisit-1~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp,family=poisson)
  tpc.ods[i,]=c(coef(retail.tpc)) ##,coef(summary(retail.tpc))[,2])
  llk.tpcods[i]=-logLik(retail.tpc)

   par0<-c(-21.313, 0.049,-0.609,0.152,-0.002,0.133,-0.095,21.968)
  retailrs.tnb<-fun.fitnb(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp, par0=par0, dist="Tnb") 
  tnb.ods[i,]=c(retailrs.tnb$par) 
  llk.tnbods[i]=retailrs.tnb$value

   par0<-c(-21.313, 0.049,-0.609,0.152,-0.002,0.133,-0.095,21.968)
  retailrs.tnbc<-fun.fitnb(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp, par0=par0, dist="OSnb") 
  tnbc.ods[i,]=c(retailrs.tnbc$par) 
  llk.tnbcods[i]=retailrs.tnbc$value


  retailtp$totalvisitc= (retailtp$totalvisit-mean(retailtp$totalvisit))/sd(retailtp$totalvisit)
  retailrs.sor=fitsor(y=retailtp$totalvisitc, x=retailtp[,c("incomec", "distc", "distc2", "agec", "MARRIEDc","kidsc")])
  nlg=length(retailrs.sor$lambdagamma)
   uy= sort(unique(retailtp$totalvisit))
  lambda=as.matrix(retailrs.sor$lambdagamma)[1,1:(nlg-6)]
  gamma = as.matrix(retailrs.sor$lambdagamma)[1,(nlg-5):nlg]
  retailrs.sor= Rfitsor(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailtp, par0=c(lambda,par0[2:7],0.578), hessian=F)
  p=c(retailrs.sor$par[1:(nlg-6)],0)
  p=exp(p)/sum(exp(p))
  p=p/uy
  p=p/sum(p)
  p0= exp(fun.dnb(0, lalpha=retailrs.tnb$par[length(retailrs.tnbc$par)], lmu=retailrs.tnbc$par[1]))  ## To make the intercept parameter comparable 
   int=sum(uy*p)*(1-p0)
  sor.ods[i,1:8]= c(log(int), retailrs.sor$par[(length(retailrs.sor$par)-6):length(retailrs.sor$par)]) 
  llk.sorods[i]=retailrs.sor$value
  nlambda[i]= length(lambda)
}

apply(tpois.ods,2,mean)
apply(tpc.ods,2,mean)
apply(tnb.ods,2,mean)
apply(tnbc.ods,2,mean)
apply(sor.ods,2,mean)
mean(-sor.ods[,1]+sor.ods[,8])

apply(tpois.ods,2,sd)
apply(tpc.ods,2,sd)
apply(tnb.ods,2,sd)
apply(tnbc.ods,2,sd)
apply(sor.ods,2,sd)
sd(-sor.ods[,1]+sor.ods[,8])

det(cov(tnbc.ods[,2:7]))^(1/6)
det(cov(sor.ods[,2:7]))^(1/6)




## Evaluation assuming SOR distribution

## True value
## Using the gamma (log-odds ratio ) estimates for SOR reported in Table 4 and 
## the lambda estimates for baseline function and unique values of visits  reported in Online Appendix I.
gamma <- c(0.012898, -0.317767,  0.073386, -0.002332,  0.070468, -0.023261) 
lambda <-  c(22.84205, 21.32392, 20.38286, 19.67065, 19.37290, 18.68119, 18.33266, 17.99773, 17.54111, 17.33210, 16.58061, 16.92434,
  16.74171,  15.28760, 15.62247, 16.24222, 15.76093, 14.98963, 14.13334, 13.78327, 14.29084, 13.40592, 13.10757, 13.00548, 12.69196,
   12.36838) 
meany=  1.983182; sy=3.412023      
uyc <- (c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,24,26,28,31,32,35,38, 131)-meany)/sy ## The uniquely observed values of standardized store visits.
pmat=pred.sor(uy=uyc, Xgrid=retail[,c("incomec", "distc","distc2", "agec","MARRIEDc","kidsc")], lambda=lambda, gamma=gamma)
pmat=t(apply(pmat, 1, function(i) i/sum(i)))

nsim=500
tpois.srs=tpc.srs=tpois.ods=tpc.ods=matrix(0,nsim,7); tnb.srs=tnbc.srs=tnb.ods=tnbc.ods=sor.srs=sor.ods=matrix(0,nsim,8)
d.tpoissrs=d.tpcsrs=d.tnbsrs=d.tnbcsrs=d.sorsrs=d.tpoisods=d.tpcods=d.tnbods=d.tnbcods=d.sorods=numeric(nsim)
llk.tpoissrs=llk.tpcsrs=llk.tnbsrs=llk.tnbcsrs=llk.sorsrs=llk.tpoisods=llk.tpcods=llk.tnbods=llk.tnbcods=llk.sorods=numeric(nsim)

for (i in 1:nsim){
   print(i)
  ## random sampling 
  retailsor=retail
  retailsor$totalvisit=apply(pmat, 1,  function(i) c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,24,26,28,31,32,35,38, 131)[as.vector(rmultinom(1,size=1, i))==1])
  retailsor$totalvisitc= (retailsor$totalvisit-meany)/sy

  ## retailsor$totalvisit=round(retailsor$totalvisitc*sd(retail$totalvisit)+mean(retail$totalvisit),0)
  ##retailsor$totalvisit=round(retailsor$totalvisitc*3.412023+1.983182,0)
  retailrs.Tpois<-vglm(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailsor,family=pospoisson())
  tpois.srs[i,]=c(coef(retailrs.Tpois)) ## , coef(summary(retailrs.Tpois))[,2])
  llk.tpoissrs[i]=-logLik(retailrs.Tpois)

  retail.tpc<-glm(totalvisit-1~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailsor,family=poisson)
  tpc.srs[i,]=c(coef(retail.tpc)) ##,coef(summary(retail.tpc))[,2])
  llk.tpcsrs[i]=-logLik(retail.tpc)

   par0<-c(-21.313, 0.049,-0.609,0.152,-0.002,0.133,-0.095,21.968)
  retailrs.tnb<-fun.fitnb(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailsor, par0=par0, dist="Tnb") 
  tnb.srs[i,]=c(retailrs.tnb$par)   ##, sqrt(diag(estvcov))) 
  llk.tnbsrs[i]=retailrs.tnb$value

   par0<-c(-21.313, 0.049,-0.609,0.152,-0.002,0.133,-0.095,21.968)
  retailrs.tnbc<-fun.fitnb(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailsor, par0=par0, dist="OSnb") 
  tnbc.srs[i,]=c(retailrs.tnbc$par) ##, sqrt(diag(estvcov))) 
  llk.tnbcsrs[i]=retailrs.tnbc$value

  retailsor.sor=fitsor(y=retailsor$totalvisitc, x=retailsor[,c("incomec", "distc", "distc2","agec", "MARRIEDc","kidsc")])
  nlg=length(retailsor.sor$lambdagamma)
   uy= sort(unique(retailsor$totalvisit))
  p=c(as.matrix(retailsor.sor$lambdagamma)[1,1:(nlg-6)],0)
  p=exp(p)/sum(exp(p))
   int=sum(uy*p) 
  sor.srs[i,1:7]= c(log(int), as.matrix(retailsor.sor$lambdagamma)[1,(nlg-5):nlg]) 
   llk.sorsrs[i]=retailsor.sor$negllk
}

apply(tpois.srs,2,mean)
apply(tpc.srs,2,mean)
apply(tnb.srs,2,mean)
apply(tnbc.srs,2,mean)
apply(sor.srs,2,mean)


apply(tpois.srs,2,sd)
apply(tpc.srs,2,sd)
apply(tnb.srs,2,sd)
apply(tnbc.srs,2,sd)
apply(sor.srs,2,sd)
sd(sor.srs[,8]-sor.srs[,1])

det(cov(sor.srs[,2:7]))^(1/6)


### onsite sampling
gamma <- c(0.012898, -0.317767,  0.073386, -0.002332,  0.070468, -0.023261) 
lambda <-  c(22.84205, 21.32392, 20.38286, 19.67065, 19.37290, 18.68119, 18.33266, 17.99773, 17.54111, 17.33210, 16.58061, 16.92434,
  16.74171,  15.28760, 15.62247, 16.24222, 15.76093, 14.98963, 14.13334, 13.78327, 14.29084, 13.40592, 13.10757, 13.00548, 12.69196,
   12.36838)     
uuy=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,24,26,28,31,32,35,38, 131)
meany=  1.983182; sy=3.412023
uyc <- (uuy-meany)/sy  ## The uniquely observed values of standardized store visits.
pmat=pred.sor(uy=uyc, Xgrid=retail[,c("incomec", "distc","distc2", "agec","MARRIEDc","kidsc")], lambda=lambda, gamma=gamma)
pmat=t(apply(pmat, 1, function(i) i*uuy/sum(i*uuy)))

for (i in 1:nsim){
   print(i)
   ## Onsite sampling 
   retailsor=retail
  retailsor$totalvisit=apply(pmat, 1,  function(i) uuy[as.vector(rmultinom(1,size=1, i))==1])
  retailsor$totalvisitc= (retailsor$totalvisit-meany)/sy

  ##retailsor$totalvisit=round(retailsor$totalvisitc*sd(retail$totalvisit)+mean(retail$totalvisit),0)
  retailrs.Tpois<-vglm(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailsor,family=pospoisson())
  tpois.ods[i,]=c(coef(retailrs.Tpois)) ##,coef(summary(retailrs.Tpois))[,2])
  llk.tpoisods[i]=-logLik(retailrs.Tpois)

  retail.tpc<-glm(totalvisit-1~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailsor,family=poisson)
  tpc.ods[i,]=c(coef(retail.tpc)) ##,coef(summary(retail.tpc))[,2])
  llk.tpcods[i]=-logLik(retail.tpc)

  par0<-c(-21.313, 0.049,-0.609,0.152,-0.002,0.133,-0.095,21.968)
  retailrs.tnb<-fun.fitnb(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailsor, par0=par0, dist="Tnb") 
  tnb.ods[i,]=c(retailrs.tnb$par)  ##, sqrt(diag(estvcov))) 
  llk.tnbods[i]=retailrs.tnb$value

   par0<-c(-21.313, 0.049,-0.609,0.152,-0.002,0.133,-0.095,21.968)
  retailrs.tnbc<-fun.fitnb(totalvisit~incomec+distc+distc2+agec+MARRIEDc+kidsc, data=retailsor, par0=par0, dist="OSnb") 
  tnbc.ods[i,]=c(retailrs.tnbc$par) ##, sqrt(diag(estvcov))) 
  llk.tnbcods[i]=retailrs.tnbc$value

  retailsor.sor=fitsor(y=retailsor$totalvisitc, x=retailsor[,c("incomec", "distc", "distc2","agec", "MARRIEDc","kidsc")])
  nlg=length(retailsor.sor$lambdagamma)
  uy= sort(unique(retailsor$totalvisit))
  p=c(as.matrix(retailsor.sor$lambdagamma)[1,1:(nlg-6)],0)
  p=exp(p)/sum(exp(p))
  p=p/uy
  p=p/sum(p)
   int=sum(uy*p) 
  sor.ods[i,1:7]= c(log(int), as.matrix(retailsor.sor$lambdagamma)[1,(nlg-5):nlg]) 
   llk.sorods[i]=retailsor.sor$negllk 
}


apply(tpois.ods,2,mean)
apply(tpc.ods,2,mean)
apply(tnb.ods,2,mean)
apply(tnbc.ods,2,mean)
apply(sor.ods,2,mean)


apply(tpois.ods,2,sd)
apply(tpc.ods,2,sd)
apply(tnb.ods,2,sd)
apply(tnbc.ods,2,sd)
apply(sor.ods,2,sd)



det(cov(sor.ods[,2:7]))^(1/6)


#######################################################################################
############ 3.3 Managerial Implications
#####  Expected Shopping Frequencies of Shoppers and Proportions of Returning 
#####    Shoppers by Shopping Distance.
#######################################################################################

dist95=1.974289; dist50=-0.3650186 ; dist05=-1.008669 

## True value
## Using the gamma (log-odds ratio ) estimates for SOR reported in Table 4 and 
## the lambda estimates for baseline function and unique values of visits  reported in Online Appendix I.
gamma <- c(0.012898, -0.317767,  0.073386, -0.002332,  0.070468, -0.023261) 
lambda <-  c(22.84205, 21.32392, 20.38286, 19.67065, 19.37290, 18.68119, 18.33266, 17.99773, 17.54111, 17.33210, 16.58061, 16.92434,
  16.74171,  15.28760, 15.62247, 16.24222, 15.76093, 14.98963, 14.13334, 13.78327, 14.29084, 13.40592, 13.10757, 13.00548, 12.69196,
   12.36838)       
uy <-  c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,24,26,28,31,32,35,38, 131) ## The uniquely observed values of store visits
uyc <- (uy-1.983182)/3.412023  ## The uniquely observed values of standardized store visits. 
sor.freq=pred.sor(uy=uyc, Xgrid=retail[,c("incomec", "distc","distc2", "agec","MARRIEDc","kidsc")], lambda=lambda, gamma=gamma)
sor.fitp=pred.sor(uy=uyc, Xgrid=rbind(c(0,dist95,dist95^2,0,0,0), c(0,dist50,dist50^2,0,0,0), c(0,dist05,dist05^2,0,0,0)), 
  lambda=lambda, gamma=gamma)
## The population true value of Expected shopping frequences by shopping Distance
sor.fitp%*%uy
## The population true value of the proportions of returning shoppers (P(Y>1)) by shopping distance.
1-sor.fitp[,1]


### TP-random 
par.tpc=c(0.224,0.043,-0.515, 0.138,-0.006,0.107,-0.057) 

mean_ztpois(lambda=exp(sum(c(1,0,dist95,dist95^2,0,0,0)*par.tpc)))
mean_ztpois(lambda=exp(sum(c(1,0,dist50,dist50^2,0,0,0)*par.tpc)))
mean_ztpois(lambda=exp(sum(c(1,0,dist05,dist05^2,0,0,0)*par.tpc)))

1-dztpois(1, lambda=exp(sum(c(1,0,dist95,dist95^2,0,0,0)*par.tpc)))
1-dztpois(1, lambda=exp(sum(c(1,0,dist50,dist50^2,0,0,0)*par.tpc)))
1-dztpois(1, lambda=exp(sum(c(1,0,dist05,dist05^2,0,0,0)*par.tpc)))

## TNB Random
par.tnbc=c(-21.313, 0.049,-0.609,0.152,-0.002,0.133,-0.095)
lalpha=21.968

mean_ztnbinom(mu=exp(sum(c(1,0,dist95,dist95^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))
mean_ztnbinom(mu=exp(sum(c(1,0,dist50,dist50^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))
mean_ztnbinom(mu=exp(sum(c(1,0,dist05,dist05^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))

1-dztnbinom(1, mu=exp(sum(c(1,0,dist95,dist95^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))
1-dztnbinom(1, mu=exp(sum(c(1,0,dist50,dist50^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))
1-dztnbinom(1, mu=exp(sum(c(1,0,dist05,dist05^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))



## TP-onsite
par.tpc=c(1.129,0.102,-1.091, 0.305,-0.012,0.217,-0.137)

mean_ztpois(lambda=exp(sum(c(1,0,dist95,dist95^2,0,0,0)*par.tpc)))
mean_ztpois(lambda=exp(sum(c(1,0,dist50,dist50^2,0,0,0)*par.tpc)))
mean_ztpois(lambda=exp(sum(c(1,0,dist05,dist05^2,0,0,0)*par.tpc)))

1-dztpois(1, lambda=exp(sum(c(1,0,dist95,dist95^2,0,0,0)*par.tpc)))
1-dztpois(1, lambda=exp(sum(c(1,0,dist50,dist50^2,0,0,0)*par.tpc)))
1-dztpois(1, lambda=exp(sum(c(1,0,dist05,dist05^2,0,0,0)*par.tpc)))

## TPA-onsite
par.tpc=c(0.823,0.115,-1.270, 0.353,-0.015,0.253,-0.157)

mean_ztpois(lambda=exp(sum(c(1,0,dist95,dist95^2,0,0,0)*par.tpc)))
mean_ztpois(lambda=exp(sum(c(1,0,dist50,dist50^2,0,0,0)*par.tpc)))
mean_ztpois(lambda=exp(sum(c(1,0,dist05,dist05^2,0,0,0)*par.tpc)))

1-dztpois(1, lambda=exp(sum(c(1,0,dist95,dist95^2,0,0,0)*par.tpc)))
1-dztpois(1, lambda=exp(sum(c(1,0,dist50,dist50^2,0,0,0)*par.tpc)))
1-dztpois(1, lambda=exp(sum(c(1,0,dist05,dist05^2,0,0,0)*par.tpc)))

## TNB-onsite
par.tnbc=c(-20.539, 0.098,-1.041,0.293,-0.005,0.198,-0.110)
lalpha=22.738

mean_ztnbinom(mu=exp(sum(c(1,0,dist95,dist95^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))
mean_ztnbinom(mu=exp(sum(c(1,0,dist50,dist50^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))
mean_ztnbinom(mu=exp(sum(c(1,0,dist05,dist05^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))

1-dztnbinom(1, mu=exp(sum(c(1,0,dist95,dist95^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))
1-dztnbinom(1, mu=exp(sum(c(1,0,dist50,dist50^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))
1-dztnbinom(1, mu=exp(sum(c(1,0,dist05,dist05^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))

## TNB-A-onsite
par.tnbc=c(-21.122, 0.089,-0.920,0.260,-0.005,0.184,-0.104)
lalpha=22.195

mean_ztnbinom(mu=exp(sum(c(1,0,dist95,dist95^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))
mean_ztnbinom(mu=exp(sum(c(1,0,dist50,dist50^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))
mean_ztnbinom(mu=exp(sum(c(1,0,dist05,dist05^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))

1-dztnbinom(1, mu=exp(sum(c(1,0,dist95,dist95^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))
1-dztnbinom(1, mu=exp(sum(c(1,0,dist50,dist50^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))
1-dztnbinom(1, mu=exp(sum(c(1,0,dist05,dist05^2,0,0,0)*par.tnbc)), theta=1/exp(lalpha))




