# This file is part of combination of R-files that contain code to perform a simulation study
# The code was used for a simulation study for which the results are presented in a currently unpublished manuscript

########################################
#### Functions to generate data     ####
#### Last modified: Feb 8 2018      ####
#### Author: M van Smeden           ####
########################################

# BEGIN CODE ---------------------------------------
  gen.independentbinX <- function(n,margp.x){
    no.X <- length(margp.x)
    SIMx <- matrix(ncol=no.X,nrow=n)
    for(i in 1:no.X) SIMx[,i]<-rbinom(n,1,margp.x[i])
    SIMx
  }
  
  gen.MVNX <- function(n,mu,sigma){
     mvrnorm(n=n,mu=mu,Sigma=sigma)
  }

  gen.uniformX <- function(n,no.X,min=0,max=1){
    SIMx <- matrix(ncol=no.X,nrow=n)
    for(i in 1:no.X) SIMx[,i]<-runif(n,min=min,max=max)
    SIMx
  }
  
  gen.binY <- function(SIMx,dgm.par){
    no.X <- ncol(SIMx)
    design.mat <- cbind(1,SIMx)
    p <- exp(design.mat%*%dgm.par)/(1+exp(design.mat%*%dgm.par))
    y <- rbinom(length(p),1,p)
    SIMxy <-data.frame(x=SIMx,y=y)
    colnames(SIMxy)[1:no.X]<-paste("X",1:no.X,sep="")  
    SIMxy
  }

  random.sampling <- function(simlist,seed){ 
    set.seed(seed)
    SIMx <- do.call(what=simlist$datagen,args=simlist$args)
    SIMxy <- gen.binY(SIMx=SIMx,dgm.par=simlist$dgm.par)
    SIMxy
  }

  retrospect.sampling <- function(simlist,seed){ 
    set.seed(seed)
    SIMx <- do.call(what=simlist$datagen,args=simlist$args)
    SIMxy <- gen.binY(SIMx=SIMx,dgm.par=simlist$dgm.par)
    
    simlist$args$n<-2
    while(sum(SIMxy$y)<simlist$y1.n | sum(1-SIMxy$y)<simlist$y0.n){
      TEMPx <- do.call(what=simlist$datagen,args=simlist$args)
      SIMxy <- rbind(SIMxy,gen.binY(SIMx=TEMPx,dgm.par=simlist$dgm.par))
    }
    
    rbind(SIMxy[order(SIMxy$y),][1:simlist$y0.n,],SIMxy[order(1-SIMxy$y),][1:simlist$y1.n,])
  }

  data.arg <- function(condition){
    if(condition$datagen=="gen.MVNX"){
      sigma <- matrix(as.numeric(unlist(as.character(condition$datacor))),ncol=condition$P,nrow=condition$P)
      diag(sigma)<-1
      args <- list(mu=rep(0,condition$P),sigma=sigma)}
    if(condition$datagen=="gen.independentbinX"){args<-list(margp.x=rep(as.numeric(unlist(as.character(condition$margp.x))),condition$P))}
    if(condition$datagen=="gen.uniform"){
      args <-list(no.X=condition$P,min=-3,max=3)}
    list(datagen=as.character(condition$datagen),args=args)
  }


# END CODE ---------------------------------------
