# This file is part of combination of R-files that contain code to perform a simulation study
# The code was used for a simulation study for which the results are presented in a currently unpublished manuscript

##################################################
#### Functions to estimate performance        ####
#### Last modified: Feb 8 2018                ####
#### Author: M van Smeden                     ####
##################################################

# BEGIN CODE ---------------------------------------

# predict outcome ---
  lp.penal.lrm <- function(glmnet.obj,SIMxy){
    as.vector(predict(object = glmnet.obj,type="link",newx=as.matrix(SIMxy[,-which(colnames(SIMxy)=="y")]),s=glmnet.obj$lambda[1]))
  }  

  lp.firth.lrm <- function(logistf.obj,SIMxy){
    SIMx <- cbind("(Intercept)"=1,SIMxy[,-which(colnames(SIMxy)=="y")])
    unlist(as.matrix(SIMx[,names(coef(logistf.obj))])%*%data.matrix(coef(logistf.obj)))
  }

  predict.lrm <- function(lp){  
    as.vector(exp(lp)/(1+exp(lp)))
  }  

# shrinkage ---
  shrinkage.factor <- function(glmnet.obj){
    glmnet.obj$lambda
  }
  
# prediction statistics ---
 
 c.stat2 <- function(preds, outcome){
   preds <- as.matrix(preds)
   cats <- sort(unique(outcome))
   n_cat <- length(cats)
   n0   <- sum(outcome == cats[2])
   n1   <- length(outcome) - n0
   r <- rank(preds[,1])
   S0 <- sum(as.numeric(r[outcome == cats[2]]))
   (S0 - n0 * (n0 + 1)/2)/(as.numeric(n0) * as.numeric(n1))
 }
 
 
 brier <- function(p,SIMxy){
   sum((SIMxy[,"y"]-p)^2)/length(p)
 }

 MSPE <- function(p,SIMxy,dgm.par){
   design.mat <- as.matrix(cbind(1,SIMxy[,-which(colnames(SIMxy)=="y")]))
   p_true <- exp(design.mat%*%dgm.par)/(1+exp(design.mat%*%dgm.par))
   mean((p_true-p)^2)
 }
 
 MAPE <- function(p,SIMxy,dgm.par){
   design.mat <- as.matrix(cbind(1,SIMxy[,-which(colnames(SIMxy)=="y")]))
   p_true <- exp(design.mat%*%dgm.par)/(1+exp(design.mat%*%dgm.par))
   mean(abs(p_true-p))
 }
  
 bias <- function(p,SIMxy,dgm.par){
   design.mat <- as.matrix(cbind(1,SIMxy[,-which(colnames(SIMxy)=="y")]))
   p_true <- exp(design.mat%*%dgm.par)/(1+exp(design.mat%*%dgm.par))
   c("avr.p.true"=mean(p_true),"avr.p.est"=mean(p),"sd.p.true"=sd(p_true),"sd.p.est"=sd(p), "bias"=mean(p-p_true))
 }
 
 cond.error <- function(p,SIMxy,dgm.par){
   design.mat <- as.matrix(cbind(1,SIMxy[,-which(colnames(SIMxy)=="y")]))
   p_true <- exp(design.mat%*%dgm.par)/(1+exp(design.mat%*%dgm.par))
   cond.err.mat <- matrix(,nrow=7,ncol=2,dimnames=list(c("MSPE","MAPE","avr.p.true","avr.p.est","sd.p.true","sd.p.est","bias"),c("0","1")))
   
   for(class in 0:1){
      cond.err.mat[,class+1] <- c(MSPE(p[SIMxy[,"y"]==class],SIMxy[SIMxy[,"y"]==class,],dgm.par), MAPE(p[SIMxy[,"y"]==class],SIMxy[SIMxy[,"y"]==class,],dgm.par),
        bias(p[SIMxy[,"y"]==class],SIMxy[SIMxy[,"y"]==class,],dgm.par))
   }
   
   list(cond.err.mat=cond.err.mat,discrimination.slope=cond.err.mat["avr.p.est","1"]-cond.err.mat["avr.p.est","0"])
  }
 
 calibration <- function(lp,SIMxy){
    coef(glm(SIMxy[,"y"]~lp, family = "binomial"))
  }


# coefficients ---
  coef.penal.lrm <- function(glmnet.obj){
   t(coef(glmnet.obj))
  }

 coef.firth.lrm <- function(logistf.obj, se.crit=sqrt(5000)){
     list("est" = logistf.obj$coefficients,
          "se" = sqrt(diag(logistf.obj$var)),
          "ci.lo" = logistf.obj$ci.lower,
          "ci.up" = logistf.obj$ci.upper,
          "se.crit.pos" = any(sqrt(diag(logistf.obj$var))>se.crit))  
 }
  
# pseudo-R2 ---
 LL <- function(p, SIMxy){
     sum(SIMxy[,"y"]*log(p)+(1-SIMxy[,"y"])*log(1-p))
 }
   
 pseudo.Rsqrs <- function(p, SIMxy){ 
     LL_fit  <- LL(p=p, SIMxy=SIMxy) 
     LL_null <- LL(p=mean(SIMxy[,"y"]), SIMxy=SIMxy)   
       cox <- 1-exp(-(LL_fit-LL_null)*2/nrow(SIMxy)) 
       cox_max <- 1 - exp(2 * nrow(SIMxy) ^ (-1) * LL_null)
    c("cox"=cox,"nagelkerke"=cox/cox_max,"mcfadden"=1-LL_fit/LL_null)
  }

# obtain stats list ---  
  stats.glmnet.lrm <- function(glmnet.obj,SIMxy,dgm.par){
    lp <- lp.penal.lrm(glmnet.obj,SIMxy)
    p <- predict.lrm(lp)
    
    list("C"=c.stat2(p,SIMxy$y),"brier"=brier(p,SIMxy),"MSPE"=MSPE(p,SIMxy,dgm.par), "MAPE"=MAPE(p,SIMxy,dgm.par),"bias.p"=bias(p,SIMxy,dgm.par),
         "conditional.error"= cond.error(p,SIMxy,dgm.par),"calibration"=calibration(lp, SIMxy),"coef"=coef.penal.lrm(glmnet.obj), 
         "Rsqrs"= pseudo.Rsqrs(p, SIMxy),"shrinkage" = shrinkage.factor(glmnet.obj))
  }

  stats.firth.lrm <- function(logistf.obj,SIMxy,dgm.par,alpha=.10,se.crit=sqrt(5000)){
    lp <- lp.firth.lrm(logistf.obj,SIMxy)
    p <- predict.lrm(lp)
    
     list("C"=c.stat2(p,SIMxy$y),"brier"=brier(p,SIMxy),"MSPE"=MSPE(p,SIMxy,dgm.par), "MAPE"=MAPE(p,SIMxy,dgm.par),"bias.p"=bias(p,SIMxy,dgm.par),
          "conditional.error"= cond.error(p,SIMxy,dgm.par),"calibration"=calibration(lp, SIMxy),"coef"=coef.firth.lrm(logistf.obj), 
          "Rsqrs"= pseudo.Rsqrs(p, SIMxy))
  }



# END CODE ---------------------------------------

  
## CHECK CODE -----------------------------------------
#       source("programs/sim_gendata.R")
#       source("programs/sim_LRMs.R")
#       source("programs/sim_dependencies.R")
#    
#       load("old/Pilot/Pilot2_4var/EPV3/simlist.RData")
#       require(rms)
#        
#        simlist$args$n <- 80000
#        simlist$y1.n <- 40000
#        simlist$y0.n <- simlist$args$n-simlist$y1.n
#        
#        DATA <- random.sampling(simlist,seed=2017)
#        DATA2 <- random.sampling(simlist,seed=2011)
#       
#        SIMxy <- DATA
#        
#        (Firthfit <- firth.lrm.plusbw(DATA,simlist))
#        (GLMnetfit <- penal.cv.lrm2(DATA,nfolds=10))
#        (Hsfit <-  heuristicshrink.lrm(DATA))
#        
#        lp <- lp.firth.lrm(Firthfit$FirthFALSE$value,DATA)
#        p <- predict.lrm(lp)
#          
#        system.time(c.stat2(p,DATA$y))
#      
#       system.time(lasso.stats1 <- stats.glmnet.lrm(glmnet.obj=GLMnetfit$lasso$fit$value,SIMxy=DATA,dgm.par=simlist$dgm.par))
#       system.time(lasso.stats2 <- stats.glmnet.lrm(glmnet.obj=GLMnetfit$lasso$fit$value,SIMxy=DATA2,dgm.par=simlist$dgm.par))
#       system.time(ridge.stats1 <- stats.glmnet.lrm(glmnet.obj=GLMnetfit$ridge$fit$value,SIMxy=DATA,dgm.par=simlist$dgm.par))
#       system.time(ridge.stats2 <- stats.glmnet.lrm(glmnet.obj=GLMnetfit$ridge$fit$value,SIMxy=DATA2,dgm.par=simlist$dgm.par))
#       
#       system.time(Firth.stats1 <- lapply(Firthfit,function(x) stats.firth.lrm(x$value,DATA,dgm.par=simlist$dgm.par)))  
#       system.time(Firth.stats2 <- lapply(Firthfit,function(x) stats.firth.lrm(x$value,DATA2,dgm.par=simlist$dgm.par)))  
#       
#       system.time(HS.stats1 <- lapply(Hsfit$OUT,function(x) stats.firth.lrm(x$value,DATA,dgm.par=simlist$dgm.par)))  
#       system.time(HS.stats2 <- lapply(Hsfit$OUT,function(x) stats.firth.lrm(x$value,DATA2,dgm.par=simlist$dgm.par)))  
#       
## --------------------------------------------------------
# Depreciated functions
#       
#       ECI <- function(lp,p,SIMxy){
#         p.smooth <- predict(vgam(SIMxy$y~s(lp,df=2),binomialff),type="response")
#         mean((p-as.vector(p.smooth))^2)*100
#       }
#       