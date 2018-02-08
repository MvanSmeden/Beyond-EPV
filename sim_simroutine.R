# This file is part of combination of R-files that contain code to perform a simulation study
# The code was used for a simulation study for which the results are presented in a currently unpublished manuscript

##################################################
#### Execute single sim run and store results ####
#### Last modified: Feb 8 2018                ####
#### Author: M van Smeden                     ####
##################################################

# BEGIN CODE ---------------------------------------

single.run <- function(simlist, simlist.validation, seeds,condno,condID,iter,...){ 
  OUT <- list()

  SIMxy <- do.call(what=simlist$sampling,args=list(simlist=simlist,seed=seeds[1]))
  SIMxy.validate <- do.call(what=simlist.validation$sampling,args=list(simlist=simlist.validation,seed=seeds[2]))
   
  OUT$SIMxy <- list(seeds = seeds, develop = list(degenerate = any(apply(SIMxy,2,var)==0), develop.y = table(SIMxy$y), N = sum(table(SIMxy$y)),
                    means = apply(SIMxy,2,mean), sd = apply(SIMxy,2,sd)),validate = list(degenerate = any(apply(SIMxy.validate,2,var)==0), 
                    validate.y = table(SIMxy.validate$y),N = sum(table(SIMxy.validate$y)),  means = apply(SIMxy.validate,2,mean), sd= apply(SIMxy.validate,2,sd)),
                    lasso.boundarylambda = NA, ridge.boundarylambda = NA, nfolds=NA)
 
  if(!any(apply(SIMxy,2,var)==0)&min(table(SIMxy$y))>2){
    TEMP.logistf <- firth.lrm.plusbw(SIMxy=SIMxy, simlist=simlist)
    TEMP.hs <- heuristicshrink.lrm(SIMxy)
    set.seed(seeds[1])
    
    nfolds <- ifelse(any(table(SIMxy$y)<8),nrow(SIMxy),10)
    TEMP.glmnet <- penal.cv.lrm2(SIMxy,nfolds=nfolds)
    
    OUT$logistf <- list(develop=lapply(TEMP.logistf,function(x) stats.firth.lrm(logistf.obj=x$value,SIMxy=SIMxy,dgm.par=simlist$dgm.par,alpha=simlist$alpha,se.crit=simlist$se.crit)))   
    OUT$hs <- list(develop=lapply(TEMP.hs$OUT,function(x) stats.firth.lrm(logistf.obj=x$value,SIMxy=SIMxy,dgm.par=simlist$dgm.par,alpha=simlist$alpha,se.crit=simlist$se.crit)),
                   hs.lambda = TEMP.hs$hs.lambda)   
    OUT$glmnet <- list(develop=list(ridge=stats.glmnet.lrm(glmnet.obj=TEMP.glmnet$ridge$fit$value,SIMxy=SIMxy,dgm.par=simlist$dgm.par),
                                    lasso=stats.glmnet.lrm(glmnet.obj=TEMP.glmnet$lasso$fit$value,SIMxy=SIMxy,dgm.par=simlist$dgm.par)),
                       lambdaseq=list(ridge.default=TEMP.glmnet$ridge$default.lambda,lasso.default=TEMP.glmnet$lasso$default.lambda))
    
    OUT$logistf$validate <- lapply(TEMP.logistf,function(x) stats.firth.lrm(logistf.obj=x$value,SIMxy=SIMxy.validate,dgm.par=simlist$dgm.par,alpha=simlist$alpha,se.crit=simlist$se.crit))   
    OUT$hs$validate <- lapply(TEMP.hs$OUT,function(x) stats.firth.lrm(logistf.obj=x$value,SIMxy=SIMxy.validate,dgm.par=simlist$dgm.par,alpha=simlist$alpha,se.crit=simlist$se.crit))
    OUT$glmnet$validate <- list(ridge=stats.glmnet.lrm(glmnet.obj=TEMP.glmnet$ridge$fit$value,SIMxy=SIMxy.validate,dgm.par=simlist$dgm.par),
                                lasso=stats.glmnet.lrm(glmnet.obj=TEMP.glmnet$lasso$fit$value,SIMxy=SIMxy.validate,dgm.par=simlist$dgm.par))
    
	OUT$SIMxy$lasso.boundarylambda = TEMP.glmnet$lasso$boundary.lambda
	OUT$SIMxy$ridge.boundarylambda = TEMP.glmnet$ridge$boundary.lambda
	OUT$SIMxy$nfolds=nfolds
   
  }
 
  dir.create("OUT", showWarnings = FALSE)
  path <- file.path("OUT",condID)
  dir.create(path,showWarnings = FALSE)
  saveRDS(OUT,file=file.path(path,paste(paste(condID,sprintf("%05.0f",iter),sep="_"),".rds",sep="")))   
}


# END CODE ---------------------------------------

  
