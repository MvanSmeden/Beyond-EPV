# This file is part of combination of R-files that contain code to perform a simulation study
# The code was used for a simulation study for which the results are presented in a currently unpublished manuscript

########################################
#### Functions to estimate models   ####
#### Last modified: Feb 8 2018      ####
#### Author: M van Smeden           ####
########################################

# BEGIN CODE ---------------------------------------
  tryCatch.W.E <- function(expr){
    W <- NULL
    w.handler <- function(w){
      W <<- w
      invokeRestart("muffleWarning")}
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler),warning = W)
  }

  BackwardFR <- function (object, sls,...){
    working <- object
    istep <- 0
    
    scope <- attr(terms(working), "term.labels")
    while ( working$df >= 1) {
      istep <- istep + 1
      mat <- drop1(working)
      if (all(mat[, 3] < sls)) 
        break
      inscope <- match(scope, rownames(mat))
      inscope <- inscope[!is.na(inscope)]
      removal <- rownames(mat)[mat[, 3] == max(mat[inscope, 3])]
      newform = as.formula(paste("~.-", removal))
      if (working$df == 1 | working$df == mat[mat[, 3] == max(mat[,3]), 2]) 
        working <- update(working, formula = newform, pl = FALSE,data=object$data)
      else working <- update(working, formula = newform,data=object$data)      
    }  
    return(working)
  }
  
  penal.cv.lrm2 <- function(SIMxy,nfolds = nrow(SIMxy),type.measure="deviance",nlambda.ridge=200){
    
    OUT <- list("lasso"=list(),"ridge"=list())
    
      ridge.defaultlambda <- glmnet(x=as.matrix(SIMxy[,-which(colnames(SIMxy)=="y")]),y=factor(SIMxy[,"y"]),family="binomial",alpha=0)$lambda
        steps <- log(ridge.defaultlambda[2])-log(ridge.defaultlambda[1])

      OUT$ridge$default.lambda <- c("min"=min(ridge.defaultlambda),"max"=max( ridge.defaultlambda),"length"=length(ridge.defaultlambda))
      OUT$ridge$lambda.sequence <- exp(seq(from=log(max(ridge.defaultlambda)),by=steps, length.out=nlambda.ridge))
      
      OUT$ridge$cv <- tryCatch.W.E(cv.glmnet(x=as.matrix(SIMxy[,-which(colnames(SIMxy)=="y")]),
                            y=factor(SIMxy[,"y"]),family="binomial",
                            nfolds = nfolds, type.measure=type.measure,
                            alpha = 0,lambda=OUT$ridge$lambda.sequence))
      
      OUT$lasso$cv <-  tryCatch.W.E(cv.glmnet(x=as.matrix(SIMxy[,-which(colnames(SIMxy)=="y")]),
                             y=factor(SIMxy[,"y"]),family="binomial",
                             nfolds = nfolds, type.measure=type.measure,
                             alpha = 1))
      
      OUT$lasso$default.lambda <- c("min"=min(OUT$lasso$cv$value$lambda),"max"=max(OUT$lasso$cv$value$lambda),"length"=length(OUT$lasso$cv$value$lambda))
      
      OUT$lasso$fit <-  tryCatch.W.E(glmnet(x=as.matrix(SIMxy[,1:(ncol(SIMxy)-1)]),y=factor(SIMxy$y),
                                            family="binomial",lambda=OUT$lasso$cv$value$lambda.min, alpha = 1))
      
      OUT$ridge$fit <-  tryCatch.W.E(glmnet(x=as.matrix(SIMxy[,1:(ncol(SIMxy)-1)]),y=factor(SIMxy$y),
                                family="binomial",lambda=OUT$ridge$cv$value$lambda.min, alpha = 0))
      
      OUT$lasso$boundary.lambda <- c("default.min"=as.numeric(OUT$lasso$default.lambda["min"])==OUT$lasso$cv$value$lambda.min,
                                     "default.max"=as.numeric(OUT$lasso$default.lambda["max"])==OUT$lasso$cv$value$lambda.min)
      
      OUT$ridge$boundary.lambda <- c("default.min"=as.numeric(OUT$ridge$default.lambda["min"])>=OUT$ridge$cv$value$lambda.min,
                                     "default.max"=as.numeric(OUT$ridge$default.lambda["max"])==OUT$ridge$cv$value$lambda.min,
                                     "new.min" = as.numeric(min(OUT$ridge$lambda.sequence))==OUT$ridge$cv$value$lambda.min)
    
      OUT
  }
  

  firth.lrm.plusbw <- function(SIMxy, simlist){
    formu <- as.formula(paste("y~", paste(colnames(SIMxy)[-which(colnames(SIMxy)=="y")],collapse="+"),sep=""))
    
    TEMPF <- list(FirthTRUE=tryCatch.W.E(logistf(formu,data=SIMxy,firth = TRUE,alpha = simlist$alpha, dataout = T)),
                  FirthFALSE=tryCatch.W.E(logistf(formu,data=SIMxy,firth = FALSE,alpha = simlist$alpha, dataout = T)))
    
    if(!is.null(simlist$bw.p)){
        for(bwiter in 1:length(simlist$bw.p)){    
          TEMPF[[paste("FirthTRUEbw",simlist$bw.p[bwiter],sep="")]] <-  tryCatch.W.E(BackwardFR(TEMPF$FirthTRUE$value,sls=simlist$bw.p[bwiter],  simlist=simlist,firth=TRUE))
          TEMPF[[paste("FirthFALSEbw",simlist$bw.p[bwiter],sep="")]] <- tryCatch.W.E(BackwardFR(TEMPF$FirthFALSE$value,sls=simlist$bw.p[bwiter], simlist=simlist,firth=FALSE)) 
        }
    }   
    TEMPF 
  }

  heuristicshrink.lrm <- function(SIMxy){
    formu <- as.formula(paste("y~", paste(colnames(SIMxy)[-which(colnames(SIMxy)=="y")],collapse="+"),sep=""))
    SIMx <- SIMxy[,-ncol(SIMxy)]
    y <- SIMxy[,ncol(SIMxy)]
    
    TEMP.heur <- list()
    unpen.fit <- glm(formu,data=SIMxy,family=binomial(link = "logit"))
    TEMP.heur$int.re.est <- tryCatch.W.E(logistf(formu,data=SIMxy,firth = FALSE,alpha = simlist$alpha, dataout = T))
      chisq <-unpen.fit$null.deviance-unpen.fit$deviance
        s <- (chisq-TEMP.heur$int.re.est$value$df)/chisq
    
    if(chisq!=0){
    A <- TEMP.heur$int.re.est$value$coefficients["(Intercept)"]
    B <- TEMP.heur$int.re.est$value$coefficients[-which(names(TEMP.heur$int.re.est$value$coefficients)=="(Intercept)")]
    
    TEMP.heur$int.re.est$value$coefficients[-which(names(TEMP.heur$int.re.est$value$coefficients)=="(Intercept)")] <- B*s
    
        TEMP.heur$int.approx <- TEMP.heur$int.re.est    
          TEMP.heur$int.approx$value$coefficients["(Intercept)"] <- (1-s)*mean(SIMxy$y)+s*A
    
        off <- data.matrix(SIMx)%*%data.matrix(B*s)
     
          TEMP.heur$int.re.est$value$coefficients["(Intercept)"] <- coefficients(glm(SIMxy$y~offset(off), family = binomial(link = "logit")))["(Intercept)"]
      }else{TEMP.heur$int.approx <- TEMP.heur$int.re.est}
        
     list(OUT=TEMP.heur,hs.lambda=s)
  }

# END CODE ---------------------------------------------
