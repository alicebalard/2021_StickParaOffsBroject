### Custom_PQLseq_Function.R ###

# Created by: Charley
# Date: 20th Mar 2024

# Removed usage of detectCores() to set no. of cores
# This caused job to take use too many cores and set off ITS node alarms, as it ignores how many cores you set as nslots?
# https://www.jottr.org/2022/12/05/avoid-detectcores/ 

# From Sam at ITS, 30/5/24:
# The default numCore value is 1 so hard coding to 6 wont fix this problem unfortunately. Within the function, you'll find:
# cl <- makeCluster(numCore)
# This line spawns threads outside of the R job. By specifying the following instead:
# cl <- makeCluster(numCore, type = "FORK", outfile = "par_log.txt")
# it should refrain from alarming the nodes. You will also have the opportunity to check on any parallel errors with the specified log file.
# As a side note, feel free to set numCores to the number of requested cores (nslots), this should not be affected by the use of detectCores()

########################################################################################################################################

custom_pqlseq <- function(RawCountDataSet, Phenotypes, Covariates=NULL, RelatednessMatrix=NULL, LibSize=NULL, 
                   fit.model="PMM", fit.method = "AI.REML", fit.maxiter=500, fit.tol=1e-5, numCore=nslots, 
                   filtering=TRUE, verbose=FALSE, ...) {
  # specify the number of cores we want to use
  # if(numCore > 1){
  #   if(numCore>detectCores()){warning("PQLseq:: the number of cores you're setting is larger than detected cores!");numCore = detectCores()-1}
  # }
  # 
  # registerDoParallel(numCore)
  
  cl <- makeCluster(numCore,  type = "FORK", outfile = "par_log.txt")
  registerDoParallel(cl,cores=numCore)
  # on.exit(stopCluster(cl))
  
  # filtering genes/sites
  if (filtering & fit.model == "PMM"){
    unfilterIdx <- apply(RawCountDataSet, 1, function(x){length(x[x>5])>=2} )
    CountData   <- RawCountDataSet[unfilterIdx,]
  }else{
    CountData   <- RawCountDataSet
  }
  rm(RawCountDataSet)
  
  numVar <- dim(CountData)[1]
  numIDV <- dim(CountData)[2]
  
  # remove the intercept
  if(length(unique(Covariates[,1])) == 1){
    Covariates<- Covariates[,-1]
  }
  
  if(is.null(Covariates)){
    numCov <- 0
  }else{
    numCov     <- dim(Covariates)[2]
    Covariates <- as.matrix(Covariates)
  }
  
  cat(paste("## number of total individuals: ", numIDV,"\n"))
  cat(paste("## number of total genes/sites: ", numVar,"\n"))
  cat(paste("## number of adjusted covariates: ", numCov,"\n"))
  
  
  CountData  <- as.matrix(CountData)
  Phenotypes <- as.matrix(Phenotypes)
  
  
  if(is.null(RelatednessMatrix)){
    stop("PQLseq::please input relatedness matrix!")
  }else{
    RelatednessMatrix <- as.matrix(RelatednessMatrix)
    scalerM           <- diag(numIDV)-(rep(1,numIDV)%*%t(rep(1,numIDV)))/numIDV
    eig               <- eigen(RelatednessMatrix)
    eigval            <- eig$value
    eigvector         <- eig$vectors
    if(any(eigval<1e-10)){ 
      warning("PQLseq::the relatedness matrix is singular, it has been modified!")
      RelatednessMatrix <- as.matrix(nearPD(RelatednessMatrix,corr=T)$mat)	
    }
    rm(scalerM)
    rm(eig)
    rm(eigval)
    rm(eigvector)
  }
  
  RelatednessMatrix <- list(RelatednessMatrix, diag(numIDV))
  
  #***********************************#
  #       Poisson Mixed Model         #
  #***********************************#
  if(fit.model == "PMM"){
    cat("# fitting Poisson mixed model ... \n")
    if(is.null(LibSize)){
      LibSize <- apply(CountData, 2, sum)
      LibSize <- as.matrix(LibSize)
    }else{
      LibSize <- as.matrix(t(LibSize))
    }
    
    
    # do parallel using foreach function
    iVar   <- NULL
    resPMM <-foreach(iVar=1:numVar,.combine=rbind)%dopar%{
      numAnalysis <- beta <- tau1 <- tau2 <- se_beta <- pvalue <- converged <- h2 <- sigma2 <- overdisp <- NA
      if(numCov==0){
        model0 <- try(glm(formula = CountData[iVar,]~Phenotypes + offset(log(LibSize)), family = poisson(link="log")))
        idx   <- match(rownames(model.frame(formula = CountData[iVar,]~Phenotypes + offset(log(LibSize)), na.action = na.omit)),
                       rownames(model.frame(formula = CountData[iVar,]~Phenotypes + offset(log(LibSize)), na.action = na.pass)))
      }else{
        model0 <- try(glm(formula = CountData[iVar,]~Covariates + Phenotypes + offset(log(LibSize)), family = poisson(link="log")))
        idx   <- match(rownames(model.frame(formula = CountData[iVar,]~Covariates + Phenotypes + offset(log(LibSize)), na.action = na.omit)),
                       rownames(model.frame(formula = CountData[iVar,]~Covariates + Phenotypes + offset(log(LibSize)), na.action = na.pass)))
      }
      
      if(verbose) {cat(paste("NO. Gene = ",iVar,"\n"))}
      
      tmpRelatednessMatrix <- RelatednessMatrix
      if(is(tmpRelatednessMatrix,"matrix")) {
        tmpRelatednessMatrix <- tmpRelatednessMatrix[idx, idx]
      }else {
        for(ik in seq_len(length(tmpRelatednessMatrix)) ) {tmpRelatednessMatrix[[ik]] <- tmpRelatednessMatrix[[ik]][idx, idx]}
      }
      
      names(tmpRelatednessMatrix) <- paste("kins", 1:length(tmpRelatednessMatrix), sep="")
      
      if(class(model0)[1]!="try-error"){
        # t1 <- system.time(model1 <- try(PQLseq.fit(model0, tmpRelatednessMatrix)))
        model1 <- try(PQLseq.fit(model0, tmpRelatednessMatrix))
      }else{
        model1 <- NULL
      }
      
      if(!is.null(model1)&(class(model1)!="try-error")){				
        if(verbose){cat(paste("PQLseq::PMM::tau = ", model1$theta,"\n"))}
        numAnalysis <- length(idx)
        beta        <- model1$coefficients[length(model1$coefficients)]
        se_beta     <- sqrt(diag(model1$cov)[length(model1$coefficients)] )
        pvalue      <- pchisq( (beta/se_beta)^2, 1, lower.tail = F)
        sigma2      <- model1$theta[2]+model1$theta[3]
        h2          <- model1$theta[2]/(sigma2)
        tau1        <- model1$theta[2]
        tau2        <- model1$theta[3]
        converged   <- model1$converged
      }else{converged <- FALSE}
      
      res <- data.frame(numIDV = numAnalysis, beta = beta, se_beta = se_beta, 
                        pvalue = pvalue, h2 = h2, sigma2 = sigma2, 
                        converged = converged) 
    }# end for iVar, parallel
    rm(iVar)
    # closeAllConnections()
    # if(nrow(showConnections())!=0){closeAllConnections()}
    # cons <- suppressWarnings(showConnections(all = TRUE))
    # rm(cons)
    
    parallel::stopCluster(cl)
    
    rownames(resPMM) <- rownames(CountData)
    return(resPMM)
  }# end PMM 
  #***********************************#
  #       Binomial Mixed Model        #
  #***********************************#
  if(fit.model == "BMM"){ 
    cat("# fitting binomial mixed model ... \n")
    if(is.null(LibSize)){
      stop("PQLseq::BMM::ERROR: please input the LibSize (total counts) file!!")
    }else{
      LibSize <- as.matrix(LibSize)
    }
    
    ratio               <- CountData/LibSize
    ratio[is.na(ratio)] <- 0
    flag                <- ratio>1.0
    sumflag             <- apply(flag,1, sum)
    idx                 <- which(sumflag>0)
    
    if (length(idx)>0){
      CountData <- CountData[-idx,]
      LibSize   <- LibSize[-idx,]
    }else{
      CountData <- CountData
      LibSize   <- LibSize
    }
    
    numVar <- dim(CountData)[1]
    numIDV <- dim(CountData)[2]	
    iVar   <- NULL
    
    # do parallel
    resBMM <- foreach(iVar=1:numVar,.combine=rbind)%dopar%{
      numAnalysis <- beta <- tau1 <- tau2 <- se_beta <- pvalue <- converged <- h2 <- sigma2 <- overdisp <- NA
      if(verbose){cat(paste("NO. Gene/Site = ",iVar,"\n"))}
      if(sum(dim(LibSize)==dim(CountData)) != 2){
        stop("PQLseq::BMM::ERROR: the dimensions of read counts and total read counts do not match!")
      }
      
      LibSize <- as.matrix(LibSize)
      
      if(numCov == 0){
        model0 <- glm(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, family = binomial(link = "logit"), weights = LibSize[iVar,])
        idx    <- match(rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, na.action = na.omit)),
                        rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, na.action = na.pass)))
      }else{
        model0 <- glm(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, family = binomial(link = "logit"), weights = LibSize[iVar,] )
        idx    <- match(rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, na.action = na.omit)),
                        rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, na.action = na.pass)))
      }
      
      model0$numTotal <- LibSize[iVar,idx]
      model0$numSucc  <- CountData[iVar,idx]
      
      redflag <- FALSE
      for( ierr in c(2:dim(model.matrix(model0))[2])){
        if(length(unique(model.matrix(model0)[,ierr])) == 1){
          warning(paste("PQLseq::BMM::the ",ierr-1,"-th column of covariates are the same for gene/site ",rownames(CountData)[iVar],"!",sep = "") )
          redflag <- TRUE
        }
      }
      if(!redflag){
        
        tmpRelatednessMatrix <- RelatednessMatrix
        if(is(tmpRelatednessMatrix,"matrix")) {
          tmpRelatednessMatrix <- tmpRelatednessMatrix[idx, idx]
        }else {
          for(ik in seq_len(length(tmpRelatednessMatrix)) ) {
            tmpRelatednessMatrix[[ik]] <- tmpRelatednessMatrix[[ik]][idx, idx]
          }
        }
        names(tmpRelatednessMatrix) <- paste("kins", 1:length(tmpRelatednessMatrix), sep="")
        
        # t1 <- system.time(model1 <- try( PQLseq.fit(model0, tmpRelatednessMatrix) ))
        model1 <- try(PQLseq.fit(model0, tmpRelatednessMatrix))
        
        if(class(model1) != "try-error"&!is.null(model1)){
          if(verbose){cat(paste("PQLseq::BMM::tau = ", model1$theta,"\n"))}
          numAnalysis <- length(idx)
          beta        <- model1$coefficients[ length(model1$coefficients) ]# the last one
          se_beta     <- sqrt( diag(model1$cov)[ length(model1$coefficients) ] )
          pvalue      <- pchisq( (beta/se_beta)^2, 1, lower.tail = F)
          sigma2      <- model1$theta[2]+model1$theta[3]
          h2          <- model1$theta[2]/(sigma2)
          tau1        <- model1$theta[2]
          tau2        <- model1$theta[3]
          converged   <- model1$converged
        }else{converged <- FALSE}
        
        res <- data.frame(numIDV = numAnalysis, beta = beta, se_beta = se_beta, 
                          pvalue = pvalue, h2 = h2, sigma2 = sigma2, 
                          converged = converged)
      }# end for iVar, parallel
      
    }
    rm(iVar)
    
    # if(nrow(showConnections())!=0){closeAllConnections()}
    # closeAllConnections()
    # cons <- suppressWarnings(showConnections(all = TRUE))
    # rm(cons)
    
    parallel::stopCluster(cl)
    
    rownames(resBMM) <- rownames(CountData)
    return(resBMM)
  }# end BMM
  
}# end function PQLseq