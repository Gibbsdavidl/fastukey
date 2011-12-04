# better to pass x
# or construct small matrices in parallel?

parallelBiWt <- function
### parallel tukey estimate call
(x,        ##<< matrix
 isandjs,  ##<< matrix of i and j!
 r,        ##<< the breakdown
 med.init, ##<< the median init value
 full.init,##<< should pairwise init be used?
 median,   ##<< the logic flag for median or covMcd
 originalFlavor = F,  ##<< use the original biwt package estimator?
 id   ##<< the job id
 ) {
  cat("Parallel job ", id, " at work!")
  corr <- vector("numeric", length(isandjs))
  idx <- 1
  for (k in 1:nrow(isandjs)) {
    i <- isandjs[k,1]
    j <- isandjs[k,2]
    if (full.init!=TRUE) { 
      if (median!=TRUE) {
        med.init<-covMcd(cbind(x[i,],x[j,]))
      } else {
        med.init<-list()
        med.init$cov <- diag(1,2)*(apply(cbind(x[i,],x[j,]),2,mad,na.rm=TRUE))^2
        med.init$center <- apply(cbind(x[i,],x[j,]),2,median,na.rm=TRUE)
      }
    }
    if (originalFlavor) {
      biwt <- tryCatch(biwt.est(x=rbind(x[i,],x[j,]),r=r,med.init=med.init),
                       error=function(e){NA})

    } else {
      biwt <- tryCatch(biwt.est.2vecs(cbind(x[i,],x[j,]),r,med.init),
                       error=function(e){NA})
    }
    if (any(is.na(biwt))) {
      corr[idx] <- NA
    } else{
      corr[idx] <- biwt$biwt.sig[1,2]/sqrt(biwt$biwt.sig[1,1]*biwt$biwt.sig[2,2])
    }
    idx <- idx+1
  }
  corr
  ### the tukey correlation vector
}



fullInitFun <- function
(x,         ##<< data matrix
 full.init, ##<< logic indicator for init using full data matrix
 median     ##<< logic indicator for median vs covMcd
 ) {
  if (full.init==TRUE){    
    rand.samp <-x[sample(1:nrow(x),2),]
    if (median != TRUE) {
      med.init <- covMcd(t(rand.samp))
    } else {
      med.init <- list()
      med.init$cov <- diag(1, 2) * (apply(rand.samp, 1, mad, na.rm = TRUE))^2
      med.init$center <- c(1, 1) * apply(rand.samp, 1, median, na.rm = TRUE)
    }
  } else {
    med.init <- NA
  }
  med.init
}



partukeycor <- function
### fast tukey biweight correlation matrix
(x,    ##<< The matrix of data 
 r=.2, ##<< The breakdown
 output="matrix",  ##<< The output format
 median=TRUE,      ##<< use of median (T)? or covMcd (F)?
 full.init=TRUE,   ##<< init using all data (T)? or pairwise (F)?
 og=FALSE, ##<< use the biwt package estimator?
 cores=1           ##<< the number of cores to use
 ){
  require(biwt)
  require(multicore)
  g <- nrow(x)           # number of variables
  n <- (g*(g-1))/2       # upper triangle
  med.init <- fullInitFun(x,full.init, median)
  
  # compute the list of pairs
  print("Computing all pairs ... ")
  x1 <- matrix(seq(2,g))
  x2 <- matrix(1:(g-1))
  p1 <- apply(x1, 2, function(i)rep(i,(i-1)))
  p2 <- unlist(apply(x2, 1, function(i) 1:i))
  isandjs <- matrix(c(p1,p2), ncol=2)

  # break vector into list with "cores" number of pieces
  cat("Breaking up work for ", cores, " number of jobs\n")
  if(cores > 1) {
    b <- cut(1:n, breaks=cores, labels=FALSE)
    ijList <- list()
    for (i in 1:cores) {
      ijList <- c(ijList, list(isandjs[which(b == i),]))
    }
  } else {
    ijList <- list(isandjs)
  }

  # send to parallel jobs.
  print("Starting parallel jobs... ")
  jobs <- lapply(1:cores, function(idx) {
    parallel(parallelBiWt(x,ijList[[idx]],r,med.init, full.init,
                          median, originalFlavor=og, idx),
             name=idx, mc.set.seed = T)
  })
  corr <- unlist(collect(jobs))
  if(output == "matrix") {
    corr.mat <- vect2diss(corr)
    diag(corr.mat) <- 1
  }
  corr.mat
  ### Returns the tukey biweight correlation matrix.
}




