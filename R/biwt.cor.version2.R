`fastukeycor` <- function
### fast tukey biweight correlation matrix
(x,    ##<< The matrix of data 
 r=.2, ##<< The breakdown
 output="matrix",  ##<< The output format
 median=TRUE,      ##<< use of median (T)? or covMcd (F)?
 full.init=TRUE,   ##<< init using all data (T)? or pairwise (F)?
 absval=TRUE       ##<< abs value of cors
 ){
  require(biwt)
  rowNums <- nrow(x)
  n <- (rowNums*(rowNums-1))/2 # upper triangle
  corr <- vector("numeric",n)
  g <- dim(x)[1]
  k <- 1

  # init 
  if (full.init==TRUE){    
    rand.samp <-x[sample(1:nrow(x),2),]
    if (median != TRUE) {
      med.init <- covMcd(t(rand.samp))
    } else {
      med.init <- list()
      med.init$cov <- diag(1, 2) * (apply(rand.samp, 1, mad, na.rm = TRUE))^2
      med.init$center <- c(1, 1) * apply(rand.samp, 1, median, na.rm = TRUE)
    }
  }

  # compute the correlations
  for(i in 1:g){
    j <- 1
    while(j < i){
      if (full.init !=TRUE){ 
        if (median!=TRUE) {
          med.init<-covMcd(cbind(x[i,],x[j,]))
        } else {
          med.init<-list()
          med.init$cov <- diag(1,2)*(apply(cbind(x[i,],x[j,]),2,mad,na.rm=TRUE))^2
          med.init$center <- apply(cbind(x[i,],x[j,]),2,median,na.rm=TRUE)
        }
      }
      # cbind here instead of rbind to make the matrix veritcal?
      biwt <- tryCatch(biwt.est.2vecs(cbind(x[i,],x[j,]),r,med.init),
                       error=function(e){NA})
      if (any(is.na(biwt))) {
        corr[k] <- NA
      } else {
        corr[k] <- biwt$biwt.sig[1,2]/sqrt(biwt$biwt.sig[1,1]*biwt$biwt.sig[2,2])
      }
      j<-j+1
      k<-k+1
    }  
  }

  # format the output
  if (output=="matrix") {
    corr.mat <- vect2diss(corr)
    diag(corr.mat) <- 1
  } else if (output=="distance"){   
    if(absval==TRUE){
      corr.mat <- vect2diss(1 - abs(corr))
    } else {
      corr.mat <- vect2diss(1 - corr)
    }    
    diag(corr.mat) <- 0
  } else {
    corr.mat <- corr
  }
  corr.mat
  ### Returns the tukey biweight correlations.
}




