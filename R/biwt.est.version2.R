
smallmatmul <- function
### Multiplies a 2x1 and 1x2 matrix.
(v ##<< a vector with length 2
 ) {
  matrix(c(v[1]*v[1],v[1]*v[2],
           v[2]*v[1],v[2]*v[2]),
         ncol=2, byrow=F)
  ### 2x2 matrix
}


`biwt.est.2vecs` <- function
### tukey biweight estimation of just two vectors
(x,   ##<< matrix with two vectors
 r=.2, ##<< the breakdown
 med.init=covMcd(x) ##<< init
 ) {
  require(biwt)
  p<-2
  n <- dim(x)[1] # number of pairs
  c1<-rejpt.bw(p=2,r)[1]
  b0<-erho.bw(p=2,c1)[1]  
  d <- sqrt(mahalanobis(x,med.init$center,med.init$cov))
  k <- fastksolve(d,p,c1,b0)
  if(is.na(k)) {  
    med.init <- covMcd(x)
    d <- sqrt(mahalanobis(x,med.init$center,med.init$cov))
    k <- fastksolve(d,p,c1,b0)
  }                                      # MCD is a more robust estimate of the center/shape

  eps <- 1e-5
  crit <- 100
  iter <- 1
  while (crit > eps & iter < 100) {
    d <- d/k
    biwt.mu <- apply(wtbw(d,c1)*x,2,sum,na.rm=TRUE) / sum (wtbw(d,c1),na.rm=TRUE)
    cent <- array(dim=c(n,p,p))
    for (i in 1:n){
      cent[i,,] <- smallmatmul(x[i,]-biwt.mu)
    }
    biwt.sig <- apply(X=cent*wtbw(d,c1),MARGIN=c(2,3),FUN=sum,na.rm=TRUE)/
      sum(vbw(d,c1),na.rm=TRUE)  
    d2 <- sqrt(mahalanobis(x,biwt.mu,biwt.sig))
    k <- fastksolve(d2,p,c1,b0)
    crit <- max(abs(d-(d2/k)),na.rm=TRUE)
    d <- d2
    iter <-  iter+1
  }
  return(list(biwt.mu=biwt.mu,biwt.sig=biwt.sig))
  ### returns the biweight means and cov matrix.
}

