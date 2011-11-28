
fastksolve <- function
### fast version of ksolve
(d,  ##<< vector of distances
 p,  ##<< p
 c1, ##<< c1
 b0  ##<< b0
 ) {
  try (
       result <- .C("fastksolvec",
                    as.double(d),
                    as.double(p),
                    as.double(c1),
                    as.double(b0),
                    as.integer(length(d)),
                    output=double(1),
                    PACKAGE="fastukeycor")
       )
  if(is.null(result$output)) {result$output<-NA}
  ### Returns the fitness function result.
  return(result$output)
}


testfastksolve <- function
### test fast ksolve with a matrix
(m ##<< matrix with two columns
 ){
  med.init=covMcd(x)
  r <- 0.2
  p<-2
  n <- dim(x)[2]
  c1<-rejpt.bw(p=2,r)[1]
  b0<-erho.bw(p=2,c1)[1]  
  d <- sqrt(mahalanobis(x,med.init$center,med.init$cov))
  k <- fastksolve(d,p,c1,b0)
  k2 <- ksolve(d,p,c1,b0)
  print(k)
  print(k2)
  print(k-k2)
  ### prints test
}
  
