
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
                    PACKAGE="fastukey")
       )
  if(is.null(result$output)) {result$output<-NA}
  ### Returns the fitness function result.
  return(result$output)
}


