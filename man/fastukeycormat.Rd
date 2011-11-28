\name{fastukeycormat}
\alias{fastukeycormat}
\title{fastukeycormat}
\description{fast tukey biweight correlation matrix}
\usage{
fastukeycormat(x,
r = 0.2,
output = "matrix",
median = TRUE,
full.init = TRUE,
absval = TRUE)
}
\arguments{
  \item{x}{The matrix of data}
  \item{r}{The breakdown}
  \item{output}{The output format}
  \item{median}{use of median (T)? or covMcd (F)?}
  \item{full.init}{init using all data (T)? or pairwise (F)?}
  \item{absval}{abs value of cors}
}

\value{Returns the tukey biweight correlations.}

\author{David L Gibbs}





