\name{partukeycor}
\alias{partukeycor}
\title{partukeycor}
\description{fast tukey biweight correlation matrix}
\usage{
partukeycor(x,
r = 0.2,
output = "matrix",
median = FALSE,
full.init = FALSE,
absval = TRUE,
cores = 1)
}
\arguments{
  \item{x}{The matrix of data}
  \item{r}{The breakdown}
  \item{output}{The output format}
  \item{median}{use of median (T)? or covMcd (F)?}
  \item{full.init}{init using all data (T)? or pairwise (F)?}
  \item{absval}{abs value of cors}
  \item{cores}{the number of cores to use}
}

\value{Returns the tukey biweight correlations.}

\author{David L Gibbs}





