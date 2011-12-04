\name{partukeycor}
\alias{partukeycor}
\title{partukeycor}
\description{fast tukey biweight correlation matrix}
\usage{
partukeycor(x,
r = 0.2,
output = "matrix",
median = TRUE,
full.init = FALSE,
og = FALSE,
cores = 1)
}
\arguments{
  \item{x}{The matrix of data}
  \item{r}{The breakdown}
  \item{output}{The output format}
  \item{median}{use of median (T)? or covMcd (F)?}
  \item{full.init}{init using all data (T)? or pairwise (F)?}
  \item{og}{use the biwt package estimator?}
  \item{cores}{the number of cores to use}
}

\value{Returns the tukey biweight correlation matrix.}

\author{David L Gibbs}





