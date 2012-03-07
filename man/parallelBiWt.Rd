\name{parallelBiWt}
\alias{parallelBiWt}
\title{parallelBiWt}
\description{parallel tukey estimate call}
\usage{
parallelBiWt(x,
isandjs,
r,
med.init,
full.init,
median,
originalFlavor = F,
id)
}
\arguments{
  \item{x}{matrix}
  \item{isandjs}{matrix of i and j!}
  \item{r}{the breakdown}
  \item{med.init}{the median init value}
  \item{full.init}{should pairwise init be used?}
  \item{median}{the logic flag for median or covMcd}
  \item{originalFlavor}{use the original biwt package estimator?}
  \item{id}{the job id}
}

\value{the tukey correlation vector}

\author{David L Gibbs}





