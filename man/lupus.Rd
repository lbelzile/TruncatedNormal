\name{lupus}
\alias{lupus}
\docType{data}
\title{
Latent Membranous Lupus Nephritis Dataset
}
\description{
The data represents
two clinical measurements (covariates), which are used to predict the occurrence of
latent membranous lupus nephritis. The dataset consists of measurements
on 55 patients of which 18 have been diagnosed with latent membranous lupus.
}
\usage{data("lupus")}
\format{
  The format is:
 "response" "const" "x1" "x2"
}
\details{
 The data was transcribed manually into this format 
 from Table 1, page 22, of Dyk and  Meng (2001).
}

\seealso{
 See Also as \code{\link{mvrandn}}, which uses this dataset.
}


\references{
D. A. van Dyk and X.-L. Meng (2001) \emph{The art of data augmentation (with discussion)}. 
Journal of Computational and  Graphical Statistics, volume 10, pages 1-50.
}
\examples{
 data("lupus")
}
\keyword{datasets}
\keyword{lupus}
