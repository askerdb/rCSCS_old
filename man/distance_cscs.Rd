\name{distance_cscs}
\alias{distance_cscs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate CSCS distance
}
\description{
Returns a CSCS distance matrix
}
\usage{
distance_cscs(features, css, cosine_threshold = 0.6, weighted = T, verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{features}{
Matrix of feature intensities in each sample (pxn)
  }
  \item{css}{
  Square matrix of cosine similarities of features (pxp)
  }
  \item{cosine_threshold}{
  Only include features below this threshold
  }
  \item{weighted}{
  Weight features by intensity (TRUE) or absence/presence (FALSE)
  }
  \item{verbose}{
  Print extra details
  }
}
\details{
The value of cosine_threshold is 0.6 in Sedio et al. but for certain applications other values might be better.
}
\value{
A matrix of CSCS distances
}
\references{
Sedio et al, 2016
}
\author{
Asker Brejnrod
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Calculate CSCS for GlobalEuphorbia dataset
distance_cscs(GEfeatures, GEcss)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
