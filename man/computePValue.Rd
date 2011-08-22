\name{computePValue}
\alias{computePValue}

\title{ Computes the p-value for a Polya-Aeppli or geometric Poisson
  distribution }

\description{
  Given the two parameters of the Polya-Aeppli distribution and the
  observed motif counts, returns the p-value. These two parameters are
  obtained from the mean \emph{E} and the variance \emph{V} of the counts :
  \itemize{
    \item {%
      \emph{a = ( V[ N(m) ] - E[ N(m) ] ) / ( V[ N(m) ] + E[ N(m) ] )}
    }
    \item{%
       \emph{lambda = (1 - a) / E[ N(m) ]} 
     }
   }
   where \emph{N(m)} is the count distribution of the considered motif
   \emph{m}.

   \code{NA} is returned if Polya-Aeppli parameters \code{lambda, a}
   have bad values (they have to be \emph{lambda > 0} and \emph{0 < a < 1} ).
    
 }

\usage{ computePValue( lambda, a, counts ) }

\arguments{
  \item{lambda}{ Polya-Aeppli parameter, \emph{lambda > 0} }
  \item{a}{ Polya-Aeppli parameter, \emph{0 < a < 1} }
  \item{counts}{ counts of motif occurrences in the network}
}

\value{ Returns the p-value. }

\author{G. Nuel, S. Schbath}

\references{
  [JSCS2006] Nuel G., \emph{Cumulative distribution function of a
    geometric Poisson distribution}. J. Stat. Comp. Sim. 2006 epub. DOI:1080/629360600997371.
} 

\keyword{stats}