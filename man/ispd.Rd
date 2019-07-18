\name{ispd}
\alias{ispd}
\title{Incomplete split-plot design for given number of blocks, number of main plot treatments, number of subplot treatments, number of main plot treatments in blocks and / or number of subplot treatments in main plots}
\description{This function generates an incomplete split-plot design with given number of 
main plot treatments(v1), number of subplot treatments (v2), number of blocks(b) and block size(k). The incomplete split-plot design may be one of the three kinds:
(i) blocks are complete with respect to main plot treatments and mainplots are incomplete with respect to subplot treatments, (ii)  blocks are incomplete with 
respect to main plot treatments and mainplots are complete with respect to subplot treatments and (iii) blocks are incomplete with respect to main plot treatments and 
mainplots are also incomplete with respect to subplot treatments.}
\usage{ispd(v1,v2,b,k1 = NULL,k2 = NULL)}
\arguments{
	\item{v1}{number of main plot treatments}
	\item{v2}{number of subplot treatments}
	\item{b}{number of blocks}
	\item{k1}{number of main plot treatments in each block. If k1 is not specified, it is assumed that k1 = v1}
	\item{k2}{number of subplot treatments in each main plot. If k2 is not specified, it is assumed that k2 = v2}
	}
\value{A list containing parameters, design layout and column layout of design}
\examples{
ispd(v1 = 3, v2 = 4, b = 3, k1 = 2)
ispd(v1 = 3, v2 = 3, b = 3, k2 = 2)
ispd(v1 = 4, b = 6, k1 = 2, v2 = 3, k2 = 2)
}

\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{incomplete split-plot design}
\keyword{main plot}
\keyword{subplot}
\keyword{whole plot}
\keyword{analysis}