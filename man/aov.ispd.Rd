\name{aov.ispd}
\alias{aov.ispd}
\title{Analysis of variance of data from an incomplete split-plot design}
\description{This function performs analysis of variance of data from experiments using
 an incomplete split-plot design for three types of situations namely (i) blocks are complete with
 respect to main plot treatments and mainplots are incomplete with respect to subplot treatments, 
 (ii)  blocks are incomplete with respect to main plot treatments and mainplots are complete with 
 respect to subplot treatments and (iii) blocks are incomplete with respect to main plot treatments and 
mainplots are also incomplete with respect to subplot treatments.}
\usage{aov.ispd(obs, block, mp, sp, y, incomplete.block = FALSE, incomplete.mp = TRUE)}
\arguments{
	\item{obs}{observation numbers} 
	\item{block}{block}
	\item{mp}{main plot treatment}
	\item{sp}{subplot treatment}
	\item{y}{response variable}
	\item{incomplete.block}{Are blocks incomplete? Default is FALSE}
	\item{incomplete.mp}{Are main plots incomplete? Default is TRUE}
} 
\value{Returns ANOVA table of incomplete split-plot design}
\examples{
data(cmis)
with(cmis, aov.ispd(obs, block, mp, sp, y, incomplete.block = FALSE, incomplete.mp = TRUE))
data(imcs)
with(imcs, aov.ispd(obs, block, mp, sp, y, incomplete.block = TRUE, incomplete.mp = FALSE))
data(imis)
with(imis, aov.ispd(obs, block, mp, sp, y, incomplete.block = TRUE, incomplete.mp = TRUE))
}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{incomplete split-plot design}
\keyword{analysis of variance}
\keyword{main plot}
\keyword{subplot}
\keyword{whole plot}
