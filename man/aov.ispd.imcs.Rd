\name{aov.ispd.imcs}
\alias{aov.ispd.imcs}
\title{Analysis of variance of data from an incomplete split-plot design with incomplete blocks and complete main plots}
\description{This function performs analysis of variance of data from experiments using
 an incomplete split-plot design for the situation when blocks are incomplete with
 respect to main plot treatments and mainplots are complete with respect to subplot treatments}
\usage{aov.ispd.imcs(obs, block, mp, sp, y)}
\arguments{
	\item{obs}{observation numbers} 
	\item{block}{block}
	\item{mp}{main plot treatment}
	\item{sp}{subplot treatment}
	\item{y}{response variable}
	}
\value{Returns ANOVA table of incomplete split-plot design}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{internal}