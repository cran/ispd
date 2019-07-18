\name{ispd.imis}
\alias{ispd.imis}
\title{Incomplete split-plot design for given number of main plot treatments, number of subplot treatments, number of blocks, number of main plot treatments in each block and number of subplot treatments in each main plot}
\description{This function generate an incomplete split-plot design for given number of main plot treatments, number of subplot treatments, number of blocks, number of main plot treatments in each block and number of subplot treatments in each main plot}
\usage{ispd.imis(v1,b,k1,v2,k2)}
\arguments{
	\item{v1}{number of main plot treatments}
	\item{b}{number of blocks}	
	\item{k1}{number of main plot treatments in each block}
	\item{v2}{number of subplot treatments}	
	\item{k2}{number of subplot treatments in each main plot}
	}
\value{A list containing parameters, design layout and column layout of design}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{internal}