\name{HarshExt}
\alias{HarshExt}
\title{a blemish detection program for microarray chips: extended defects only}
\description{Harshlight automatically detects and masks blemishes in microarray chips of class \code{AffyBatch}}
\usage{
HarshExt(affy.object, my.ErrorImage = NULL, extended.radius = 10)
}
\arguments{
 \item{affy.object}{An AffyBatch object containing two or more chips.}
 \item{my.ErrorImage}{A batch of ErrorImages obtained through other programs. The error images must be in a matrix format, in which the first index represents each cell in the matrix and the second index represents the chip number. By default, the program calculates the error images for the batch of chips affy.object as described in Suarez-Farinas M et al., BMC Bioinformatics - 2005. If a batch of error images is provided, the affy.object is also required.}
 \item{extended.radius}{Radius of the median kernel used to identify extended defects on the chip.}
}
\value{
HarshExt is used to detect only extended defects on the surface of the chip. It does not detect compact or diffuse defects (see the help page for Harshlight).
}
\author{
 Mayte Suarez-Farinas, Maurizio Pellegrino, Knut M. Wittkwosky, Marcelo O. Magnasco
 \email{mpellegri@rockefeller.edu}
}
\references{
	\url{http://asterion.rockefeller.edu/Harshlight/}

	Harshlight: a "corrective make-up" program for microarray chips, Mayte Suarez-Farinas, Maurizio Pellegrino, Knut M Wittkowski and Marcelo O Magnasco, BMC Bioinformatics 2005 Dec 10; 6(1):294

	"Harshlighting" small blemishes on microarrays, Suarez-Farinas M, Haider A, Wittkowski KM., BMC Bioinformatics. 2005 Mar 22;6(1):65.
}
\keyword{file}
