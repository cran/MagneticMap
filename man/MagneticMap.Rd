\name{MagneticMap}
\alias{MagneticMap}
\title{
Magnetic Laplacian Matrix and Magnetic Eigenmap Visualization
}
\description{
Constructs the normalized magnetic Laplacian Matrix of a square matrix, returns the eigenvectors and visualization of magnetic eigenmap.
}
\usage{
MagneticMap(W, g, Plot2D=FALSE, Plot3D=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{W}{
  	a numeric or complex matrix whose magnetic Laplacian matrix is to be computed
}
  \item{g}{
		the electric charge parameter that quantifies the importance of the directional     information, 0 <= g <1/2
}
  \item{Plot2D}{
		a logical value indicating whether a 2-D map should be plotted
}
  \item{Plot3D}{
		a logical value indicating whether a 3-D map should be plotted
}
}
\value{
\item{matrix}{
		the normalized magnetic Laplacian matrix
}
\item{eigenvector}{
		the eigenvectors of normalized magnetic Laplacian matrix
}
}
\references{

Fanuel et al. (2016)
\emph{Magnetic Eigenmaps for Visualization of Directed Networks}, preprint arXiv:1606.08266.


}
\examples{
#Create a square matrix
M <- matrix(c(1:16), 4, 4)
#Construct normalized magnetic Laplacian matrix and plot map
m <- MagneticMap(M, g=1/3, Plot2D=TRUE, Plot3D=TRUE)
#See the eigenvectors
vector<-m$eigenvector
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.