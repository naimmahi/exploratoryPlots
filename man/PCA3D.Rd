% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca3D.R
\name{PCA3D}
\alias{PCA3D}
\title{PCA3D function}
\usage{
PCA3D(eset, logTransform = FALSE, property)
}
\arguments{
\item{eset}{}

\item{logTransform}{If the data should be log transformed? Default is FALSE.}

\item{property}{Property or the variable of interest}
}
\description{
This function allows you to draw 3D PCA plot.
}
\examples{
PCA3D(eset,logTransform=FALSE, property="ER")
}
