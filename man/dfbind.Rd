% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic.R
\name{dfbind}
\alias{dfbind}
\title{Bind a list a data frames}
\usage{
dfbind(x)
}
\arguments{
\item{x}{a list of data frames}
}
\value{
A data frame
}
\description{
Bind a list a data frames
}
\examples{

z <- list(data.frame(x=1, y=1), data.frame(x=2, y=2))
dfbind(z)

}
\seealso{
\link[base]{rbind}
}
