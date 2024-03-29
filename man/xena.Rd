% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OOP.R
\name{xena}
\alias{xena}
\title{Create a xena database object.}
\usage{
xena(x, ...)
}
\arguments{
\item{x}{a data frame.}

\item{...}{extra arguments, currently none.}
}
\value{
a xena database which inherits most of its methods from the
data.frame class.
}
\description{
Creates a XENA object on the top of a data frame.
}
\details{
The XENA object enhances functionality of tibbles to hold
mutation data. You may manipulate the XENA object with a plethora of
tidyverse tools in a very similar way as an ordinary tibble.
}
