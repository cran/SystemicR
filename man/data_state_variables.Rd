\name{data_state_variables}
\alias{data_state_variables}
\docType{data}
\title{
State variables
}
\description{
This dataset includes state variables data extracted from the FRED. Specifically, it includes data on credit spread, liquidity spread, yield spread, 3M Treasury bill and VIX.
}
\usage{data("data_state_variables")}
\format{
  A data frame with 5030 observations on the following 7 variables.
  \describe{
    \item{\code{Date}}{a date vector}
    \item{\code{CRESPR}}{a numeric vector}
    \item{\code{LIQSPR}}{a numeric vector}
    \item{\code{YIESPR}}{a numeric vector}
    \item{\code{TBR3M}}{a numeric vector}
    \item{\code{RESI}}{a numeric vector}
    \item{\code{VIX}}{a numeric vector}

  }
}

\source{
Federal Reserve Economic Data (FRED) St. Louis Fed
}
\references{
Hasse, Jean-Baptiste. "Systemic Risk: a Network Approach". AMSE Working Paper (2020)
Hasse, Jean-Baptiste, and Quentin Lajaunie. "Does the Yield Curve Signal Recessions? New Evidence from an International Panel Data Analysis." AMSE Working Paper (2020).
}
\examples{
data("data_state_variables")
head(data_state_variables)
}
\keyword{datasets}
