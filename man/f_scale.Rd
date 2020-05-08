\name{f_scale}
\alias{f_scale}
\title{
Rescale
}
\description{
This function normalizes data to 0-1 range. Specifically, this function computes linearly rescaled values from a vector of numeric values.
}
\usage{
f_scale(v_time_series)
}
\arguments{
  \item{v_time_series}{
Vector of numeric values
}
}

\value{
A vector of numeric normalized values
}
\references{
Hasse, Jean-Baptiste. "Systemic Risk: a Network Approach". AMSE Working Paper (2020)

}
\author{
Jean-Baptiste Hasse
}

\examples{
# Scale the entries of a vector to the interval [0,1]

# NOT RUN {
\donttest{

  # Generate data
  v_data <- numeric(10)
  v_data <- c(1, 5, 3, 2, 15, 12, 9, 11, 7, 13)

  # Rescale data
  v_rescaled_data <- numeric(10)
  v_rescaled_data <- f_scale(v_data)

  # print rescaled data
  print(v_rescaled_data)
}

# }
}


\keyword{Econometrics}
\keyword{Data}
