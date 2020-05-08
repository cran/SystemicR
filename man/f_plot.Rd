\name{f_plot}
\alias{f_plot}
\title{
Plot dynamic risk measures
}
\description{
This function provides a framework to plot xts time series.
}
\usage{
f_plot(xts_index_returns)
}
\arguments{
  \item{xts_index_returns}{
A xts vector
}
}

\value{No return value, called for side effects}

\author{
Jean-Baptiste Hasse
}

\examples{
# Plot a xts vector

# NOT RUN {
\donttest{

  # Generate data returns
  v_returns <- numeric(10)
  v_returns <- rnorm(10, 0, 0.01)
  v_date <- seq(from = as.Date("2019-01-01"), to = as.Date("2019-10-01"), by = "month")
  xts_returns <- xts(v_returns, order.by = v_date)

  # Plot the xts vector of simulated returns
  f_plot(xts_returns)
}

# }
}


\keyword{Plot}
\keyword{Illustration}
