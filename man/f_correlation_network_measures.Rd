\name{f_correlation_network_measures}
\alias{f_correlation_network_measures}
\title{
Dynamic systemic risk measures from correlation-based networks.
}
\description{
This function provides methods to compute dynamic systemic risk measures from correlation-based networks.
}
\usage{
f_correlation_network_measures(df_data_returns)
}
\arguments{
  \item{df_data_returns}{
A dataframe including dates and stock returns
}
}

\value{

\item{Degree}{xts vector}
\item{Closeness_Centrality}{xts vector}
\item{Eigenvector_Centrality}{xts vector}
\item{SR}{xts vector}
\item{Volatility}{xts vector}



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

  # Load data
  data("data_stock_returns")

  # Compute topological risk measures from correlation-based financial networks
  l_result <- f_correlation_network_measures(data_stock_returns)

  # Plot SR_t
  f_plot(l_result$SR)

}

# }
}




\keyword{Network}
\keyword{Econometrics}
