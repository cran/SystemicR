\name{f_CoVaR_Delta_CoVaR_i_q_t}
\alias{f_CoVaR_Delta_CoVaR_i_q_t}

\title{
Computing dynamic CoVaR and Delta CoVaR
}
\description{

This function computes the dynamic CoVaR and the Delta CoVaR of a given financial institution i for a given quantile q at time t. The dynamic and aggregate Delta CoVaR is also computed.
}
\usage{
f_CoVaR_Delta_CoVaR_i_q_t(df_data_returns, df_data_state_variables)
}

\arguments{
  \item{df_data_returns}{

A dataframe including data: dates and stock returns
}
  \item{df_data_state_variables}{

A dataframe including data: dates and macroeconomic variables

}
}

\value{

\item{CoVaR_i_q_t}{A xts matrix}
\item{Delta_CoVaR_i_q_t}{A xts matrix}
\item{Delta_CoVaR_t}{A xts vector}

%% ...
}
\references{

Adrian, Tobias, and Markus K. Brunnermeier. "CoVaR". American Economic Review 106.7 (2016): , 106, 7, 1705-1741.
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
  data("data_state_variables")


  # Compute CoVaR_i_q_t , Delta_CoVaR_i_q_t and Delta_CoVaR_t
  l_result <- f_CoVaR_Delta_CoVaR_i_q_t(data_stock_returns, data_state_variables)

  # Plot Delta_CoVaR_t
  f_plot(l_result$Delta_CoVaR_t)

}

# }
}


\keyword{Econometrics}
\keyword{Systemic Risk}
