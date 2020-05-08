\name{f_CoVaR_Delta_CoVaR_i_q}
\alias{f_CoVaR_Delta_CoVaR_i_q}

\title{
Computing static CoVaR and Delta CoVaR
}
\description{

This function computes the CoVaR and the Delta CoVaR of a given financial institution i for a given quantile q.
}
\usage{
f_CoVaR_Delta_CoVaR_i_q(df_data_returns)
}

\arguments{
  \item{df_data_returns}{

A dataframe including data: dates and stock returns
}
}
\value{

\item{CoVaR_i_q}{A numeric matrix}
\item{Delta_CoVaR_i_q}{A numeric vector}

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

  # Compute CoVaR_i_q and Delta_CoVaR_i_q
  f_CoVaR_Delta_CoVaR_i_q(data_stock_returns)

}

# }
}


\keyword{Econometrics}
\keyword{Systemic Risk}
