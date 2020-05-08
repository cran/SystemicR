#####################################################################################
#                                                                                   #
#                                       Functions                                   #
#                                                                                   #
#####################################################################################

#####################################################################################
#                                     Scale norm                                    #
#####################################################################################

f_scale <- function(v_time_series)
{
  v_scaled_time_series <- (v_time_series - min(v_time_series)) / (max(v_time_series) - min(v_time_series))
  return(v_scaled_time_series)
}

#####################################################################################
#                                        Plot                                       #
#####################################################################################

f_plot <- function(xts_index_returns)
{
  plot(xts_index_returns, type = "l", main = "SR", grid.col = NA, yaxis.right = FALSE)
}

#####################################################################################
#                             CoVaR_i_q and Delta_CoVaR_i_q                         #
#####################################################################################

f_CoVaR_Delta_CoVaR_i_q <- function(df_data_returns)
{
  Number_Institutions <- length(df_data_returns[1,]) - 2
  df_data_2 <- df_data_returns[,2:(Number_Institutions + 1)]
  m_data_returns <- as.matrix(df_data_returns[,2:(Number_Institutions + 2)])

  # Eq. 7 - Estimation of quantile regressions

  q_95 <- 0.95
  q_50 <- 0.50

  res_qr <- matrix(0, ncol = Number_Institutions, nrow = 2)

  for (cpt_institutions in 1:Number_Institutions)
  {
    tampon_qr <- quantreg::rq(m_data_returns[,1] ~ m_data_returns[,cpt_institutions + 1], data = df_data_2, tau = (1 - q_95))
    res_qr[1,cpt_institutions] <- tampon_qr$coef[[1]]
    res_qr[2,cpt_institutions] <- tampon_qr$coef[[2]]
  }

  # Eq. 8 - Estimation of VaR_i_q

  m_VaR_95_50 <- matrix(0, nrow = 2, ncol = Number_Institutions)

  for (cpt_institutions in 1:Number_Institutions)
  {
    mean <- mean(m_data_returns[,cpt_institutions + 1])
    sd <- sd(m_data_returns[,cpt_institutions + 1])
    m_VaR_95_50[1,cpt_institutions] <- qnorm((1 - q_95),mean,sd)
    m_VaR_95_50[2,cpt_institutions] <- qnorm((1 - q_50),mean,sd)

  }

  # Eq. 9 - Estimation of CoVaR_i_q

  m_CoVaR_95_50 <- matrix(0, nrow = 2, ncol = Number_Institutions)

  for (cpt_institutions in 1:Number_Institutions)
  {
    m_CoVaR_95_50[1,cpt_institutions] <- res_qr[1,cpt_institutions] + res_qr[2,cpt_institutions] * m_VaR_95_50[1,cpt_institutions]
    m_CoVaR_95_50[2,cpt_institutions] <- res_qr[1,cpt_institutions] + res_qr[2,cpt_institutions] * m_VaR_95_50[2,cpt_institutions]
  }

  # Eq. X - Estimation of Delta_CoVaR_i_q

  v_Delta_CoVaR <- vector("numeric", length = Number_Institutions)

  for (cpt_institutions in 1:Number_Institutions)
  {
    v_Delta_CoVaR[cpt_institutions] <- res_qr[2,cpt_institutions] * (m_VaR_95_50[1,cpt_institutions] - m_VaR_95_50[2,cpt_institutions])
  }

  result <- list("CoVaR_i_q" = m_CoVaR_95_50, "Delta_CoVaR_i_q" = v_Delta_CoVaR)
  return(result)
}

#####################################################################################
#                 CoVaR_i_q_t , Delta_CoVaR_i_q_t and Delta_CoVaR_t                 #
#####################################################################################

f_CoVaR_Delta_CoVaR_i_q_t <- function(df_data_returns, df_data_state_variables)
{
  #df_data_returns <- data_stock_returns
  Number_Institutions <- length(df_data_returns[1,]) - 2
  m_data_returns <- as.matrix(df_data_returns[,2:(Number_Institutions + 2)])


  #df_data_state_variables <- data_state_variables
  Number_State_Variables  <- length(df_data_state_variables[1,]) - 2
  m_data_state_variables <- as.matrix(df_data_state_variables[,2:(Number_State_Variables + 2)])

  Number_Observations <- length(m_data_returns[,1])
  Number_Observations_2 <- length(m_data_state_variables[,1])


  # State variables (contemporaneous and lagged respectively)
  # RESI   : The weekly real estate sector return in excess of the market financial sector return (from the real estate companies with SIC code 65-66)
  # VIX    : Equity volatility, which is computed as the 22-day rolling standard deviation of the daily CRSP equity market retur
  # TBR3M  : Change in the three-month yield: Change in the three-month Treasury bill rate
  # CRESPR : The change in the credit spread between Moody?s Baa-rated bonds and the ten year Treasury rate from the Federal Reserve Board?
  # LIQSPR : A short term "TED spread":  difference between the threemonth LIBOR rate and the three-month secondary market treasury bill rate. This spread measures short-term funding liquidity risk. We use the three-month LIBOR rate that is available from the British Bankers? Association, and obtain the three-month Treasury rate from the Federal Reserve Bank of New Yor
  # YIESPR : Change in the slope of the yield curve:  spread between the composite long-term bond yield and the three-month bill rate
  # The weekly market return computed from the S&P500 (not implemented)

  M_t <- matrix(0, nrow = Number_Observations, ncol = Number_State_Variables)
  M_t_1 <- matrix(0, nrow = (Number_Observations - 1), ncol = Number_State_Variables)

  M_t <- m_data_state_variables[,1:Number_State_Variables]
  M_t_1 <- M_t[1:(Number_Observations - 1),]

  # Banks' returns (with one missing observation to be regressed with lagged state variables)

  X_i_t <- matrix(0, nrow = (Number_Observations - 1), ncol = Number_Institutions)
  X_i_t <- m_data_returns[2:Number_Observations,]

  q_95 <- 0.95
  q_50 <- 0.50

  # Eq. 11a - Estimation of quantile regressions

  # quantile 95

  res_qr_2 <- matrix(0, ncol = Number_Institutions, nrow = Number_State_Variables + 1)

  for (cpt_institutions in 1:Number_Institutions)
  {
    tampon_qr_2 <- quantreg::rq(X_i_t[,cpt_institutions] ~ M_t_1[,1:Number_State_Variables], tau = (1 - q_95))
    res_qr_2[1,cpt_institutions] <- tampon_qr_2$coef[[1]]
    res_qr_2[2:(Number_State_Variables + 1),cpt_institutions] <- tampon_qr_2$coef[2:(Number_State_Variables + 1)]

  }

  # quantile 50

  res_qr_32 <- matrix(0, ncol = Number_Institutions, nrow = 1 + Number_State_Variables + 1)

  for (cpt_institutions in 1:Number_Institutions)
  {
    tampon_qr_32 <- quantreg::rq(X_i_t[,1] ~ M_t_1[,1:Number_State_Variables] + X_i_t[,cpt_institutions + 1], tau = (1 - q_50))
    res_qr_32[1,cpt_institutions] <- tampon_qr_32$coef[[1]]
    res_qr_32[2:(1 + Number_State_Variables + 1),cpt_institutions] <- tampon_qr_32$coef[2:(1 + Number_State_Variables + 1)]

  }

  # Eq. 11b - Estimation of quantile regressions

  res_qr_3 <- matrix(0, ncol = Number_Institutions, nrow = 1 + Number_State_Variables + 1)

  for (cpt_institutions in 1:Number_Institutions)
  {
    tampon_qr_3 <- quantreg::rq(X_i_t[,1] ~ M_t_1[,1:Number_State_Variables] + X_i_t[,cpt_institutions + 1], tau = (1 - q_50))
    res_qr_3[1,cpt_institutions] <- tampon_qr_3$coef[[1]]
    res_qr_3[2:(1 + Number_State_Variables + 1),cpt_institutions] <- tampon_qr_3$coef[2:(1 + Number_State_Variables + 1)]

  }

  # Eq. 12a - Estimation of VaR_q_i_t

  # quantile 95

  m_VaR_iqt <- matrix(0, nrow = (Number_Observations - 1), ncol = Number_Institutions)

  for (cpt_institutions in 1:Number_Institutions)
  {
    alpha <- res_qr_2[1,cpt_institutions]
    beta <- res_qr_2[2:(Number_State_Variables + 1),cpt_institutions]
    m_VaR_iqt[,cpt_institutions] <- t(alpha + beta %*% t(M_t_1))
  }

  # quantile 50

  m_VaR_iqt_50 <- matrix(0, nrow = (Number_Observations - 1), ncol = Number_Institutions)

  for (cpt_institutions in 1:Number_Institutions)
  {
    alpha <- res_qr_32[1,cpt_institutions]
    beta <- res_qr_32[2:(Number_State_Variables + 1),cpt_institutions]
    m_VaR_iqt_50[,cpt_institutions] <- t(alpha + beta %*% t(M_t_1))
  }

  # Eq. 12b - Estimation of CoVaR_q_i_t

  m_CoVaR_iqt <- matrix(0, nrow = (Number_Observations - 1), ncol = Number_Institutions)

  for (cpt_institutions in 1:Number_Institutions)
  {
    alpha <- res_qr_3[1,cpt_institutions]
    beta <- res_qr_3[2:(1 + Number_State_Variables + 1),cpt_institutions]
    m_CoVaR_iqt[,cpt_institutions] <- t(alpha + beta[1:Number_State_Variables] %*% t(M_t_1) + beta[(Number_State_Variables + 1)] %*% t(m_VaR_iqt[,cpt_institutions]))
  }

  # Eq. 13-14 - Estimation of Delta_CoVaR_q_i_t

  m_D_CoVaR_iqt <- matrix(0, nrow = (Number_Observations - 1), ncol = Number_Institutions)

  for (cpt_institutions in 1:Number_Institutions)
  {
    beta <- res_qr_3[2:(1 + Number_State_Variables + 1),cpt_institutions]
    m_D_CoVaR_iqt[,cpt_institutions] <- beta[(Number_State_Variables + 1)] * (m_VaR_iqt[,cpt_institutions] - m_VaR_iqt_50[,cpt_institutions])
  }

  # Delta CoVaR equi-weighted
  v_D_CoVaR_t <- vector("numeric", length = (Number_Observations - 1))

  for (cpt_time in 1:Number_Observations)
  {
    v_D_CoVaR_t[cpt_time] <- 100 * ( - mean(m_D_CoVaR_iqt[cpt_time - 1,]))
  }

  xts_CoVaR_iqt <- xts::xts(m_CoVaR_iqt, order.by=as.Date(df_data_returns[2:Number_Observations,1],"%d/%m/%Y"))
  xts_D_CoVaR_iqt <- xts::xts(m_D_CoVaR_iqt, order.by=as.Date(df_data_returns[2:Number_Observations,1],"%d/%m/%Y"))
  xts_D_CoVaR_t <- xts::xts(v_D_CoVaR_t, order.by=as.Date(df_data_returns[,1],"%d/%m/%Y"))

  result <- list("CoVaR_i_q_t" = xts_CoVaR_iqt, "Delta_CoVaR_i_q_t" = xts_D_CoVaR_iqt, "Delta_CoVaR_t" = xts_D_CoVaR_t)

  return(result)

}

#####################################################################################
#                                 Network Measures                                  #
#####################################################################################

f_correlation_network_measures <- function(df_data_returns)
{
  Number_Institutions <- length(df_data_returns[1,]) - 2
  df_data_2 <- df_data_returns[,3:(Number_Institutions + 1)]
  m_data_returns <- as.matrix(df_data_returns[,2:(Number_Institutions + 2)])
  Number_Observations <- length(m_data_returns[,1])


  #                                   Data treatment                                  #

  # Declare vectors in which we're gonna extract variables we need to expoloit data
  Test_Change_Month <- vector(mode="numeric", length=Number_Observations)
  Extract_Month <- vector(mode="numeric", length=Number_Observations)

  # Allocation and management of dates data
  data_date_brut <-   df_data_returns[,1]
  data_date <- as.Date(data_date_brut, "%d/%m/%Y")

  # - (Build a vector of dates from the initial date of the downloaded data and) -  extract the month of each dates
  for (cpt in 1:(Number_Observations))
  {
    Extract_Month[cpt] <- as.numeric(format(as.Date(data_date[cpt]), "%m"))
  }

  # Build a vector testing (and so identifying) the month changement trough the vector of dates
  for (cpt in 1:(Number_Observations - 1))
  {
    if (Extract_Month[cpt] != Extract_Month[cpt+1]){
      Test_Change_Month[cpt] <- 1
    }else{
      Test_Change_Month[cpt] <- 0
    }
  }

  # Define the number of months of the dataset in order to build the new dataset (one column per month)
  Number_Months <- sum(Test_Change_Month)

  a_data <- array(0, dim = c(31, Number_Months, Number_Institutions))

  # Loop going through institutions
  for (cpt_institutions in 2:(Number_Institutions + 1))
  {

    # Initialization of the external cursor and declaration of the new dataset
    cpt_vect <- 0

    # Extraction of the month of the corresponding date of the target cell (later we'll replace Extract_Month[cpt_vect] because what we are interested in is the return)
    for (cpt_col in 1:Number_Months)
    {
      for (cpt_row in 1:31)
      {
        cpt_vect <- cpt_vect + 1
        if (Test_Change_Month[cpt_vect] != 1){
          a_data[cpt_row,cpt_col,cpt_institutions - 1] <- df_data_returns[cpt_vect,cpt_institutions]
        }else{
          break
        }
      }
    }
  }

  #                           Estimation of Realized Volatility                       #

  # Declaration of matrices
  m_mean <- matrix(0, ncol = Number_Institutions, nrow = Number_Months)
  m_sd <- matrix(0, ncol = Number_Institutions, nrow = Number_Months)
  m_Sum_Squared_Returns <- matrix(0, ncol = Number_Institutions, nrow = Number_Months)
  m_Realized_Volatility <- matrix(0, ncol = Number_Institutions, nrow = Number_Months)

  # Loop going through variables
  for (cpt_institutions in 1:Number_Institutions)
  {
    # Estimation of realized volatility
    for (cpt_row in 1:Number_Months)
    {
      tampon <- a_data[,cpt_row,cpt_institutions]
      length_tampon <- length(tampon[tampon != 0])
      m_mean[cpt_row,cpt_institutions] <- mean(tampon[tampon != 0])
      m_sd[cpt_row,cpt_institutions] <- sd(tampon[tampon != 0])
      tampon_2 <- tampon[tampon != 0]
      m_Sum_Squared_Returns[cpt_row,cpt_institutions] <- sum((tampon_2 - m_mean[cpt_row,cpt_institutions])^2)
      m_Realized_Volatility[cpt_row,cpt_institutions] <- sqrt((1/length_tampon) * m_Sum_Squared_Returns[cpt_row,cpt_institutions])
    }
  }

  # NaN appears we should have a look later about what is going on
  m_Realized_Volatility[is.nan(m_Realized_Volatility)] <- 0

  #                            Estimation of Realized Covariance                      #

  # Declaration of covariance matrix
  a_covar <- array(0, dim = c(Number_Months, Number_Institutions, Number_Institutions))
  m_tampon_symetric <- matrix(0, ncol = Number_Institutions, nrow = Number_Institutions)

  # Loop to fullfill the covariance matrix
  for (cpt_month in 1:Number_Months)
  {
    for (cpt_institutions_1 in 1:Number_Institutions)
    {
      for (cpt_institutions_2 in 1:Number_Institutions)
      {
        tampon_1 <- a_data[,cpt_month,cpt_institutions_1]
        tampon_2 <- a_data[,cpt_month,cpt_institutions_2]
        tampon_11 <- tampon_1[tampon != 0]
        tampon_22 <- tampon_2[tampon != 0]
        mean_1 <- mean(tampon_11)
        mean_2 <- mean(tampon_22)
        T <- length(tampon_1[(tampon_1 != 0)])
        a_covar[cpt_month, cpt_institutions_1, cpt_institutions_2] <- (1/T) * sum((tampon_1 - mean_1) * (tampon_2 - mean_2))

      }
    }

    m_tampon_symetric <- a_covar[cpt_month, , ]
    m_tampon_symetric <- Matrix::forceSymmetric(m_tampon_symetric)
    a_covar[cpt_month, , ] <- as.matrix(m_tampon_symetric)
  }

  #                           Estimation of Realized Correlations                     #

  # Declaration of covariance matrix
  a_cor <- array(0, dim = c(Number_Months, Number_Institutions, Number_Institutions))

  # Loop to fullfill the covariance matrix
  for (cpt_month in 1:Number_Months)
  {
    for (cpt_institutions_1 in 1:Number_Institutions)
    {
      for (cpt_institutions_2 in 1:Number_Institutions)
      {
        a_cor[cpt_month, cpt_institutions_1, cpt_institutions_2] <- (a_covar[cpt_month, cpt_institutions_1, cpt_institutions_2]) / ((m_Realized_Volatility[cpt_month,cpt_institutions_1]) * (m_Realized_Volatility[cpt_month,cpt_institutions_2]))
      }
    }
  }

  for (cpt_month in 1:Number_Months)
  {
    diag(a_cor[cpt_month, , ]) <- 1
  }

  #                           Estimation of Realized Correlations                     #

  # Declaration of correlation matrix
  a_cor <- array(0, dim = c(Number_Months, Number_Institutions, Number_Institutions))

  # Loop to fullfill the covariance matrix
  for (cpt_month in 1:Number_Months)
  {
    for (cpt_institutions_1 in 1:Number_Institutions)
    {
      for (cpt_institutions_2 in 1:Number_Institutions)
      {
        a_cor[cpt_month, cpt_institutions_1, cpt_institutions_2] <- (a_covar[cpt_month, cpt_institutions_1, cpt_institutions_2]) / ((m_Realized_Volatility[cpt_month,cpt_institutions_1]) * (m_Realized_Volatility[cpt_month,cpt_institutions_2]))
      }
    }
  }

  for (cpt_month in 1:Number_Months)
  {
    diag(a_cor[cpt_month, , ]) <- 1
  }

  #                                   Correction PSD & NaN                            #

  a_cor[a_cor > 1] <- 0
  a_cor[a_cor < -1] <- 0

  a_cor[is.na(a_cor)] <- 0
  a_cor[is.nan(a_cor)] <- 0

  ##############################################################################################################################
  ##############################################################################################################################
  ##############################################################################################################################

  #                              Computing Correlation Networks                       #

  # Declaration of variables
  a_dist <- array(0, dim = c(Number_Months, Number_Institutions, Number_Institutions))
  m_distances <- matrix(0, ncol = Number_Institutions, nrow = Number_Institutions)
  m_dist <- matrix(0, ncol = Number_Institutions, nrow = Number_Institutions)
  Mean_Corr_Diff <- vector(mode="numeric", length=Number_Months)
  sum_corr <- vector(mode="numeric", length=Number_Months)
  Ratio_Indirect_Links <- vector(mode="numeric", length=Number_Months)
  m_tampon <- matrix(0, ncol = Number_Institutions, nrow = Number_Institutions)
  m_cor <- matrix(0, ncol = Number_Institutions, nrow = Number_Institutions)
  a_diff_corr <- array(0, dim = c(Number_Months, Number_Institutions, Number_Institutions))
  a_diff_corr_2 <- array(0, dim = c(Number_Months, Number_Institutions, Number_Institutions))
  a_corrected_corr <- array(0, dim = c(Number_Months, Number_Institutions, Number_Institutions))
  v_Degree <- vector(mode="numeric", length = Number_Months)
  v_Closeness_Centrality <- vector(mode="numeric", length = Number_Months)
  v_Eigenvector_Centrality <- vector(mode="numeric", length = Number_Months)

  v_norm_systemic_risk <- vector(mode="numeric", length = Number_Months)
  v_norm_degree <- vector(mode="numeric", length = Number_Months)
  v_norm_closeness_centrality <- vector(mode="numeric", length = Number_Months)
  v_norm_eigenvector_centrality <- vector(mode="numeric", length = Number_Months)

  # Loop building distance matrices, correlation matrices then Index Complex
  for (cpt_month in 1:Number_Months)
  {
    #	Transform Correlation Matrix a_cor into Distance Matrix m_dist     #

    # Prepare a N x N unity matrix to compute distance matrices
    m_1 <- matrix( rep(1, Number_Institutions * Number_Institutions), nrow = Number_Institutions, ncol = Number_Institutions)


    # distance = 1 - correlation
    m_tampon <- (- log ( a_cor[cpt_month, , ] * a_cor[cpt_month, , ] ) )
    m_tampon_1 <- round(m_tampon,4)
    a_dist[cpt_month, , ] <- m_tampon_1
    m_dist <- a_dist[cpt_month, , ]


    # The diagonal of a distance matrix must be equal to 0
    diag(m_dist) <- 0

    #	         Transform Distance Matrix m_dist into Graph ig            #

    m_dist[m_dist == Inf] <- 10000
    # Compute graph from distance matrix (the lattest is an adjacency matrix actually)
    ig <- igraph::graph_from_adjacency_matrix(m_dist, mode = "undirected", weighted = TRUE)

    #	         Transform Distance Matrix m_dist into Graph ig            #

    v_Closeness_Centrality[cpt_month] <- sum(igraph::closeness(ig, vids = igraph::V(ig), mode = "total", normalized = FALSE))
    v_Degree[cpt_month] <- sum(igraph::strength(ig))
    buffer_eigen_value <- igraph::eigen_centrality(ig, directed = FALSE)#, scale = TRUE)
    v_Eigenvector_Centrality[cpt_month] <- buffer_eigen_value$value


    #         Compute Shortest Paths Matrix m_distances from Graph ig        #

    tryCatch({
      m_distances <- igraph::distances(ig, v = igraph::V(ig), to = igraph::V(ig), algorithm = "bellman-ford")
    }, error=function(e){warning("Distance calculation error")})

    #                                                                        #

    # Here we should change the name of correlation matrices "m_dist_1" and "m_distances_1" because they're not distance matrices but correlation matrices




      # m_dist (distance matrix) is prepared to be compared
      m_dist[abs(m_dist) < 0.000001] <- 0
      m_dist_1 <- sqrt(exp( - m_dist))
      m_dist_1[!is.finite(m_dist_1)] <- 0
      diag(m_dist_1) <- 0

      # m_distances (shortest paths matrix) is prepared to be compared
      m_distances[abs(m_distances) < 0.000001] <- 0
      m_distances_1 <- sqrt(exp( - m_distances))
      m_distances_1[!is.finite(m_distances_1)] <- 0
      diag(m_distances_1) <- 0





    #       Compare correlation matrix and corrected correlation matrix      #

    m_corr_diff <- m_distances_1 - m_dist_1
    m_corr_diff[abs(m_corr_diff) < 0.000001] <- 0


    #                     Back to Comparison Matrices                        #

    # Load the correlation difference matrix in an array to exploit it later
    a_diff_corr[cpt_month,,] <- m_corr_diff
    diag(a_diff_corr[cpt_month,,]) <- 1

    a_corrected_corr[cpt_month,,] <- m_distances_1
    diag(a_corrected_corr[cpt_month,,]) <- 1


    #                            Prepare Results                             #

    # Compute the mean correlation differences
    tampon_m_corr_diff <- mean(m_corr_diff[m_corr_diff != 0])
    tampon_m_corr_diff[is.nan(tampon_m_corr_diff)] = 0
    Mean_Corr_Diff[cpt_month] <- tampon_m_corr_diff


  }

  #                   Results Exploitation : Part 2 - Risk Estimation                 #

  # Declare variables
  a_real_volat_port <- array(0, dim = c(Number_Months, Number_Institutions, Number_Institutions))
  real_volat_port <- vector(mode="numeric", length=Number_Months)
  real_volat_port_2 <- vector(mode="numeric", length=Number_Months)
  real_volat_port_corr <- vector(mode="numeric", length=Number_Months)
  real_volat_port_diff_corr <- vector(mode="numeric", length=Number_Months)
  v_Weights <- numeric(Number_Institutions)


  # Equal-weighted portfolio (index)
  Equal_Weights <- rep(1 / Number_Institutions, Number_Months)
  v_Weights <- Equal_Weights

  #      Compute the realized volatility of portfolio    #


  # Initialization of the realized volatility of portfolio
  real_volat_port[cpt_month] <- 0

  # Double Loop to compute the realized volatility of the portfolio
  for (cpt_month in 1:Number_Months)
  {
    for (cpt_institutions_1 in 1:Number_Institutions)
    {
      for (cpt_institutions_2 in 1:Number_Institutions)
      {
        real_volat_port[cpt_month] <- real_volat_port[cpt_month] + v_Weights[cpt_month] * v_Weights[cpt_month] * m_Realized_Volatility[cpt_month,cpt_institutions_1] * m_Realized_Volatility[cpt_month,cpt_institutions_2] * a_cor[cpt_month, cpt_institutions_1, cpt_institutions_2]
      }
    }

  }


  #         Compute the realized volatility of portfolio with CORRECTED CORRELATION       #


  # Initialization of the realized volatility of portfolio WITH CORRECTED CORRELATION
  real_volat_port_corr[cpt_month] <- 0

  # Double Loop to compute the realized volatility of the portfolio WITH CORRECTED CORRELATION
  for (cpt_month in 1:Number_Months)
  {
    for (cpt_institutions_1 in 1:Number_Institutions)
    {
      for (cpt_institutions_2 in 1:Number_Institutions)
      {
        real_volat_port_corr[cpt_month] <- real_volat_port_corr[cpt_month] + v_Weights[cpt_month] * v_Weights[cpt_month] * m_Realized_Volatility[cpt_month,cpt_institutions_1] * m_Realized_Volatility[cpt_month,cpt_institutions_2] * a_corrected_corr[cpt_month, cpt_institutions_1, cpt_institutions_2]
      }
    }
  }


  #                         Compute Systemic Risk                          #

  systemic_risk <- real_volat_port_corr - real_volat_port

  #                           Prepare results                           #

  # Normalization
  v_norm_systemic_risk <- f_scale(systemic_risk)
  v_norm_volat <- f_scale(real_volat_port)
  v_norm_degree <- f_scale(v_Degree)
  v_norm_closeness_centrality <- f_scale(v_Closeness_Centrality)
  v_norm_eigenvector_centrality <- f_scale(v_Eigenvector_Centrality)

  # Time scale
  time_scale <- seq(as.Date(df_data_returns[2,1], "%d/%m/%Y"), as.Date(df_data_returns[Number_Observations,1], "%d/%m/%Y"), by="months")
  time_scale_2 <- time_scale[2:(Number_Months+1)]

  # xts convert
  xts_SR <- xts::xts(v_norm_systemic_risk, order.by=as.Date(time_scale_2,"%Y/%m/%d"))
  xts_volat <- xts::xts(v_norm_volat , order.by=as.Date(time_scale_2,"%Y/%m/%d"))
  xts_degree <- xts::xts(v_norm_degree, order.by=as.Date(time_scale_2,"%Y/%m/%d"))
  xts_closeness_centrality <- xts::xts(v_norm_closeness_centrality, order.by=as.Date(time_scale_2,"%Y/%m/%d"))
  xts_eigenvector_centrality <- xts::xts(v_norm_eigenvector_centrality, order.by=as.Date(time_scale_2,"%Y/%m/%d"))

  # Result as a list
  result <- list("Degree" = xts_degree, "Closeness_Centrality" = xts_closeness_centrality, "Eigenvector_Centrality" = xts_eigenvector_centrality, "SR" = xts_SR, "Volatility" = xts_volat)
  return(result)
}

#####################################################################################
#                                      END                                          #
#####################################################################################
