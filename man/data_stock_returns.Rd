\name{data_stock_returns}
\alias{data_stock_returns}
\docType{data}
\title{
Financial institutions (banks, insurers and asset managers) stock returns
}
\description{
This dataset includes state variables data extracted from the FRED and Yahoo Finance. Specifically, it includes dates, MSCI STOXX Europe 600 Index returns and banks, insurers and asset managers stock returns.
}
\usage{data("data_stock_returns")}
\format{
  A data frame with 5030 observations on the following 74 variables.
  \describe{
    \item{\code{ACKB.BB.Equity}}{a numeric vector}
    \item{\code{AGN.NA.Equity}}{a numeric vector}
    \item{\code{AGS.BB.Equity}}{a numeric vector}
    \item{\code{AIBG.ID.Equity}}{a numeric vector}
    \item{\code{ALV.GY.Equity}}{a numeric vector}
    \item{\code{AV..LN.Equity}}{a numeric vector}
    \item{\code{BALN.SE.Equity}}{a numeric vector}
    \item{\code{BARC.LN.Equity}}{a numeric vector}
    \item{\code{BBVA.SQ.Equity}}{a numeric vector}
    \item{\code{BIRG.ID.Equity}}{a numeric vector}
    \item{\code{BKT.SQ.Equity}}{a numeric vector}
    \item{\code{BNP.FP.Equity}}{a numeric vector}
    \item{\code{BPE.IM.Equity}}{a numeric vector}
    \item{\code{CBG.LN.Equity}}{a numeric vector}
    \item{\code{CBK.GY.Equity}}{a numeric vector}
    \item{\code{CNP.FP.Equity}}{a numeric vector}
    \item{\code{CS.FP.Equity}}{a numeric vector}
    \item{\code{CSGN.SE.Equity}}{a numeric vector}
    \item{\code{DANSKE.DC.Equity}}{a numeric vector}
    \item{\code{DBK.GY.Equity}}{a numeric vector}
    \item{\code{DNB.NO.Equity}}{a numeric vector}
    \item{\code{Date}}{a date vector}
    \item{\code{EBS.AV.Equity}}{a numeric vector}
    \item{\code{EMG.LN.Equity}}{a numeric vector}
    \item{\code{G.IM.Equity}}{a numeric vector}
    \item{\code{GBLB.BB.Equity}}{a numeric vector}
    \item{\code{GLE.FP.Equity}}{a numeric vector}
    \item{\code{HELN.SE.Equity}}{a numeric vector}
    \item{\code{HNR1.GY.Equity}}{a numeric vector}
    \item{\code{HSBA.LN.Equity}}{a numeric vector}
    \item{\code{HSX.LN.Equity}}{a numeric vector}
    \item{\code{ICP.LN.Equity}}{a numeric vector}
    \item{\code{III.LN.Equity}}{a numeric vector}
    \item{\code{INDUA.SS.Equity}}{a numeric vector}
    \item{\code{INGA.NA.Equity}}{a numeric vector}
    \item{\code{INVEB.SS.Equity}}{a numeric vector}
    \item{\code{ISP.IM.Equity}}{a numeric vector}
    \item{\code{JYSK.DC.Equity}}{a numeric vector}
    \item{\code{KBC.BB.Equity}}{a numeric vector}
    \item{\code{KINVB.SS.Equity}}{a numeric vector}
    \item{\code{KN.FP.Equity}}{a numeric vector}
    \item{\code{KOMB.CK.Equity}}{a numeric vector}
    \item{\code{LGEN.LN.Equity}}{a numeric vector}
    \item{\code{LLOY.LN.Equity}}{a numeric vector}
    \item{\code{LUNDB.SS.Equity}}{a numeric vector}
    \item{\code{MAP.SQ.Equity}}{a numeric vector}
    \item{\code{MB.IM.Equity}}{a numeric vector}
    \item{\code{MF.FP.Equity}}{a numeric vector}
    \item{\code{MUV2.GY.Equity}}{a numeric vector}
    \item{\code{NDA.SS.Equity}}{a numeric vector}
    \item{\code{NXG.LN.Equity}}{a numeric vector}
    \item{\code{OML.LN.Equity}}{a numeric vector}
    \item{\code{PARG.SE.Equity}}{a numeric vector}
    \item{\code{PRU.LN.Equity}}{a numeric vector}
    \item{\code{RBS.LN.Equity}}{a numeric vector}
    \item{\code{RF.FP.Equity}}{a numeric vector}
    \item{\code{RSA.LN.Equity}}{a numeric vector}
    \item{\code{SAMPO.FH.Equity}}{a numeric vector}
    \item{\code{SAN.SQ.Equity}}{a numeric vector}
    \item{\code{SCR.FP.Equity}}{a numeric vector}
    \item{\code{SDR.LN.Equity}}{a numeric vector}
    \item{\code{SEBA.SS.Equity}}{a numeric vector}
    \item{\code{SHBA.SS.Equity}}{a numeric vector}
    \item{\code{SLHN.SE.Equity}}{a numeric vector}
    \item{\code{SREN.SE.Equity}}{a numeric vector}
    \item{\code{STAN.LN.Equity}}{a numeric vector}
    \item{\code{STB.NO.Equity}}{a numeric vector}
    \item{\code{STJ.LN.Equity}}{a numeric vector}
    \item{\code{SWEDA.SS.Equity}}{a numeric vector}
    \item{\code{SXXP.Index}}{a numeric vector}
    \item{\code{SYDB.DC.Equity}}{a numeric vector}
    \item{\code{UBSG.SE.Equity}}{a numeric vector}
    \item{\code{UCG.IM.Equity}}{a numeric vector}
    \item{\code{ZURN.SE.Equity}}{a numeric vector}




  }
}
\source{
Federal Reserve Economic Data (FRED) St. Louis Fed and Yahoo Finance
}
\references{
Hasse, Jean-Baptiste. "Systemic Risk: a Network Approach". AMSE Working Paper (2020)

}
\examples{
data("data_stock_returns")
head(data_stock_returns)
}
\keyword{datasets}
