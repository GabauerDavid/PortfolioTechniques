
### ENGLE, RF. (2016)
### DYNAMIC CONDITIONAL BETA
### JOURNAL OF FINANCIAL ECONOMETRICS

# KRONER, KF., AND SULTAN, J. (1993)
# TIME-VARYING DISTRIBUTIONS AND DYNAMIC HEDGING WITH FOREIGN CURRENCY FUTURES
# JOURNAL OF FINANCIAL AND QUANTITATIVE ANALYSIS
### replicated by David Gabauer (https://sites.google.com/view/davidgabauer/contact-details)

HedgeRatio = function(data, H) {
   k=ncol(data)
   NAMES = matrix(NA, ncol=k, nrow=k)
   for (i in 1:k){
      NAMES[i,]=paste0(colnames(data),"/",colnames(data)[i])
   }
   lowpre = NAMES[lower.tri(diag(k))]
   
   DM = array(NA, c(k, k, (nrow(data))))
   for (i in 1:k){
      for (j in 1:k){
         DM[i,j,] = H[i,j,] / H[j,j,]
      }
   }
   
   values = function(x){
      x = as.matrix(x)
      res = matrix(NA, nrow=ncol(x), ncol=4)
      for (i in 1:ncol(x)){
         res[i,]=matrix(c(mean(x[,i]), sd(x[,i]), quantile(x[,i],0.05), quantile(x[,i],0.95)), nrow=1)
      }
      colnames(res)=c("Mean","Std.Dev.","5%","95%")
      res
   }
   
   HR = array(1, c(k, 4, k))
   for (i in 1:k){
      for (j in 1:k){
         HR[i,,j] = values(DM[i,j,])
      }
   }
   colnames(HR) = c("Mean","Std.Dev.","5%","95%")
   dimnames(HR)[[1]] = colnames(data)
   
   for (i in 1:k){
      if (i==1){
         HRatio=HR[,,i]
      } else {
         HRatio=rbind(HRatio, HR[,,i])
      }
   }
   rownames(HRatio) = c(t(NAMES))
   
   pval = HE = matrix(NA,ncol=k,nrow=k)
   colnames(HE)=rownames(HE)=colnames(pval)=rownames(pval)=colnames(data)
   portfolio_return = cumulative_portfolio_return = cumulative_asset_return = array(NA,c(k,k,nrow(data)))
   cumulative_ratio = array(NA,c(k,k))
   for (i in 1:k) {
      for (j in 1:k) {
         portfolio_return[i,j,] = data[,i] - DM[i,j,]*data[,j]
         HE[i,j] = 1 - var(portfolio_return[i,j,])/var(data[,i])
         pval[i,j] = var.test(x=portfolio_return[i,j,],y=data[,i],ratio=1)$p.value
         cumulative_asset_return[i,j,] = cumsum(data[,i])
         cumulative_portfolio_return[i,j,] = cumsum(portfolio_return[i,j,])
         cumulative_ratio[i,j] = cumulative_portfolio_return[i, j, dim(cumulative_portfolio_return)[3]] / cumulative_asset_return[i,j,dim(cumulative_portfolio_return)[3]]
      }
   }
   table=cbind(round(HRatio,2),c(t(HE)),c(t(pval)))
   table=table[-which(table[,1]==1),]
   colnames(table)=c("Mean","Std.Dev.","5%","95%","HE","p-value")
   
   return = list(hedge_ratio=DM, summary=format(round(table,2),nsmall=2), portfolio_return=portfolio_return, cumulative_portfolio_return=cumulative_portfolio_return, cumulative_asset_return=cumulative_asset_return, hedging_effectiveness=HE, pvalue=pval)
}

library("quantmod")
TICKERS = c("^GSPC,T,F,HON")
TICKERS = unlist(strsplit(TICKERS,","))
START = "2010-01-01"
FREQUENCY = "daily"

getSymbols(c(TICKERS), from=START, periodicity=FREQUENCY)
ind = which(substr(TICKERS,1,1)=="^")
if (length(ind)>0) {
   TICKERS[ind] = substr(TICKERS[ind],2,10)
} else {
   TICKERS = TICKERS
}

data = NULL
for (i in 1:length(TICKERS)) {
   data = cbind(data, 100 * diff(log(get(TICKERS[i])[,1])))
}
data = na.omit(data)
date = index(data)
colnames(data) = unlist(strsplit(colnames(data), ".Open"))
k = ncol(data)

library("rmgarch")
ugarch = ugarchspec(variance.model=list(garchOrder=c(1, 1), model="sGARCH"),
                    mean.model=list(armaOrder=c(0, 0)),
                    distribution.model="norm")
mgarch = multispec( replicate(k, ugarch) )
dccgarch_spec = cgarchspec(uspec=mgarch, dccOrder=c(1,1), asymmetric=FALSE,
                           distribution.model=list(copula="mvt", method="Kendall", time.varying=TRUE, transformation="parametric"))
dcc_fit = cgarchfit(dccgarch_spec, data=data, solver="solnp", fit.control=list(eval.se=TRUE) )
H = rcov(dcc_fit)

### HEDGE RATIOS / CAPM BETA
hr = HedgeRatio(data, H)
hr$summary # FYI: hedge ratios across stock indices and assets are considered as CAPM betas
hr$hedging_effectiveness # variance reduction (Ederington, 1979)
hr$pvalue # significance level (alpha), (Antonakakis et al., 2018)

par(mfcol = c(ceiling((k-1)*k/4), 2), oma = c(.5,.5,0,0) + 0.1, mar = c(.5,.5,0,0) + 1, mgp=c(.5, .5, 0))
for (i in 1:k) {
   for (j in 1:k) {
      if (j<i) {
         plot(date, hr$hedge_ratio[i,j,], type="l", col="steelblue4", las=1, ylim=c(min(hr$hedge_ratio), max(hr$hedge_ratio)), main=paste0(colnames(data)[i], "/", colnames(data)[j]), xlab="", ylab="", xaxs="i", yaxs="i")
         grid(); box()
         lines(date, hr$hedge_ratio[i,j,], col="steelblue4")
         lines(date, hr$hedge_ratio[j,i,], col="steelblue1")
         legend("topleft",c(paste0(colnames(data)[i], "/", colnames(data)[j]), paste0(colnames(data)[j], "/", colnames(data)[i])), fill=c("steelblue4","steelblue1"), bty='n',)
      }
   }
}

par(mfcol = c(ceiling((k-1)*k/4), 2), oma = c(.5,.5,0,0) + 0.1, mar = c(.5,.5,0,0) + 1, mgp=c(.5, .5, 0))
for (i in 1:k) {
   for (j in 1:k) {
      if (j<i) {
         plot(date, hr$cumulative_asset_return[i,i,], type="l", col="grey", las=1, ylim=c(min(hr$cumulative_portfolio_return), max(hr$cumulative_asset_return)), main=paste0(colnames(data)[i], "/", colnames(data)[j]), xlab="", ylab="", xaxs="i", yaxs="i")
         lines(date, hr$cumulative_asset_return[j,j,],col="grey", lty=3)
         lines(date, hr$cumulative_portfolio_return[i,j,], col="steelblue4")
         lines(date, hr$cumulative_portfolio_return[j,i,], col="steelblue1")
         grid(); box()
         abline(h=0, lty=3)
         legend("topleft",c(paste0(colnames(data)[i], "/", colnames(data)[j]), paste0(colnames(data)[j], "/", colnames(data)[i])), fill=c("steelblue4","steelblue1"), bty='n',)
      }
   }
}

### SHARPE RATIO
apply(data, 2, mean) / apply(data, 2, sd) # individual assets
apply(hr$portfolio_return, 1:2, mean) / apply(hr$portfolio_return, 1:2, sd) # hedged positions

### END
