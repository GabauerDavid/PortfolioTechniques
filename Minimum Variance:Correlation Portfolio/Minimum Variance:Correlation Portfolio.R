
### CHRISTOFFERSEN, P., ERRUNZA, V., JACOBS, K., AND JIN, X. (2014)
### CORRELATION DYNAMICS AND INTERNATIONAL DIVERSIFICATION BENEFITS
### INTERNATIONAL JOURNAL OF FORECASTING

# MARKOWITZ, H. (1952) 
# PORTFOLIO SELECTION
# JOURNAL OF FINANCE
### replicated by David Gabauer (https://sites.google.com/view/davidgabauer/contact-details)

MinimumVariancePortfolio = function(data, H, type="long") {
   k = ncol(data)
   n = nrow(data)
   I = matrix(1,k,1)
   
   values = function(x){
      x = as.matrix(x)
      res = matrix(NA, nrow=ncol(x), ncol=4)
      for (i in 1:ncol(x)){
         res[i,]=matrix(c(mean(x[,i]), sd(x[,i]), quantile(x[,i],0.05), quantile(x[,i],0.95)), nrow=1)
      }
      colnames(res)=c("Mean","Std.Dev.","5%","95%")
      res
   }
   
   portfolio_weights = matrix(NA,nrow=n,ncol=k)
   for (i in 1:n) {
      Vinv = solve(H[,,i])
      portfolio_weights[i,] = (Vinv%*%I) / c(t(I)%*%Vinv%*%I)
      if (type=="long") {
         portfolio_weights[i,] = ifelse(portfolio_weights[i,]<0,0,portfolio_weights[i,])
         portfolio_weights[i,] = ifelse(portfolio_weights[i,]>1,1,portfolio_weights[i,])
         portfolio_weights[i,] = portfolio_weights[i,] / sum(portfolio_weights[i,])
      }
   }
   summary = NULL
   for (i in 1:k) {
      summary = rbind(summary, values(portfolio_weights[,i]))
   }
   rownames(summary)=colnames(data)
   
   cumulative_ratio = f_statistics = pval = HE = matrix(NA,ncol=1,nrow=k)
   rownames(cumulative_ratio)=rownames(f_statistics)=rownames(pval)=rownames(HE)=colnames(data)
   cumulative_asset_return = matrix(NA,ncol=k,nrow=n)
   portfolio_return = matrix(NA,ncol=1,nrow=n)
   for (i in 1:n) {
      portfolio_return[i,] = sum(portfolio_weights[i,]*as.numeric(data[i,]))
   }
   cumulative_portfolio_return = cumsum(portfolio_return)
   
   for (i in 1:k) {
      HE[i,] = 1 - var(portfolio_return)/var(data[,i])
      xxx = rbind(cbind(portfolio_return,1),cbind(data[,i],2))
      dfpval = data.frame(val=xxx[,1],group=as.character(xxx[,2]))
      pval[i,] = fligner.test(val~group,dfpval)$p.val
      #    pval[i,] = var.test(x=portfolio_return,y=data[,i],ratio=1)$p.value
      cumulative_asset_return[,i] = cumsum(data[,i])
      cumulative_ratio[i,] = cumulative_portfolio_return[length(cumulative_portfolio_return)]/cumulative_asset_return[nrow(cumulative_asset_return),i]
   }
   table = cbind(summary,HE,pval)
   colnames(table)=c("Mean","Std.Dev.","5%","95%","HE","p-value")
   
   return = list(portfolio_weights=portfolio_weights, summary=format(round(table,2),nsmall=2), hedging_effectiveness=HE, pvalue=pval, portfolio_return=portfolio_return, cumulative_portfolio_return=cumulative_portfolio_return, cumulative_asset_return=cumulative_asset_return)
}

library("quantmod")
TICKERS = c("^GSPC,^GDAXI,^FTSE,^N225")
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
R = rcor(dcc_fit)

### MINIMUM VARIANCE PORTFOLIO
mvp = MinimumVariancePortfolio(data, H)
mvp$summary
mvp$hedging_effectiveness # variance reduction (Ederington, 1979)
mvp$pvalue # significance level (alpha), (Antonakakis et al., 2018)

### MINIMUM CORRELARION PORTFOLIO
mcp = MinimumVariancePortfolio(data, R)
mcp$summary
mcp$hedging_effectiveness # variance reduction (Ederington, 1979)
mcp$pvalue # significance level (alpha), (Antonakakis et al., 2018)

par(mfcol = c(ceiling(k/2), 2), oma = c(.5,.5,0,0) + 0.1, mar = c(.5,.5,0,0) + 1, mgp=c(.5, .5, 0))
for (i in 1:k) {
   plot(date, mvp$portfolio_weights[,i], type="l", col="steelblue4", las=1, ylim=c(0, max(mvp$portfolio_weights)), main=paste(colnames(data)[i], "weight"), xlab="", ylab="", xaxs="i", yaxs="i")
   grid(); box()
   lines(date, mvp$portfolio_weights[,i], col="steelblue4")
   lines(date, mcp$portfolio_weights[,i], col="steelblue1")
}

par(mfrow=c(1,1))
plot(date, mvp$cumulative_portfolio_return, type="l", col="steelblue4", las=1, ylim=c(min(mvp$cumulative_asset_return), max(mvp$cumulative_asset_return)), main="", xlab="", ylab="", xaxs="i", yaxs="i")
for (i in 1:k) {
   lines(date, mvp$cumulative_asset_return[,i], col="grey", lty=i)
   grid(); box()
   abline(h=0, lty=3)
}
lines(date, mvp$cumulative_portfolio_return, col="steelblue4")
lines(date, mcp$cumulative_portfolio_return, col="steelblue1")

### SHARPE RATIO
apply(data, 2, mean) / apply(data, 2, sd) # individual assets
apply(mvp$portfolio_return, 2, mean) / apply(mvp$portfolio_return, 2, sd) # minimum variance portfolio
apply(mcp$portfolio_return, 2, mean) / apply(mcp$portfolio_return, 2, sd) # minimum correlation portfolio

### END
