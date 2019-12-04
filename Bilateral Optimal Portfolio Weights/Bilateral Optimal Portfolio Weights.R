
### Kroner, KF., and NG, VK (1998)
### MODELING ASYMMETRIC COMOVEMENTS OF ASSET RETURNS
### REVIEW OF FINANCIAL STUDIES
### replicated by David Gabauer (https://sites.google.com/view/davidgabauer/contact-details)

PortfolioWeights = function(data, H, type="long") {
   k=ncol(data)
   NAMES=matrix(NA, ncol=k, nrow=k)
   for (i in 1:k){
      NAMES[i,]=paste0(colnames(data),"/",colnames(data)[i])
   }
   lowpre=NAMES[lower.tri(diag(k))]
   
   values = function(x){
      x = as.matrix(x)
      res = matrix(NA, nrow=ncol(x), ncol=4)
      for (i in 1:ncol(x)){
         res[i,]=matrix(c(mean(x[,i]), sd(x[,i]), quantile(x[,i],0.05), quantile(x[,i],0.95)), nrow=1)
      }
      colnames(res)=c("Mean","Std.Dev.","5%","95%")
      res
   }
   col=sum(lower.tri(diag(k)))
   sel=which(lower.tri(diag(ncol(data)))==TRUE, arr.ind=T)
   
   summary = NULL
   portfolio_weights = array(1, c(k, k, (nrow(data))))
   for (i in 1:k){
      for (j in 1:k){
         if (i==j) {
            pwpre = c(10,10,10,10)
            summary = rbind(summary,pwpre)
         } else {
            RR = (H[j,j,]-H[i,j,])/(H[i,i,]-2*H[i,j,]+H[j,j,])
            if (type=="long") {
               RR = ifelse(RR>1,1,RR)
               RR = ifelse(RR<0,0,RR)
            }
            pwpre = values(RR)
            summary = rbind(summary,pwpre)
            colnames(summary)=c("Mean","Std.Dev.","5%","95%")
            portfolio_weights[i,j,] = RR
         }
      }
   }
   rownames(summary) = c(NAMES)
   
   pval = cumulative_ratio = HE = matrix(NA,ncol=k,nrow=k)
   rownames(cumulative_ratio)=colnames(cumulative_ratio)=rownames(HE)=colnames(HE)=colnames(data)
   portfolio_return = cumulative_portfolio_return = cumulative_asset_return = array(NA,c(k,k,nrow(data)))
   for (i in 1:k) {
      for (j in 1:k) {
         portfolio_return[i,j,] = portfolio_weights[i,j,]*data[,i] + (1-portfolio_weights[i,j,])*data[,j]
         HE[i,j] = 1 - var(portfolio_return[i,j,])/var(data[,i])
         pval[i,j] = var.test(x=portfolio_return[i,j,],y=data[,i],ratio=1)$p.value
         cumulative_asset_return[i,j,] = cumsum(data[,i])
         cumulative_portfolio_return[i,j,] = cumsum(portfolio_return[i,j,])
         cumulative_ratio[i,j] = cumulative_portfolio_return[i,j,dim(cumulative_portfolio_return)[3]]/cumulative_asset_return[i,j,dim(cumulative_asset_return)[3]]
      }
   }
   table = round(cbind(summary,c(t(HE)),c(t(pval))),2)
   rownames(table)=c(NAMES)
   table = table[-which(table[,1]==10),]
   colnames(table)=c("Mean","Std.Dev.","5%","95%","HE","p-value")
   
   return = list(portfolio_weights=portfolio_weights, summary=format(round(table,2),nsmall=2), portfolio_return=portfolio_return, cumulative_portfolio_return=cumulative_portfolio_return, cumulative_asset_return=cumulative_asset_return, hedging_effectiveness=HE, pvalue=pval)
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

### BILATERAL PORTFOLIO WEIGHTS
bpw = PortfolioWeights(data, H)
bpw$summary
bpw$hedging_effectivenss # variance reduction (Ederington, 1979)
bpw$pvalue # significance level (alpha), (Antonakakis et al., 2018)

par(mfcol = c(ceiling((k-1)*k/4), 2), oma = c(.5,.5,0,0) + 0.1, mar = c(.5,.5,0,0) + 1, mgp=c(.5, .5, 0))
for (i in 1:k) {
   for (j in 1:k) {
      if (j<i) {
         plot(date, bpw$portfolio_weights[i,j,], type="l", col="steelblue4", las=1, ylim=c(0, max(bpw$portfolio_weights)), main=paste0(colnames(data)[i], "/", colnames(data)[j]), xlab="", ylab="", xaxs="i", yaxs="i")
         grid(); box()
         lines(date, bpw$portfolio_weights[i,j,], col="steelblue4")
         lines(date, bpw$portfolio_weights[j,i,], col="steelblue1")
         legend("topleft",c(paste0(colnames(data)[i], "/", colnames(data)[j]), paste0(colnames(data)[j], "/", colnames(data)[i])), fill=c("steelblue4","steelblue1"), bty='n',)
      }
   }
}

par(mfcol = c(ceiling((k-1)*k/4), 2), oma = c(.5,.5,0,0) + 0.1, mar = c(.5,.5,0,0) + 1, mgp=c(.5, .5, 0))
for (i in 1:k) {
   for (j in 1:k) {
      if (j<i) {
         plot(date, bpw$cumulative_asset_return[i,i,], type="l", col="grey", las=1, ylim=c(min(hr$cumulative_portfolio_return), max(hr$cumulative_asset_return)), main=paste0(colnames(data)[i], "/", colnames(data)[j]), xlab="", ylab="", xaxs="i", yaxs="i")
         lines(date, bpw$cumulative_asset_return[j,j,],col="grey", lty=3)
         lines(date, bpw$cumulative_portfolio_return[i,j,], col="steelblue4")
         grid(); box()
         abline(h=0, lty=3)
      }
   }
}

### SHARPE RATIO
apply(data, 2, mean) / apply(data, 2, sd) # individual assets
apply(bpw$portfolio_return, 2, mean) / apply(bpw$portfolio_return, 2, sd)

### END
