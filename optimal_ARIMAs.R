library(tsoutliers)

best.ARIMA=function(series,search=2){
	n.col=5
	n.row=search^3
	df=data.frame(matrix(rep(0,n.col*n.row),n.row)) #pre-allocating in memory
	names(df)=c("p","d","q","AIC","VALIDITY")
	i=0
	for(d in 0:search){
		if(d!=0){
			adf.p=adf.test(diff(series,d))$p.value	# Dickey-Fuller stationarity test
		} else {
			adf.p=adf.test(series)$p.value
		}
		for(p in 0:search){
			for(q in 0:search){
				i=i+1
				skip_to_next <- FALSE
				tryCatch(print(arima(series, order = c(p,d,q), method="ML")), # estimation errors not breaking loop
								 error = function(e) { skip_to_next <<- TRUE})  
				# problems: we have to estimate every model TWICE, and each model will be printed the first time it's estimated
				if(skip_to_next) { next }
				mod=arima(series, order = c(p,d,q), method="ML")   # fitting model and saving it in memory
				jb.p=jarque.bera.test(mod$residuals)$p.value    # Jarque-Bera test for normality of residuals
				box.p=Box.test(mod$residuals, type="Ljung-Box")$p.value     # Box-Ljung test for independence of residuals
				df[i,]=c(p,d,q,mod$aic,
								 (sum(pnorm(c(abs(mod$coef)/sqrt(diag(mod$var.coef))),
								 					 mean=0, sd=1, lower.tail=FALSE)>0.05)==0)*(jb.p>0.05)*(box.p>0.05)*(adf.p<=0.05)) 
				# list with p, d, q, AIC and wether the model satisfies the necessary assumptions at the 95% confidence (not corrected)
			}
		}
	}
	df=df[df$VALIDITY==1,] # get valid models only
	df=na.omit(df[order(df$AIC),]) # get best valid model
	return(df[1,]) # the function outputs the hyperparameters of the valid model with highest Akaike Information Criteria
}

Stocks <- read.csv("~/Datasets/StocksT.csv")
stocks=ts(Stocks[,-1])
parameters=data.frame(matrix(rep(0,5*nrow(stocks)),nrow(stocks)))
names(parameters)=c("p","d","q","AIC","VALIDITY")
for(i in 1:ncol(stocks)){
	parameters[i,]=suppressWarnings(best.ARIMA(stocks[,i],search=3))  # applying best.ARIMA to all stocks
}
parameters=cbind("ticker"=names(Stocks)[-1],parameters)
View(parameters)   # results
