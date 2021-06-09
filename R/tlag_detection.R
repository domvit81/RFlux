
tlag_detection <- function (x, y, mfreq, x.model = ar.res, AIC=FALSE, LAG.MAX=mfreq*10, show.plot=FALSE) 
{
	
	local_max <- function(x) which(x - shift(x, 1) > 0  & x - shift(x, 1, type='lead') > 0)
	local_min <- function(x) which(x - shift(x, 1) < 0  & x - shift(x, 1, type='lead') < 0)
	## https://stackoverflow.com/questions/34205515/finding-local-maxima-and-minima-in-r


	cross_cov <- ccf(x = x, y = y, lag.max=LAG.MAX, plot=FALSE, type="covariance")
	scross_cov <- rollapply(cross_cov$acf, width=3, FUN="mean", fill=NA) ## smoothed cross-covariance function
	
    filter.mod = function(x, model) {
    	x <- x - mean(x, na.rm=TRUE)
        if (length(model$Delta) >= 1) 
            x = stats::filter(x, filter = c(1, -model$Delta), method = "convolution", sides = 1)
        if (length(model$theta) >= 1 && any(model$theta != 0)) 
            x = stats::filter(na.omit(x), filter = -model$theta, method = "recursive", sides = 1)
        if (length(model$phi) >= 1 && any(model$phi != 0)) 
            x = stats::filter(x, filter = c(1, -model$phi), method = "convolution", sides = 1)
        x
    }
    if (!missing(x.model)) {
        x = filter.mod(x, model = x.model$model)
        y = filter.mod(y, model = x.model$model)
    }
    else {
        if (AIC==TRUE) ar.res = ar.ols(x, aic=TRUE)
        if (AIC==FALSE) ar.res = ar.ols(x, aic=FALSE, order.max=10*log(length(x),10))
        x = stats::filter(x, filter = c(1, -ar.res$ar), method = "convolution", sides = 1)
        y = stats::filter(y, filter = c(1, -ar.res$ar), method = "convolution", sides = 1)
    }
    cross_cor = ccf(x = x, y = y, na.action = na.omit, plot=FALSE, lag.max=LAG.MAX)
    
    scross_cor <- rollapply(cross_cor$acf, width=13, FUN="mean", fill=NA) ## smoothed cross-correlation function 
  
	tl0 <- which.max(abs(scross_cor))
 
    tlag_pw <- tl0 - LAG.MAX - 1
    tlag_max <- local_max(scross_cov[max(1,tl0-12):min(tl0+12, LAG.MAX*2)]) + tl0 - 13 - LAG.MAX - 1
    tlag_min <- local_min(scross_cov[max(1,tl0-12):min(tl0+12, LAG.MAX*2)]) + tl0 - 13 - LAG.MAX - 1
    
   if(length(tlag_max)==0 & length(tlag_min)==0)  tlag_opt <- tlag_pw
   if(length(tlag_max)==1 & length(tlag_min)==0) tlag_opt <- tlag_max 
   if(length(tlag_max)==0 & length(tlag_min)==1) tlag_opt <- tlag_min
   if(length(tlag_max) >= 1 | length(tlag_min) >= 1) tlag_opt <- tlag_pw
    corr_est <- cross_cor$acf[tlag_opt+LAG.MAX+1] 
   
    if (show.plot==TRUE) {
	par(mfrow=c(2,1), mar=c(4,4,2,1), oma=c(1,1,1,1), las=0)
	plot(cross_cor$lag, cross_cor$acf, ylab="cross-correlation", xlab="Lag (timesteps)", type="h", col="grey68")    
    lines(cross_cor$lag, scross_cor, col=1, lwd=2)
    abline(h=1.96/sqrt(length(x)), col=4, lty=3)
    abline(h=-1.96/sqrt(length(x)), col=4, lty=3)
    abline(h=2.57/sqrt(length(x)), col=4, lty=2)
    abline(h=-2.57/sqrt(length(x)), col=4, lty=2)
    abline(v=tlag_opt, col="red", lty=3)
 
   	if (abs(corr_est)>=2.57/sqrt(length(x))) mtext(side=3, adj=0, paste0("Detected peak at ", tlag_opt, " timesteps stat. sign. at 0.01 level")) 
    if (abs(corr_est)<2.57/sqrt(length(x)) && abs(corr_est)>1.96/sqrt(length(x))) mtext(side=3, adj=0, paste0("Detected peak at ", tlag_opt, " timesteps stat. sign. at 0.05 level")) 
    if (abs(corr_est)<1.96/sqrt(length(x))) mtext(side=3, adj=1, paste0("Detected peak not stat. significant")) 
	
	plot(cross_cov$lag, cross_cov$acf, ylab="cross-covariance", xlab="Lag (timesteps)", type="h", col="grey68")    
    lines(cross_cov$lag, scross_cov, col="cyan")
    abline(v=tlag_opt, col="red", lty=3)
    mtext(side=3, adj=0, paste0("Optimal Time Lag at ", tlag_opt, " timesteps")) 
    }
 
    return(list(ccf = cross_cor, model = x.model, "opt_tlag" = tlag_opt, "tlag_pw"= tlag_pw, "tlag_lmax" = tlag_max, "tlag_lmin"= tlag_min, "corr_est" = corr_est, "cv5pct"= 1.96/sqrt(length(x)),"cv1pct"= 2.57/sqrt(length(x))))
}
