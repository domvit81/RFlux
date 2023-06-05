
tlag_detection <- function (scalar_var, tsonic_var, w_var, mfreq, model = "ar", AIC=TRUE, method=c("pearson", "spearman", "boot"), LAG.MAX=mfreq*10, j=100, Rboot, cM, plot.it=FALSE,...) 
{
	local_max <- function(x) which(x - shift(x, 1) > 0 & x - shift(x, 1, type='lead') > 0) 								   
	local_min <- function(x) which(x - shift(x, 1) < 0  & x - shift(x, 1, type='lead') < 0) 
	## https://stackoverflow.com/questions/34205515/finding-local-maxima-and-minima-in-r
	
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

	set <- na.omit(cbind(na.approx(scalar_var, na.rm=FALSE), na.approx(tsonic_var, na.rm=FALSE), na.approx(w_var, na.rm=FALSE)))
	kid_w <- kurtosis(diff(set[,3]), na.rm=TRUE)
    	kid_ts <- kurtosis(diff(set[,2]), na.rm=TRUE)
    	kid_c <- kurtosis(diff(set[,1]), na.rm=TRUE)
	
	if (model=="ar"){
	if (egcm::bvr.test(set[,1])$p.val<0.01 & bvr.test(set[,2])$p.val<0.01 & bvr.test(set[,3])$p.val<0.01) {
	x <- set[,1];
	y <- set[,2];
	z <- set[,3]
	};
	if (egcm::bvr.test(set[,1])$p.val>=0.01 | bvr.test(set[,2])$p.val>=0.01 | bvr.test(set[,3])$p.val>=0.01) {
	x <- diff(set[,1]);
	y <- diff(set[,2]);
	z <- diff(set[,3])
	};
	if (AIC==TRUE) ar.resx = ar(x, aic=TRUE, order.max=floor(j*log10(length(x))));
    	if (AIC==FALSE) ar.resx = ar(x, aic=FALSE, order.max=floor(j*log10(length(x))));
	if (AIC==TRUE) ar.resy = ar(y, aic=TRUE, order.max=floor(j*log10(length(x))));
    	if (AIC==FALSE) ar.resy = ar(y, aic=FALSE, order.max=floor(j*log10(length(x))));
	if (AIC==TRUE) ar.resz = ar(z, aic=TRUE, order.max=floor(j*log10(length(x))));
    	if (AIC==FALSE) ar.resz = ar(z, aic=FALSE, order.max=floor(j*log10(length(x))));
    	x1 = stats::filter(x, filter = c(1, -ar.resx$ar), method = "convolution", sides = 1);
    	y1 = stats::filter(y, filter = c(1, -ar.resx$ar), method = "convolution", sides = 1);
    	z1 = stats::filter(z, filter = c(1, -ar.resx$ar), method = "convolution", sides = 1);
    	x2 = stats::filter(x, filter = c(1, -ar.resy$ar), method = "convolution", sides = 1);
    	y2 = stats::filter(y, filter = c(1, -ar.resy$ar), method = "convolution", sides = 1)
    	z2 = stats::filter(z, filter = c(1, -ar.resy$ar), method = "convolution", sides = 1)
    	x3 = stats::filter(x, filter = c(1, -ar.resz$ar), method = "convolution", sides = 1);
    	y3 = stats::filter(y, filter = c(1, -ar.resz$ar), method = "convolution", sides = 1);
    	z3 = stats::filter(z, filter = c(1, -ar.resz$ar), method = "convolution", sides = 1);
   	}


	if (model=="arima"){
	if (egcm::bvr.test(set[,1])$p.val<0.01 & bvr.test(set[,2])$p.val<0.01 & bvr.test(set[,3])$p.val<0.01) {
	x <- set[,1];
	y <- set[,2];
	z <- set[,3]
	};
	if (egcm::bvr.test(set[,1])$p.val>=0.01 | bvr.test(set[,2])$p.val>=0.01 | bvr.test(set[,3])$p.val>=0.01) {
	x <- diff(set[,1]);
	y <- diff(set[,2]);
	z <- diff(set[,3])
	};
	mod1 <- forecast::auto.arima(x, d=0);
	mod2 <- forecast::auto.arima(y, d=0);
	mod3 <- forecast::auto.arima(z, d=0);
	x1 = filter.mod(x, model = mod1$model);
	y1 = filter.mod(y, model = mod1$model);
	z1 = filter.mod(z, model = mod1$model);
	x2 = filter.mod(x, model = mod2$model);
 	y2 = filter.mod(y, model = mod2$model);
    	x3 = filter.mod(x, model = mod3$model);
 	z3 = filter.mod(z, model = mod3$model)
 	}
    
## Cross-Covariance Function    
	cross_cov_w <- as.vector(ccf(x = scalar_var, y = w_var, lag.max=LAG.MAX, plot=FALSE, type="covariance", na.action=na.pass)$acf)
	cross_cov_ts <- as.vector(ccf(x = scalar_var, y = tsonic_var, lag.max=LAG.MAX, plot=FALSE, type="covariance", na.action=na.pass)$acf)
       
## PW - Pearson
    if (method=="pearson") {
    	ccf_ct <- as.vector(ccf(x = x1, y = y1, na.action = na.omit, plot=FALSE, lag.max=LAG.MAX)$acf);
    	ccfs_ct <- as.vector(zoo::na.locf(zoo::na.locf(rollapply(ccf_ct, width=13, FUN="mean", fill=NA) , na.rm=FALSE), fromLast=TRUE)); ## smoothed cross-correlation function 
    	ccf_cw <- as.vector(ccf(x = x1, y = z1, na.action = na.omit, plot=FALSE, lag.max=LAG.MAX)$acf);
    	ccfs_cw <- as.vector(zoo::na.locf(zoo::na.locf(rollapply(ccf_cw, width=13, FUN="mean", fill=NA) , na.rm=FALSE), fromLast=TRUE)); ## smoothed cross-correlation function 
    	ccf_tc <- as.vector(ccf(x = x2, y = y2, na.action = na.omit, plot=FALSE, lag.max=LAG.MAX)$acf);
    	ccfs_tc <- as.vector(zoo::na.locf(zoo::na.locf(rollapply(ccf_tc, width=13, FUN="mean", fill=NA) , na.rm=FALSE), fromLast=TRUE)); ## smoothed cross-correlation function 
    	ccf_wc <- as.vector(ccf(x = x3, y = z3, na.action = na.omit, plot=FALSE, lag.max=LAG.MAX)$acf);
    	ccfs_wc <- as.vector(zoo::na.locf(zoo::na.locf(rollapply(ccf_wc, width=13, FUN="mean", fill=NA) , na.rm=FALSE), fromLast=TRUE)) ## smoothed cross-correlation function 
	}
 
 ## PW - Spearman
    if (method=="spearman") {
    xrank <- rank(x);
    yrank <- rank(y, ties.method="average");
    attributes(xrank) <- attributes(x);
    attributes(yrank) <- attributes(y);
    cross_cor <- as.vector(ccf(xrank, yrank, na.action = na.omit, plot=FALSE, lag.max=LAG.MAX)$acf);
    cross_cor_s <- NULL;
    scross_cor <- as.vector(zoo::na.locf(zoo::na.locf(rollapply(cross_cor, width=13, FUN="mean", fill=NA) , na.rm=FALSE), fromLast=TRUE)); ## smoothed cross-correlation function 
    scross_cor3 <- as.vector(zoo::na.locf(zoo::na.locf(rollapply(cross_cor, width=3, FUN="mean", fill=NA) , na.rm=FALSE), fromLast=TRUE)); ## smoothed cross-correlation function 
    scross_cor5 <- as.vector(zoo::na.locf(zoo::na.locf(rollapply(cross_cor, width=5, FUN="mean", fill=NA) , na.rm=FALSE), fromLast=TRUE)) ## smoothed cross-correlation function 
   }

## PW - BOOTSTRAPPING (pearson)
 if (method=="boot") {
	bootccf_ct <- tsboot(cbind(x1,y1), function(x) ccf(x[,1],x[,2],na.action=na.pass, plot=FALSE, lag.max=LAG.MAX)$acf, R=Rboot, sim="fixed", l=LAG.MAX*2, parallel="snow");
	ccf_ct <- unlist(apply(bootccf_ct$t, MARGIN=2, function(x) mean(x, na.rm=TRUE)));
	ccfs_ct <- unlist(zoo::na.locf(zoo::na.locf(rollapply(ccf_ct, width=13, FUN="mean", fill=NA) , na.rm=FALSE), fromLast=TRUE));
 	ccf_ct1 <- unlist(apply(abs(bootccf_ct$t), MARGIN=1, function(x) which.max(x)));
 	ccfs_ct1 <- unlist(apply(bootccf_ct$t, MARGIN=1, function(x) which.max(abs(as.vector(zoo::na.locf(zoo::na.locf(rollapply(x, width=13, FUN="mean", fill=NA), na.rm=FALSE), fromLast=TRUE))))));
	hdi_ct <- HDInterval::hdi(ccf_ct1, credMass=cM);
 	hdis_ct <- HDInterval::hdi(ccfs_ct1, credMass=cM);
 	bootccf_cw <- tsboot(cbind(x1,z1), function(x) ccf(x[,1],x[,2],na.action=na.pass, plot=FALSE, lag.max=LAG.MAX)$acf, R=Rboot, sim="fixed", l=LAG.MAX*2, parallel="snow");
 	ccf_cw <- unlist(apply(bootccf_cw$t, MARGIN=2, function(x) mean(x, na.rm=TRUE)));
 	ccfs_cw <- unlist(zoo::na.locf(zoo::na.locf(rollapply(ccf_cw, width=13, FUN="mean", fill=NA) , na.rm=FALSE), fromLast=TRUE));
 	ccf_cw1 <- unlist(apply(abs(bootccf_cw$t), MARGIN=1, function(x) which.max(x)));
 	ccfs_cw1 <- unlist(apply(bootccf_cw$t, MARGIN=1, function(x) which.max(abs(as.vector(zoo::na.locf(zoo::na.locf(rollapply(x, width=13, FUN="mean", fill=NA), na.rm=FALSE), fromLast=TRUE))))));
	hdi_cw <- HDInterval::hdi(ccf_cw1, credMass=cM);
 	hdis_cw <- HDInterval::hdi(ccfs_cw1, credMass=cM);
 	bootccf_tc <- tsboot(cbind(x2,y2), function(x) ccf(x[,1],x[,2],na.action=na.pass, plot=FALSE, lag.max=LAG.MAX)$acf, R=Rboot, sim="fixed", l=LAG.MAX*2, parallel="snow");
 	ccf_tc <- unlist(apply(bootccf_tc$t, MARGIN=2, function(x) mean(x, na.rm=TRUE)));
	ccfs_tc <- unlist(zoo::na.locf(zoo::na.locf(rollapply(ccf_tc, width=13, FUN="mean", fill=NA) , na.rm=FALSE), fromLast=TRUE));
	ccf_tc1 <- unlist(apply(abs(bootccf_tc$t), MARGIN=1, function(x) which.max(x)));
 	ccfs_tc1 <- unlist(apply(bootccf_tc$t, MARGIN=1, function(x) which.max(abs(as.vector(zoo::na.locf(zoo::na.locf(rollapply(x, width=13, FUN="mean", fill=NA), na.rm=FALSE), fromLast=TRUE))))));
 	hdi_tc <- HDInterval::hdi(ccf_tc1, credMass=cM);
	hdis_tc <- HDInterval::hdi(ccfs_tc1, credMass=cM);
	bootccf_wc <- tsboot(cbind(x3,z3), function(x) ccf(x[,1],x[,2],na.action=na.pass, plot=FALSE, lag.max=LAG.MAX)$acf, R=Rboot, sim="fixed", l=LAG.MAX*2, parallel="snow");
 	ccf_wc <- unlist(apply(bootccf_wc$t, MARGIN=2, function(x) mean(x, na.rm=TRUE)));
    	ccfs_wc <- unlist(zoo::na.locf(zoo::na.locf(rollapply(ccf_wc, width=13, FUN="mean", fill=NA) , na.rm=FALSE), fromLast=TRUE))
	ccf_wc1 <- unlist(apply(abs(bootccf_wc$t), MARGIN=1, function(x) which.max(x)));
 	ccfs_wc1 <- unlist(apply(bootccf_wc$t, MARGIN=1, function(x) which.max(abs(as.vector(zoo::na.locf(zoo::na.locf(rollapply(x, width=13, FUN="mean", fill=NA), na.rm=FALSE), fromLast=TRUE))))));
 	hdi_wc <- HDInterval::hdi(ccf_wc1, credMass=cM);
 	hdis_wc <- HDInterval::hdi(ccfs_wc1, credMass=cM);
	map <- round(c(bayestestR::map_estimate(abs(ccf_ct1+rnorm(length(ccf_ct1),0,0.0001))), bayestestR::map_estimate(abs(ccf_cw1+rnorm(length(ccf_cw1),0,0.0001))), bayestestR::map_estimate(abs(ccf_tc1+rnorm(length(ccf_tc1),0,0.0001))), bayestestR::map_estimate(abs(ccf_wc1+rnorm(length(ccf_wc1),0,0.0001)))),0);
	maps <- round(c(bayestestR::map_estimate(abs(ccfs_ct1+rnorm(length(ccf_ct1),0,0.0001))), bayestestR::map_estimate(abs(ccfs_cw1+rnorm(length(ccf_cw1),0,0.0001))), bayestestR::map_estimate(abs(ccfs_tc1+rnorm(length(ccf_tc1),0,0.0001))), bayestestR::map_estimate(abs(ccfs_wc1+rnorm(length(ccf_wc1),0,0.0001)))),0)
	}     
 
## STATS
	mc_peak_w <- which.max(abs(cross_cov_w))
	mc_peak_ts <- which.max(abs(cross_cov_ts))
 	corr_est <- c(ccf_ct[map[1]], ccf_cw[map[2]], ccf_tc[map[3]], ccf_wc[map[4]])
	corr_est_s <- c(ccfs_ct[maps[1]], ccfs_cw[maps[2]], ccfs_tc[maps[3]], ccfs_wc[maps[4]]) 
	corr_ind <- which.max(abs(corr_est_s))
    	corr_max <- corr_est_s[corr_ind]
	peak_ref <- maps[corr_ind]
	if (corr_ind==1) hdis <- hdis_ct
	if (corr_ind==2) hdis <- hdis_cw
	if (corr_ind==3) hdis <- hdis_tc
	if (corr_ind==4) hdis <- hdis_wc
    
 	cov_est_w <- cross_cov_w[mc_peak_w]
	cov_est_ts <- cross_cov_w[mc_peak_ts]
	cov_est <- c(cross_cov_w[map[1]], cross_cov_w[map[2]], cross_cov_w[map[3]], cross_cov_w[map[4]])
	cov_est_s <-  c(cross_cov_w[maps[1]], cross_cov_w[maps[2]], cross_cov_w[maps[3]], cross_cov_w[maps[4]])

## PLOT
    if (plot.it) {
	par(mfrow=c(3,2), mar=c(5,4,2,1), oma=c(1,1,5,0.5), las=0, cex.axis=1.3, cex.lab=1.3)
	plot((-LAG.MAX:LAG.MAX), cross_cov_w, ylab="cross-cov (c,w)", xlab="Lag (sec)", type="h", col="grey68", ylim=c(min(cross_cov_w*1.05,0),max(cross_cov_w*1.05,0)), xaxt="n")
	axis(side=1, at=seq(-LAG.MAX,LAG.MAX,40), labels=seq(-LAG.MAX,LAG.MAX,40)/mfreq)  
 	polygon(x=c(hdis[1]:hdis[2], hdis[2]:hdis[1])-LAG.MAX-1, y=c(cross_cov_w[hdis[1]:hdis[2]],rep(0,hdis[2]-hdis[1]+1)), col="lightblue", border="lightblue")
	segments(x0=peak_ref - LAG.MAX - 1, y0=0, y1=cross_cov_w[peak_ref], col=2, lwd=2)
 	segments(x0=mc_peak_w - LAG.MAX - 1, y0=0, y1=cross_cov_w[mc_peak_w], col=1, lwd=1)
	mtext(side=3, line=.5, adj=0, paste0("Peak at ", (mc_peak_w-LAG.MAX-1)/mfreq, " sec"), cex=1.1) 
	mtext(side=3, line=.5, adj=1, "a", cex=1.5, font=2) 
	box(lwd=1.5)
	mtext(side=3, line=3, paste0("Time lag (MAP) at: ", (peak_ref-LAG.MAX-1)/20, " sec"), cex=1.25, col=2)
	plot((-LAG.MAX:LAG.MAX), cross_cov_ts, ylab="cross-cov (c,t)", xlab="Lag (sec)", type="h", col="grey68", ylim=c(min(cross_cov_ts*1.05,0),max(cross_cov_ts*1.05,0)), xaxt="n")    
	axis(side=1, at=seq(-LAG.MAX,LAG.MAX,40), labels=seq(-LAG.MAX,LAG.MAX,40)/mfreq)  
	polygon(x=c(hdis[1]:hdis[2], hdis[2]:hdis[1])-LAG.MAX-1, y=c(cross_cov_ts[hdis[1]:hdis[2]],rep(0,hdis[2]-hdis[1]+1)), col="lightblue", border="lightblue")
	segments(x0=peak_ref - LAG.MAX - 1, y0=0, y1=cross_cov_ts[peak_ref], col=2, lwd=2)
	segments(x0=mc_peak_ts - LAG.MAX - 1, y0=0, y1=cross_cov_ts[mc_peak_ts], col=1, lwd=1)
	mtext(side=3, line=.5, adj=0, paste0("Peak at ", (mc_peak_ts-LAG.MAX-1)/mfreq, " sec"), cex=1.1) 
 	mtext(side=3, line=.5, adj=1, "b", cex=1.5, font=2) 
	box(lwd=1.5)
	mtext(side=3, line=3, paste0("Unc (",cM*100, "% HDI): [", (hdis[1]-LAG.MAX-1)/20, "; ", (hdis[2]-LAG.MAX-1)/20, "] sec"), cex=1.25, col=2)

	plot((-LAG.MAX:LAG.MAX), ccf_cw, ylab="pwb cross-cor (c,w)", xlab="Lag (sec)", type="h", col="grey68", ylim=c(min(ccf_cw,-4/sqrt(length(x))),max(ccf_cw,4/sqrt(length(x)))), xaxt="n")    
    	lines((-LAG.MAX:LAG.MAX), ccfs_cw, col=1, lwd=2)
	axis(side=1, at=seq(-LAG.MAX,LAG.MAX,40), labels=seq(-LAG.MAX,LAG.MAX,40)/mfreq)  
    	abline(h=c(-3.291,3.291)/sqrt(length(x)*13), col=4, lty=2, lwd=2)
    	points(x=maps[2] - LAG.MAX - 1, y=min(ccf_cw,-4/sqrt(length(x))), pch=24, col=1, cex=1.25, bg="red")
    	mtext(side=3, line=1, adj=1, "c", cex=1.5, font=2) 
   	mtext(side=3, line=.5, adj=0, paste0("Peak at ", (maps[2]-LAG.MAX-1)/mfreq ," sec"), cex=1.1) 
	box(lwd=1.5)

	plot((-LAG.MAX:LAG.MAX), ccf_ct, ylab="pwb cross-cor (c,t)", xlab="Lag (sec)", type="h", col="grey68", ylim=c(min(ccf_ct,-4/sqrt(length(x))),max(ccf_ct,4/sqrt(length(x)))), xaxt="n")    
    	lines((-LAG.MAX:LAG.MAX), ccfs_ct, col=1, lwd=2)
	axis(side=1, at=seq(-LAG.MAX,LAG.MAX,40), labels=seq(-LAG.MAX,LAG.MAX,40)/mfreq)  
    	abline(h=c(-3.291,3.291)/sqrt(length(x)*13), col=4, lty=2, lwd=2)
    	points(x=maps[1] - LAG.MAX - 1, y=min(ccf_ct,-4/sqrt(length(x))), pch=24, col=1, cex=1.25, bg="red")
	mtext(side=3, line=1, adj=1, "d", cex=1.5, font=2) 
   	mtext(side=3, line=.5, adj=0, paste0("Peak at ", (maps[1]-LAG.MAX-1)/mfreq ," sec"), cex=1.1) 
	box(lwd=1.5)

	plot((-LAG.MAX:LAG.MAX), ccf_wc, ylab="pwb cross-cor (w,c)", xlab="Lag (sec)", type="h", col="grey68", ylim=c(min(ccf_wc,-4/sqrt(length(x))),max(ccf_wc,4/sqrt(length(x)))), xaxt="n")    
    	lines((-LAG.MAX:LAG.MAX), ccfs_wc, col=1, lwd=2)
	axis(side=1, at=seq(-LAG.MAX,LAG.MAX,40), labels=seq(-LAG.MAX,LAG.MAX,40)/mfreq)  
    	abline(h=c(-3.291,3.291)/sqrt(length(x)*13), col=4, lty=2, lwd=2)
    	points(x=maps[4] - LAG.MAX - 1, y=min(ccf_wc,-4/sqrt(length(x))), pch=24, col=1, cex=1.25, bg="red")
	mtext(side=3, line=1, adj=1, "e", cex=1.5, font=2) 
   	mtext(side=3, line=.5, adj=0, paste0("Peak at ", (maps[4]-LAG.MAX-1)/mfreq ," sec"), cex=1.1) 
	box(lwd=1.5)

 	plot((-LAG.MAX:LAG.MAX), ccf_tc, ylab="pwb cross-cor (t,c)", xlab="Lag (sec)", type="h", col="grey68", ylim=c(min(ccf_tc,-4/sqrt(length(x))),max(ccf_tc,4/sqrt(length(x)))), xaxt="n")    
   	lines((-LAG.MAX:LAG.MAX), ccfs_tc, col=1, lwd=2)
	axis(side=1, at=seq(-LAG.MAX,LAG.MAX,40), labels=seq(-LAG.MAX,LAG.MAX,40)/mfreq)  
    	abline(h=c(-3.291,3.291)/sqrt(length(x)*13), col=4, lty=2, lwd=2)
    	points(x=maps[3] - LAG.MAX - 1, y=min(ccf_tc,-4/sqrt(length(x))), pch=24, col=1, cex=1.25, bg="red")
    	mtext(side=3, line=1, adj=1, "f", cex=1.5, font=2) 
   	mtext(side=3, line=.5, adj=0, paste0("Peak at ", (maps[3]-LAG.MAX-1)/mfreq ," sec"), cex=1.1) 
	box(lwd=1.5)
   }

return(list("tlag_ref"=peak_ref, "tlag_lower"= hdis[1],"tlag_upper"= hdis[2], "corr_ref"=corr_max, 
	    "mc_tlag_w"=mc_peak_w, "mc_tlag_ts"=mc_peak_ts, "corr_est" = corr_est, "corr_est_s" = corr_est_s,
	    "cov_est_w"= cov_est_w,"cov_est_ts"= cov_est_ts,"cov_est"= cov_est, "cov_est_s"= cov_est_s,
	    "kid_w"=kid_w, "kid_ts"=kid_ts, "kid_c"=kid_c,
	    "hdi"=c(hdi_ct, hdi_cw, hdi_tc, hdi_wc,hdis_ct, hdis_cw, hdis_tc, hdis_wc), "map"=map, "maps"=maps))
}
