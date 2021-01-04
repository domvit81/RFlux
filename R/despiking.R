
despiking <- function(x, mfreq, variant, wsignal, wscale, wby=1, zth=5, alpha=0.01){
	## x vector
	## mfreq the main frequency (24, 48 for hourly and half-hourly time series respectively; 10, 20 for EC raw data acquired at 10Hz or 20Hz, respectively)
	## variant flx mainly designed for (half-)hourly EC fluxes (missing values are allowed); met mainly designed for (half-)hourly meteo time series (a low percentage of missing values is allowed); raw mainly designed for high-frequency EC raw data
	## wsignal the window widht used to estimate the underlying signal (only for met and raw variants)
	## wscale the window widht used to estimate the local scale paramater (only for met and raw variants)
	## wby calculate the scale parameter at every by-th point rather than every point (default=1). Large values of wby reduce the computational time but can introduce bias in the scale parameter estimation.
	## alpha the significance level used in the outlier detection rule (only for flx variant)
	## zth  the threshold value of the z-sigma rule ((only for met and raw variants)
	N <- length(x)
	na_index <- which(is.na(x))

	if(variant=="v1"){ ## mainly for half-hourly EC fluxes
		TS <- ts(x, frequency=mfreq)
		if(N < mfreq*10) {warning(call.=FALSE, "Outlier detection is performed when data cover a period of at least 10 consecutive days")}
		if(any(apply(matrix(TS, nrow=mfreq), MARGIN=1, function(x) sum(is.na(x))==(N/mfreq-3)))) {warning(call.=FALSE, "Too missing values for some sub-series! Outlier detection cannot be performed with STL algorithm")}
		if(N >= mfreq*10 & all(apply(matrix(TS, nrow=mfreq), MARGIN=1, function(x) sum(is.na(x))<(N/mfreq-3)))){
			spike_index <- c()
			if(length(TS>(mfreq*10-1)) & all(apply(matrix(TS, nrow=mfreq), MARGIN=1, function(x) length(which(is.na(x)))) < ceiling(N/mfreq))){
				mod <- stlplus(log10(ts(TS, start=c(1,1), frequency=mfreq)+1000), s.window=7, t.window=mfreq*7+1, s.degree=0, fc.degree=c(1,0), fc.window=c(mfreq*7+1,7) , outer=20, robust=TRUE, fc.name=c("Long-Run","Short-Run"));
				signal <- as.vector(10^(apply(cbind(mod$data[,2], mod$fc[,1:2]), MARGIN=1, FUN=sum))-1000);
				res <- as.vector(TS - signal);
				alpha <- alpha;
				indq <- quantile(signal, seq(0,1,.1), na.rm=TRUE);
				for (i in 1:(length(indq)-1)){
					ind <- which(signal >= indq[i] & signal <= indq[i+1])
					if(length(na.omit(res[ind]))>9) {asd <- aout.laplace(res[ind], c(median(res[ind], na.rm=TRUE), Qn(na.omit(res[ind]))), alpha); spike_index <- c(spike_index, ind[which(asd$is.outlier)])}
					if(length(na.omit(res[ind]))<10) {spike_index <- c(spike_index, ind[which(!is.na(res[ind]))])}## to prevent cases where missing values unable the application of aaut.laplace
				}
			}
		}
	}
	
	if(variant=="v2"){ ## mainly for meteo time series
		TS <- na.locf(na.locf(x, na.rm=FALSE), fromLast=TRUE)
		spike_tmp <- c()
		TSd <- c(NA,diff(TS))
		th_local <- as.numeric(rollapply(as.ts(TSd), width=wscale, by=wby, function(x) zth*Qn(na.omit(x)), fill=NA))
		spike_tmp <- which(abs(TSd) > th_local)
				
		spike_index <- c()
		estimate <- c()
		if (length(spike_tmp)>0) {
			for (i in 1:length(spike_tmp)){
				if (spike_tmp[i] < (wsignal-1)/2 | spike_tmp[i] >= N-(wsignal-1)/2) {spike <- c(spike, spike_tmp[i])}
				if (spike_tmp[i] > (wsignal-1)/2 & spike_tmp[i] < N-(wsignal-1)/2) {
					ind_window_start <- spike_tmp[i]-(wsignal-1)/2;
					ind_window_end <- spike_tmp[i]+(wsignal-1)/2;
					TSw <- TS[ind_window_start:ind_window_end];
					signal <- hybrid.filter(TSw, width=wsignal, method=c("MED", "RM", "MH", "MMH", "PRMH", "PRMMH"));		
					res0 <- as.vector(TSw-signal$level)[(wsignal+1)/2,];
					ind <- which.min(abs(res0));
					res <- res0[ind];
					if(abs(res) > th_local[spike_tmp[i]]) {spike_index <- c(spike_index, spike_tmp[i])}
					}
				}
			}
		}

	if(variant=="v3"){ ## mainly for EC raw data
		TS <- na.locf(na.locf(x, na.rm=FALSE), fromLast=TRUE)
		x_raw <- TS + rnorm(N, 0, abs(0.0001*TS))
		
		mod0 <- rm.filter(x_raw, width=wsignal)
		
		ifelse(length(which(is.na(mod0$level$RM))) > 1, 
			{res0 <- c(NA, diff(x_raw));
			lbound <- Qn(res0);
			v00 <- na.locf(na.locf(na.approx(rollapply(res0, width=wscale, by=wby, function(x) Qn(x), fill=NA), na.rm=FALSE), na.rm=FALSE), fromLast=TRUE);
			v0 <- replace(v00, which(v00 < lbound), lbound);
			spike_index <- which(abs(res0) > zth*v0);
			n_spike <- length(spike_index);
			x_cl <- na.locf(na.locf(replace(x_raw, spike_index, NA), na.rm=FALSE), fromLast=TRUE)
			},
			{res0 <- TS - mod0$level$RM;
			lbound <- Qn(res0);
			v00 <- na.locf(na.locf(na.approx(rollapply(res0, width=wscale, by=wby, function(x) Qn(x), fill=NA), na.rm=FALSE), na.rm=FALSE), fromLast=TRUE);
			v0 <- replace(v00, which(v00 < lbound), lbound)
			spike_index <- which(abs(res0) > zth*v0);
			n_spike <- length(spike_index);
			x_cl <- na.locf(na.locf(replace(TS, spike_index, NA), na.rm=FALSE), fromLast=TRUE)
			}) ## per preservare sudden shift in media non usiamo approssimazioni lineari
	}
	return(list("despiked_ts"=replace(x, spike_index, NA), "spike_loc"=spike_index))
}
