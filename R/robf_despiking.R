


robf_despiking <- function(x, mfreq, file_length){
	# Window width selection
	nL <- mfreq*60*file_length
	x <- x[1:nL]
	mod_rlm <- rlm(as.vector(x)~poly(seq(1,nL,1), 5, raw = FALSE))
	ind <- which(abs(as.vector(residuals(mod_rlm))) > 3*Qn(as.vector(residuals(mod_rlm))))
	out_patch1 <- statsNA(replace(x, ind, NA), print_only=FALSE)$longest_na_gap
	out_patch2 <- max(rollapply(replace(x, ind, NA), width=mfreq*30+1, by=1, function(x) length(which(is.na(x)))))
	wl0 <- max(c(mfreq*5+1, out_patch1*4, out_patch2*4), na.rm=TRUE)
	ifelse(odd(wl0), wl <- wl0, wl <- wl0+1)
	
	if(wl > mfreq*60+1){
		dspk <- despiking(x, mfreq=mfreq, variant="v3", wsignal=wl, wscale=wl, wby=1, zth=5);
		spike_loc <- dspk$spike_loc;
		ts_cleaned <- dspk$ts_cleaned
		}
	
	if(wl <= mfreq*60+1){
		mm <- matrix(x, nrow=nL/6);
		plan(multicore);
		dspk <- future_apply(mm, MARGIN=2, function(x) despiking(x, mfreq=mfreq, variant="v3", wsignal=wl, wscale=wl, wby=1, zth=5));
		spike_loc <- na.omit(as.vector(c(dspk[[1]]$spike_loc, dspk[[2]]$spike_loc + mfreq*60*5, dspk[[3]]$spike_loc + mfreq*60*5*2, dspk[[4]]$spike_loc + mfreq*60*5*3, dspk[[5]]$spike_loc + mfreq*60*5*4, dspk[[6]]$spike_loc + mfreq*60*5*5)));
		ts_cleaned <- as.vector(c(dspk[[1]]$ts_cleaned, dspk[[2]]$ts_cleaned,dspk[[3]]$ts_cleaned, dspk[[4]]$ts_cleaned, dspk[[5]]$ts_cleaned, dspk[[6]]$ts_cleaned))
		}
	
	return(list("ts_cleaned"=ts_cleaned, "spike_loc"=spike_loc, "wl"=wl))
}

