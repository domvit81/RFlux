


parallel_despiking <- function(x, mfreq, wsignal = mfreq*30+1, wscale = mfreq*30+1, wby=1, zth=5){
	mm <- matrix(x, nrow=mfreq*60*5)
	dspk <- future_apply(mm, MARGIN=2, function(x) despiking(x, variant="v3", wsignal=wsignal, wscale=wscale, wby=wby, zth=zth))
	spike_index <- as.vector(c(unlist(dspk[[1]]$spike_loc), 
							   unlist(dspk[[2]]$spike_loc)+mfreq*60*5,
							   unlist(dspk[[3]]$spike_loc)+mfreq*60*5*2,
							   unlist(dspk[[4]]$spike_loc)+mfreq*60*5*3,
							   unlist(dspk[[5]]$spike_loc)+mfreq*60*5*4,
							   unlist(dspk[[6]]$spike_loc)+mfreq*60*5*5))
	return(list("despiked_ts"=replace(x, spike_index, NA), "spike_loc"=spike_index))
	}
